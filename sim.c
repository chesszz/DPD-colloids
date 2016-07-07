#include "utils.h"
#include "sim.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define SKIN 1.0
#define R_SS 1.0

void refold_positions(Dyn_Vars *dyn_vars, Inputs in, double *pos_list, 
                    double *vel_list, int num_obj);

void calculate_acc(Dyn_Vars *dyn_vars, Inputs in, int *neigh_list_p, 
                    int *neigh_list_w, int *neigh_list_pw, double strength);

void get_rel_vector(Inputs in, int t, double *pos_list_1, double *pos_list_2, 
                    int i, int j, double *Rij);

void update_neigh_list(Inputs in, int t, double *pos_list, double cutoff_dist, 
                    int *neigh_list, int num_obj);

void update_neigh_list_pw(Inputs in, int t, 
                        double *pos_list_p, double *pos_list_w, 
                        double cutoff_dist, int *neigh_list_pw, 
                        int num_p, int num_w);

int check_update_req(double *disp_list_w, double *disp_list_p, int n_w, int n_p);
double calc_temp(Dyn_Vars *dyn_vars, Inputs in);


/* Evolves the system from the given initial conditions and input parameters.
 * Will print the variables' time evolution into a file.
 */ 
void evolve_system(Dyn_Vars *dyn_vars, Inputs in) {

    #if DEBUG
        int list_recomp_count = 0;
    #endif

    /* Used to accumulate the displacements each time step. */
    double disp;

    double dt = in.TIME_STEP;
    int num_inc = 26;
    int t_stab = 50;

    int checkpoint = in.N_STEPS / 10;
    FILE *wat_out = fopen("water.out", "w");
    FILE *part_out = fopen("particles.out", "w");
    FILE *temp_out = fopen("temp.out", "w");

    int update_req = 1;

    /* List of displacements for each water since the last update. */
    double *disp_list_w = (double *) calloc(3 * in.N_WATER, sizeof(double));
    /* List of displacements for each particle since the last update. */
    double *disp_list_p = (double *) calloc(3 * in.N_PARTICLES, sizeof(double));

    /* List of water-water neighbours. Boundaries are marked by -1. Note
     * that because we only consider pairs of waters i, j where i < j, thus
     * we have at most N(N-1)/2 pairs. However, we also need ~N dividers, hence
     * we have a total of N(N+1)/2 entries. */
    int *neigh_list_w = (int *) calloc(in.N_WATER * (in.N_WATER + 1) / 2, sizeof(int));
    /* List of particle-particle neighbours. */
    int *neigh_list_p = (int *) calloc(in.N_PARTICLES * (in.N_PARTICLES + 1) / 2, sizeof(int));
    /* List of particle-water neighbours. Particles are the "host". We can have 
     * N_PARTICLES * N_WATER possible combinations, and then add another 
     * N_PARTICLES for the dividers. */
    int *neigh_list_pw = (int *) calloc(in.N_PARTICLES * (in.N_WATER + 1), sizeof(int));

    if (disp_list_w  == NULL || disp_list_p == NULL  || 
        neigh_list_w == NULL || neigh_list_p == NULL || neigh_list_pw == NULL) {
        fprintf(stderr, "Error in allocating memory.\n");
        exit(1);
    }

    fprintf(wat_out, "Time Step, Water Number, Position, Velocity, Acceleration\n");
    fprintf(part_out, "Time Step, Particle Number, Position, Velocity, Acceleration\n");

    /* STABILISING LOOP. */
    for (int inc = 1; inc < num_inc; inc++) {

        /* inc ranges from 1 ~ (num_inc - 1)
         * strength ranges from 1 / num_inc ~ (num_inc - 1) / num_inc */
        double strength = 1.0 / num_inc * inc;

        /* Simulate t_stab time steps for each incremement. */
        for (int t = 0; t < t_stab; t++) {

            (dyn_vars->t)++;

            /* Velocity Verlet 1: V0(t+0.5*dt) */
            for (int i = 0; i < 3*in.N_WATER; i++) {
                dyn_vars->watvel[i] += 0.5 * dt * dyn_vars->watacc[i];
            }
            for (int i = 0; i < 3*in.N_PARTICLES; i++) {
                dyn_vars->partvel[i] += 0.5 * dt * dyn_vars->partacc[i];
            }

            /* Velocity Verlet 2: R(t+dt) */
            for (int i = 0; i < 3*in.N_WATER; i++) {
                /* displacement = v * dt */
                disp = dt * dyn_vars->watvel[i];
                /* Add displacement to both the position and displacement list. */
                dyn_vars->watpos[i] += disp;
                /* x-Displacement gets extra term from the shear velocity. x
                 * coordinate is a multiple of 3, so the logical test gives 1, 
                 * other axes give 0  for the test and so have no shear 
                 * displacement. */
                disp_list_w[i]      += (disp + (i % 3 == 0) * in.V_SHEAR * dt);
            } 
            for (int i = 0; i < 3*in.N_PARTICLES; i++) {
                disp = dt * dyn_vars->partvel[i];
                /* Add displacement to both the position and displacement list. */
                dyn_vars->partpos[i] += disp;
                /* Displacement gets extra term from the shear velocity. */
                disp_list_p[i]       += (disp + (i % 3 == 0) * in.V_SHEAR * dt);
            }

            /* Wrap the objects back to the box if they escaped. Handle both water
             * and particles. */
            refold_positions(dyn_vars, in, dyn_vars->watpos, dyn_vars->watvel, in.N_WATER);
            refold_positions(dyn_vars, in, dyn_vars->partpos, dyn_vars->partvel, in.N_PARTICLES);

            /* Update neighbour list if required. */
            if (update_req) {
                update_neigh_list(in, t, dyn_vars->watpos, R_SS, neigh_list_w, in.N_WATER);
                update_neigh_list(in, t, dyn_vars->partpos, in.R_CC, neigh_list_p, in.N_PARTICLES);
                update_neigh_list_pw(in, t, dyn_vars->partpos, dyn_vars->watpos, 
                                    in.R_SC, neigh_list_pw, 
                                    in.N_PARTICLES, in.N_WATER);

                /* Set the displacement lists to be 0. */
                for (int i = 0; i < 3*in.N_WATER; i++) {
                    disp_list_w[i] = 0;
                }
                for (int i = 0; i < 3*in.N_PARTICLES; i++) {
                    disp_list_p[i] = 0;
                }

                update_req = 0;       
            }            

            /* Velocity Verlet 3: F(t+dt) */
            calculate_acc(dyn_vars, in, neigh_list_p, neigh_list_w, neigh_list_pw, strength);

            /* Velocity Verlet 4: V(t+dt) */
            for (int i = 0; i < 3*in.N_WATER; i++) {
                dyn_vars->watvel[i] += 0.5 * dt * dyn_vars->watacc[i];
            }
            for (int i = 0; i < 3*in.N_PARTICLES; i++) {
                dyn_vars->partvel[i] += 0.5 * dt * dyn_vars->partacc[i];
            }

            /* Check if the objects have moved too much since the last update
             * that requires a new neighbour list to be created.
             */
            update_req = check_update_req(disp_list_w, disp_list_p, in.N_WATER, in.N_PARTICLES);
        }

        fprintf(stderr, "Finished increment loop %d. Time is now %d.\n", inc, dyn_vars->t);
    }
    
    fprintf(stderr, "Finished all increment loops. Time is now %d.\n", dyn_vars->t);
    int t_after_settle = dyn_vars->t;

    /* -------------------------------------------------------------------- */

    /* MAIN TIME LOOP */
    for (int t = t_after_settle; t < t_after_settle + in.N_STEPS; t++) {

        (dyn_vars->t)++;

        /* Velocity Verlet 1: V0(t+0.5*dt) */
        for (int i = 0; i < 3*in.N_WATER; i++) {
            dyn_vars->watvel[i] += 0.5 * dt * dyn_vars->watacc[i];
        }
        for (int i = 0; i < 3*in.N_PARTICLES; i++) {
            dyn_vars->partvel[i] += 0.5 * dt * dyn_vars->partacc[i];
        }

        /* Velocity Verlet 2: R(t+dt) */
        for (int i = 0; i < 3*in.N_WATER; i++) {
            /* displacement = v * dt */
            disp = dt * dyn_vars->watvel[i];
            /* Add displacement to both the position and displacement list. */
            dyn_vars->watpos[i] += disp;
            /* x-Displacement gets extra term from the shear velocity. x
             * coordinate is a multiple of 3, so the logical test gives 1, 
             * other axes give 0  for the test and so have no shear 
             * displacement. */
            disp_list_w[i]      += (disp + (i % 3 == 0) * in.V_SHEAR * dt);
        } 
        for (int i = 0; i < 3*in.N_PARTICLES; i++) {
            disp = dt * dyn_vars->partvel[i];
            /* Add displacement to both the position and displacement list. */
            dyn_vars->partpos[i] += disp;
            /* Displacement gets extra term from the shear velocity. */
            disp_list_p[i]       += (disp + (i % 3 == 0) * in.V_SHEAR * dt);
        }

        /* Wrap the objects back to the box if they escaped. Handle both water
         * and particles. */
        refold_positions(dyn_vars, in, dyn_vars->watpos, dyn_vars->watvel, in.N_WATER);
        refold_positions(dyn_vars, in, dyn_vars->partpos, dyn_vars->partvel, in.N_PARTICLES);

        /* Update neighbour list if required. */
        if (update_req) {
            update_neigh_list(in, t, dyn_vars->watpos, R_SS, neigh_list_w, in.N_WATER);
            update_neigh_list(in, t, dyn_vars->partpos, in.R_CC, neigh_list_p, in.N_PARTICLES);
            update_neigh_list_pw(in, t, dyn_vars->partpos, dyn_vars->watpos, 
                                in.R_SC, neigh_list_pw, 
                                in.N_PARTICLES, in.N_WATER);

            #if DEBUG
                fprintf(stderr, "List recomputed: %d\n", list_recomp_count);
                list_recomp_count++;
            #endif       

            /* Set the displacement lists to be 0. */
            for (int i = 0; i < 3*in.N_WATER; i++) {
                disp_list_w[i] = 0;
            }
            for (int i = 0; i < 3*in.N_PARTICLES; i++) {
                disp_list_p[i] = 0;
            }

            update_req = 0;       
        }

        /* Velocity Verlet 3: F(t+dt) */
        calculate_acc(dyn_vars, in, neigh_list_p, neigh_list_w, neigh_list_pw, 1.0);

        /* Velocity Verlet 4: V(t+dt) */
        for (int i = 0; i < 3*in.N_WATER; i++) {
            dyn_vars->watvel[i] += 0.5 * dt * dyn_vars->watacc[i];
        }
        for (int i = 0; i < 3*in.N_PARTICLES; i++) {
            dyn_vars->partvel[i] += 0.5 * dt * dyn_vars->partacc[i];
        }

        /* Check if the objects have moved too much since the last update
         * that requires a new neighbour list to be created.
         */
        update_req = check_update_req(disp_list_w, disp_list_p, in.N_WATER, in.N_PARTICLES);

        /* Print out the positions and velocities of everything for each t. */
        #if VERBOSE
            for (int i = 0; i < in.N_WATER; i++) {
                fprintf(stderr, "Time: %d, Water    %d, position: (%f, %f, %f), velocity: (%f, %f, %f), acc: (%f, %f, %f)\n", 
                    t, i, 
                    dyn_vars->watpos[3*i], dyn_vars->watpos[3*i+1], dyn_vars->watpos[3*i+2], 
                    dyn_vars->watvel[3*i], dyn_vars->watvel[3*i+1], dyn_vars->watvel[3*i+2],
                    dyn_vars->watacc[3*i], dyn_vars->watacc[3*i+1], dyn_vars->watacc[3*i+2]);
            }
            for (int i = 0; i < in.N_PARTICLES; i++) {
                fprintf(stderr, "Time: %d, Particle %d, position: (%f, %f, %f), velocity: (%f, %f, %f), acc: (%f, %f, %f)\n", 
                    t, i, 
                    dyn_vars->partpos[3*i], dyn_vars->partpos[3*i+1], dyn_vars->partpos[3*i+2], 
                    dyn_vars->partvel[3*i], dyn_vars->partvel[3*i+1], dyn_vars->partvel[3*i+2],
                    dyn_vars->partacc[3*i], dyn_vars->partacc[3*i+1], dyn_vars->partacc[3*i+2]);
            }
        #endif

        /* At each time step, print the positions, velocities, and 
         * accelerations for each particle.
         */
        for (int i = 0; i < in.N_WATER; i++) {
            fprintf(wat_out, "%d, %d, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", 
                t, i,
                dyn_vars->watpos[3*i], dyn_vars->watpos[3*i+1], dyn_vars->watpos[3*i+2], 
                dyn_vars->watvel[3*i], dyn_vars->watvel[3*i+1], dyn_vars->watvel[3*i+2],
                dyn_vars->watacc[3*i], dyn_vars->watacc[3*i+1], dyn_vars->watacc[3*i+2]);
        }
        for (int i = 0; i < in.N_PARTICLES; i++) {
            fprintf(part_out, "%d, %d, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", 
                t, i, 
                dyn_vars->partpos[3*i], dyn_vars->partpos[3*i+1], dyn_vars->partpos[3*i+2], 
                dyn_vars->partvel[3*i], dyn_vars->partvel[3*i+1], dyn_vars->partvel[3*i+2],
                dyn_vars->partacc[3*i], dyn_vars->partacc[3*i+1], dyn_vars->partacc[3*i+2]);
        }

        /* At each time step, print the temperature. */
        fprintf(temp_out, "%d, %f\n", t, calc_temp(dyn_vars, in));

        /* Print a message every 10% complete. */
        if (checkpoint != 0 && t % checkpoint == 0) {
            fprintf(stderr, "Time step %d / %d complete.\n", t, t_after_settle + in.N_STEPS);
        }  
    }

    fclose(wat_out);
    fclose(part_out);
    fclose(temp_out);

    free(disp_list_w);
    free(disp_list_p);
    free(neigh_list_w);
    free(neigh_list_p);
}

/* Enforce periodic BCs. */
void refold_positions(Dyn_Vars *dyn_vars, Inputs in, double *pos_list, double *vel_list, int num_obj) {

    double size = in.BOX_SIZE;
    int t = dyn_vars->t;

    /* If water escapes from the left or right side of the box, wrap it back
     * to the inside. Each particle is handled together.
     */
    for (int i = 0; i < num_obj; i++) {

        /* Indices for the i-th object. */
        int x = 3 * i;
        int y = x + 1;
        int z = y + 1;

        /* Wrap x dimension. */
        if (pos_list[x] > size) {
            pos_list[x] -= size;
        }
        else if (pos_list[x] < 0) {
            pos_list[x] += size;
        }

        /* Wrap y dimension. Implemented L-E BC here. Particles gain -V_SHEAR
         * when they pass through the top and emerge from the bottom, and they
         * gain +V_SHEAR when they pass through the bottom and emerge from the 
         * top. Moreover, passing through the top plate pushes the x-position
         * in the negative direction to account for the bottom cell moving left,
         * and vice versa for the bottom boundary. 
         */
        if (pos_list[y] > size) {
            pos_list[y] -= size;
            pos_list[x] -= t * in.V_SHEAR * in.TIME_STEP - 
                floor(t * in.V_SHEAR * in.TIME_STEP / size) * size;

            vel_list[x] -= in.V_SHEAR;
        }
        else if (pos_list[y] < 0) {
            pos_list[y] += size;
            pos_list[x] += t * in.V_SHEAR * in.TIME_STEP - 
                floor(t * in.V_SHEAR * in.TIME_STEP / size) * size;

            vel_list[x] += in.V_SHEAR;
        }

        /* Wrap z dimension. */
        if (pos_list[z] > size) {
            pos_list[z] -= size;
        }
        else if (pos_list[z] < 0) {
            pos_list[z] += size;
        }
    }
}

/* Calculate the forces acting on each water and particle. */
void calculate_acc(Dyn_Vars *dyn_vars, Inputs in, int *neigh_list_p, int *neigh_list_w, int *neigh_list_pw, double strength) {

    /* Zero out all accelerations. */
    for (int i = 0; i < 3*in.N_WATER; i++) {
        dyn_vars->watacc[i] = 0;
    }
    for (int i = 0; i < 3*in.N_PARTICLES; i++) {
        dyn_vars->partacc[i] = 0;
    }

    /* Calculate sqrt of dt for use in F_r. Calculate it just once here. */
    double sqrt_dt = sqrt(in.TIME_STEP);
    /* Also calculate sigma, the Brownian noise strength. */
    double sigma =  sqrt(2 * in.DAMP_CONST);

    double R_CC_sq = in.R_CC * in.R_CC;
    double R_SC_sq = in.R_SC * in.R_SC;

    /* -------------------------------------------------------------------- */
    /* 1. Water-Water interactions. */
    /* -------------------------------------------------------------------- */

    int i = 0; /* "host" object counter */

    /* Traverse through the array until the "host" object has reached the 
     * last object.
     */
    for (int neigh_counter = 0; i < in.N_WATER; neigh_counter++) {

        int j = neigh_list_w[neigh_counter];

        /* Check that we don't exceed the array. */
        assert(neigh_counter < in.N_WATER * (in.N_WATER + 1) / 2);
        /* We should never reach the -255 part because we should break out using 
         * the i < in.WATER condition by then. 
         */
        assert (j != -255);
        /* Each entry is either a boundary marker (-1) or is an index of a 
         * water, and so lies between 0 and N_WATER-1.
         */
        assert (j == -1 || (j >= 0 && j < in.N_WATER)); 

        #if DEBUG
            printf("counter = %d, i = %d, j = %d\n", neigh_counter, i, j);
        #endif

        /* A neighbour of -1 means that we have reached the end of this host. */
        if (j == -1) {
            i++;
            continue;
        }

        /* Rij stores the relative displacement vector. Points from j to i. */
        double Rij[3];
        get_rel_vector(in, dyn_vars->t, dyn_vars->watpos, dyn_vars->watpos, i, j, Rij);

        /* Find the squared distance between particles. */
        double dist_sq = Rij[0] * Rij[0] + Rij[1] * Rij[1] + Rij[2] * Rij[2];

        #if VERBOSE
            printf("Waters (%d, %d), vector from j to i: (%f,%f,%f), dist_sq: %f\n", 
                i, j, Rij[0], Rij[1], Rij[2], dist_sq);
        #endif
        
        /* Cutoff distance is set to be 1 by definition of the units. */
        if (dist_sq < 1) {

            /* Calculate the relative distance and the factor that models
             * how force decays (linearly) with distance. 
             */
            double dist = sqrt(dist_sq);
            double dist_weight = (1 - dist);

            assert(dist_weight > 0);
            assert(dist_weight <= 1);

            /* Normalise the separation vector. */
            Rij[0] /= dist;
            Rij[1] /= dist;
            Rij[2] /= dist;

            /* Find the relative velocity, Vi - Vj. It is the velocity of i
             * relative to / from the frame of j. 
             */
            double Vij[3];
            Vij[0] = dyn_vars->watvel[3*i]   - dyn_vars->watvel[3*j];
            Vij[1] = dyn_vars->watvel[3*i+1] - dyn_vars->watvel[3*j+1];
            Vij[2] = dyn_vars->watvel[3*i+2] - dyn_vars->watvel[3*j+2];

            double V_dot_R =    Vij[0] * Rij[0] + 
                                Vij[1] * Rij[1] + 
                                Vij[2] * Rij[2];

            #if VERBOSE
                printf("Vi_x: %f, Vj_x: %f\n", dyn_vars->watvel[3*i], dyn_vars->watvel[3*j]);
                printf("Vij_x %f\n", Vij[0]);
                printf("Mod Vij: %f\n", Vij[0] * Vij[0] + Vij[1] * Vij[1] + Vij[2] * Vij[2]);
            #endif  

            /* Conservative Force. */
            double F_c =  in.SPRING_CONST * dist_weight;
            /* Dissipative/Damping Force. */
            double F_d = - in.DAMP_CONST * dist_weight * dist_weight * V_dot_R;
            /* Random/Brownian Force.
             * Variance of a uniform distribution is 1/12 * interval^2.
             * Therefore, for mean = 0 and variance = 1, assume we have -x/2
             * to x/2, then 1/12 * x^2 = 1 ==> x = sqrt(12) ~ 3.464. Note 
             * that we first scale from 0 ~ RAND_MAX ==> 0 ~ 1 ==> 
             * -1/2 ~ 1/2 ==> -sqrt(12)/2 ~ sqrt(12)/2.
             */
            double F_r = sigma * dist_weight * 3.464 * ((rand()/(double) RAND_MAX) - 0.5) / sqrt_dt;

            #if VERBOSE
                printf("Waters (%d, %d), F_c: %f, F_d: %f, F_r: %f\n", i, j, F_c, F_d, F_r);
            #endif
            
            /* Total acceleration. There is no mass_water term because it is 
             * taken to be exactly 1 by our definition of units.*/
            double A = (F_c + F_d + F_r);

            /* Force by j on water i. */
            dyn_vars->watacc[3*i]   += A * Rij[0];
            dyn_vars->watacc[3*i+1] += A * Rij[1];
            dyn_vars->watacc[3*i+2] += A * Rij[2];

            /* Reaction force by i on water j. */
            dyn_vars->watacc[3*j]   -= A * Rij[0];
            dyn_vars->watacc[3*j+1] -= A * Rij[1];
            dyn_vars->watacc[3*j+2] -= A * Rij[2];
        }
    }

    /* -------------------------------------------------------------------- */
    /* 2. Particle-Particle interactions. */
    /* -------------------------------------------------------------------- */

    i = 0; /* "host" object counter */
    double e_cc = strength * in.E_CC;
    double s_cc = strength * in.SIGMA_CC;

    /* Traverse through the array until the "host" object has reached the 
     * last object.
     */
    for (int neigh_counter = 0; i < in.N_PARTICLES; neigh_counter++) {

        int j = neigh_list_p[neigh_counter];

        /* Check that we don't exceed the array. */
        assert(neigh_counter < in.N_PARTICLES * (in.N_PARTICLES + 1) / 2);
        /* We should never reach the -255 part because we should break out using 
         * the i < in.N_PARTICLES condition by then. 
         */
        assert (j != -255);
        /* Each entry is either a boundary marker (-1) or is an index of a 
         * particle, and so lies between 0 and N_PARTICLES-1.
         */
        assert (j == -1 || (j >= 0 && j < in.N_PARTICLES)); 

        #if DEBUG
            printf("counter = %d, i = %d, j = %d\n", neigh_counter, i, j);
        #endif

        /* A neighbour of -1 means that we have reached the end of this host. */
        if (j == -1) {
            i++;
            continue;
        }

        /* Rij stores the relative displacement vector. Points from j to i. */
        double Rij[3];
        get_rel_vector(in, dyn_vars->t, dyn_vars->partpos, dyn_vars->partpos, i, j, Rij);

        /* Find the squared distance between particles. */
        double dist_sq = Rij[0] * Rij[0] + Rij[1] * Rij[1] + Rij[2] * Rij[2];

        #if VERBOSE
            printf("Particles (%d, %d), vector from j to i: (%f,%f,%f), dist_sq: %f\n", 
                i, j, Rij[0], Rij[1], Rij[2], dist_sq);
        #endif
        
        /* Cutoff distance is R_CC. */
        if (dist_sq < R_CC_sq) {

            /* Calculate the relative distance and the factor that models
             * how force decays (linearly) with distance. 
             */
            double dist = sqrt(dist_sq);
            double dist_weight = (1 - dist/in.R_CC);

            assert(dist_weight > 0);
            assert(dist_weight <= 1);

            /* Normalise the separation vector. */
            Rij[0] /= dist;
            Rij[1] /= dist;
            Rij[2] /= dist;

            /* Find the relative velocity, Vi - Vj. It is the velocity of i
             * relative to / from the frame of j. 
             */
            double Vij[3];
            Vij[0] = dyn_vars->partvel[3*i]   - dyn_vars->partvel[3*j];
            Vij[1] = dyn_vars->partvel[3*i+1] - dyn_vars->partvel[3*j+1];
            Vij[2] = dyn_vars->partvel[3*i+2] - dyn_vars->partvel[3*j+2];

            double V_dot_R =    Vij[0] * Rij[0] + 
                                Vij[1] * Rij[1] + 
                                Vij[2] * Rij[2];

            #if VERBOSE
                printf("Vi_x: %f, Vj_x: %f\n", dyn_vars->partvel[3*i], dyn_vars->partvel[3*j]);
                printf("Vij_x %f\n", Vij[0]);
                printf("Mod Vij: %f\n", Vij[0] * Vij[0] + Vij[1] * Vij[1] + Vij[2] * Vij[2]);
            #endif  

            /* Conservative Force. */
            double F_c =  -12 * e_cc * pow(s_cc, 6) / pow(dist, 8) * (1 - pow(s_cc / dist, 6));
            /* Dissipative/Damping Force. */
            double F_d = - in.DAMP_CONST * dist_weight * dist_weight * V_dot_R;
            /* Random/Brownian Force.
             * Variance of a uniform distribution is 1/12 * interval^2.
             * Therefore, for mean = 0 and variance = 1, assume we have -x/2
             * to x/2, then 1/12 * x^2 = 1 ==> x = sqrt(12) ~ 3.464. Note 
             * that we first scale from 0 ~ RAND_MAX ==> 0 ~ 1 ==> 
             * -1/2 ~ 1/2 ==> -sqrt(12)/2 ~ sqrt(12)/2.
             */
            double F_r = sigma * dist_weight * 3.464 * ((rand()/(double) RAND_MAX) - 0.5) / sqrt_dt;

            #if VERBOSE
                printf("Particles (%d, %d), F_c: %f, F_d: %f, F_r: %f\n", i, j, F_c, F_d, F_r);
            #endif
            
            /* Total acceleration = Total F / m. */
            double A = (F_c + F_d + F_r) / in.M_PARTICLE;

            /* Force by j on particle i. */
            dyn_vars->partacc[3*i]   += A * Rij[0];
            dyn_vars->partacc[3*i+1] += A * Rij[1];
            dyn_vars->partacc[3*i+2] += A * Rij[2];

            /* Reaction force by i on particle j. */
            dyn_vars->partacc[3*j]   -= A * Rij[0];
            dyn_vars->partacc[3*j+1] -= A * Rij[1];
            dyn_vars->partacc[3*j+2] -= A * Rij[2];
        }
    }

    /* -------------------------------------------------------------------- */
    /* 3. Particle-Water interactions. */ 
    /* -------------------------------------------------------------------- */

    i = 0; /* "host" object counter. Host is particle in this case. */
    double e_sc = strength * in.E_SC;
    double s_sc = strength * in.SIGMA_SC;

    /* Traverse through the array until the "host" object has reached the 
     * last object.
     */
    for (int neigh_counter = 0; i < in.N_PARTICLES; neigh_counter++) {

        int j = neigh_list_pw[neigh_counter];

        /* Check that we don't exceed the array. */
        assert(neigh_counter < in.N_PARTICLES * (in.N_WATER + 1));
        /* We should never reach the -255 part because we should break out using 
         * the i < in.N_PARTICLES condition by then. 
         */
        assert (j != -255);
        /* Each entry is either a boundary marker (-1) or is an index of a 
         * water, and so lies between 0 and N_WATER-1.
         */
        assert (j == -1 || (j >= 0 && j < in.N_WATER)); 

        #if DEBUG
            printf("counter = %d, i = %d, j = %d\n", neigh_counter, i, j);
        #endif

        /* A neighbour of -1 means that we have reached the end of this host. */
        if (j == -1) {
            i++;
            continue;
        }

        /* Rij stores the relative displacement vector. Points from j to i. */
        double Rij[3];
        get_rel_vector(in, dyn_vars->t, dyn_vars->partpos, dyn_vars->watpos, i, j, Rij);

        /* Find the squared distance between particles. */
        double dist_sq = Rij[0] * Rij[0] + Rij[1] * Rij[1] + Rij[2] * Rij[2];

        #if VERBOSE
            printf("Particles (%d, %d), vector from j to i: (%f,%f,%f), dist_sq: %f\n", 
                i, j, Rij[0], Rij[1], Rij[2], dist_sq);
        #endif
        
        /* Cutoff distance is R_SC. */
        if (dist_sq < R_SC_sq) {

            /* Calculate the relative distance and the factor that models
             * how force decays (linearly) with distance. 
             */
            double dist = sqrt(dist_sq);
            double dist_weight = (1 - dist/in.R_SC);

            assert(dist_weight > 0);
            assert(dist_weight <= 1);

            /* Normalise the separation vector. */
            Rij[0] /= dist;
            Rij[1] /= dist;
            Rij[2] /= dist;

            /* Find the relative velocity, Vi - Vj. It is the velocity of i
             * relative to / from the frame of j. i is a particle, j  is water.
             */
            double Vij[3];
            Vij[0] = dyn_vars->partvel[3*i]   - dyn_vars->watvel[3*j];
            Vij[1] = dyn_vars->partvel[3*i+1] - dyn_vars->watvel[3*j+1];
            Vij[2] = dyn_vars->partvel[3*i+2] - dyn_vars->watvel[3*j+2];

            double V_dot_R =    Vij[0] * Rij[0] + 
                                Vij[1] * Rij[1] + 
                                Vij[2] * Rij[2];

            #if VERBOSE
                printf("Vi_x: %f, Vj_x: %f\n", dyn_vars->partvel[3*i], dyn_vars->partvel[3*j]);
                printf("Vij_x %f\n", Vij[0]);
                printf("Mod Vij: %f\n", Vij[0] * Vij[0] + Vij[1] * Vij[1] + Vij[2] * Vij[2]);
            #endif  

            /* Conservative Force. */
            double F_c =  -12 * e_sc * pow(s_sc, 6) / pow(dist, 8) * (1 - pow(s_sc / dist, 6));
            /* Dissipative/Damping Force. */
            double F_d = - in.DAMP_CONST * dist_weight * dist_weight * V_dot_R;
            /* Random/Brownian Force.
             * Variance of a uniform distribution is 1/12 * interval^2.
             * Therefore, for mean = 0 and variance = 1, assume we have -x/2
             * to x/2, then 1/12 * x^2 = 1 ==> x = sqrt(12) ~ 3.464. Note 
             * that we first scale from 0 ~ RAND_MAX ==> 0 ~ 1 ==> 
             * -1/2 ~ 1/2 ==> -sqrt(12)/2 ~ sqrt(12)/2.
             */
            double F_r = sigma * dist_weight * 3.464 * ((rand()/(double) RAND_MAX) - 0.5) / sqrt_dt;

            #if VERBOSE
                printf("Particles (%d, %d), F_c: %f, F_d: %f, F_r: %f\n", i, j, F_c, F_d, F_r);
            #endif
            
            /* Total acceleration = Total F / m. */
            double A_w = (F_c + F_d + F_r);
            double A_p = A_w / in.M_PARTICLE;

            /* Force by water j on particle i. */
            dyn_vars->partacc[3*i]   += A_p * Rij[0];
            dyn_vars->partacc[3*i+1] += A_p * Rij[1];
            dyn_vars->partacc[3*i+2] += A_p * Rij[2];

            /* Reaction force by particle i on water j. */
            dyn_vars->watacc[3*j]   -= A_w * Rij[0];
            dyn_vars->watacc[3*j+1] -= A_w * Rij[1];
            dyn_vars->watacc[3*j+2] -= A_w * Rij[2];
        }
    }

}

/* Gives a vector pointing from j to i, which is computed by Ri - Rj. 
 * Enforces periodic BC by ensuring that all separations are at most L/2. 
 */
void get_rel_vector(Inputs in, int t, double *pos_list_1, double *pos_list_2, int i, int j, double *Rij) {

    double x_shift = 0;

    /* Takes Ri - Rj for the y-axis. */
    Rij[1] = pos_list_1[3*i+1] - pos_list_2[3*j+1];

    /* Enforces periodic BC's nearest image for y-direction first. 
     * If we're folding in the y-direction, imagine that i is near the 
     * top and j is near the bottom, so Rij (pointing from j to i) is
     * positive and larger than size / 2. Hence, we fold this by 
     * bringing particle i to the box below the main box, hence the 
     * vector is reduced by size. Moreover, since the bottom box is
     * moving to the left, particle i (and therefore also the vector)
     * is shifted left (negative) by the distance travelled.
     */
    if (Rij[1] > in.BOX_SIZE/2) {
        Rij[1] -= in.BOX_SIZE;
        x_shift = - t * in.V_SHEAR * in.TIME_STEP - 
                  floor(t * in.V_SHEAR * in.TIME_STEP / in.BOX_SIZE) * in.BOX_SIZE;
    }
    else if (Rij[1] < -in.BOX_SIZE/2) {
        Rij[1] += in.BOX_SIZE;
        x_shift = t * in.V_SHEAR * in.TIME_STEP - 
                  floor(t * in.V_SHEAR * in.TIME_STEP / in.BOX_SIZE) * in.BOX_SIZE;
    }

    /* Repeat the finding of nearest neighbours for x-axis. */
    Rij[0] = pos_list_1[3*i] - pos_list_2[3*j];
    /* Put in the effect of y boundary traversal. */
    Rij[0] += x_shift;

    /* Need this while loop because the x_shift might cause Rij[0] (i.e. the x
     * component) to be more than 1 box size off from the principle box, so we
     * might need to do the correction > 1 time.
     */
    while (Rij[0] > in.BOX_SIZE/2 || Rij[0] < -in.BOX_SIZE/2) {
        if (Rij[0] > in.BOX_SIZE/2) {
        Rij[0] -= in.BOX_SIZE;
        }
        else if (Rij[0] < -in.BOX_SIZE/2) {
            Rij[0] += in.BOX_SIZE;
        }
    }
    
    /* Repeat the finding of nearest neighbours for z-axis. */
    Rij[2] = pos_list_1[3*i+2] - pos_list_2[3*j+2];

    if (Rij[2] > in.BOX_SIZE/2) {
        Rij[2] -= in.BOX_SIZE;
    }
    else if (Rij[2] < -in.BOX_SIZE/2) {
        Rij[2] += in.BOX_SIZE;
    }

    #if DEBUG
        printf("t:%d, i:%d, j:%d, %f, %f, %f\n", t, i, j, Rij[0], Rij[1], Rij[2]);
    #endif
}

void update_neigh_list(Inputs in, int t, double *pos_list, double cutoff_dist, 
                            int *neigh_list, int num_obj) {

    /* Points to the next empty slot in the neighbour list. */
    int list_pointer = 0;
    int neigh_arr_length = num_obj * (num_obj + 1) / 2;
    double neigh_dist_sq = (cutoff_dist + SKIN) * (cutoff_dist + SKIN);

    /* Set all values in the neigh_list to be a conspicuous vale of -255. */
    for (int i = 0; i < neigh_arr_length; i++) {
        neigh_list[i] = -255;
    }

    for (int i = 0; i < num_obj; i++) {
        for (int j = i+1; j < num_obj; j++) {

            /* Check that we didn't exceed the array. */
            assert(list_pointer < neigh_arr_length);

            /* Rij stores the relative displacement vector. Points from j to i. */
            double Rij[3];
            get_rel_vector(in, t, pos_list, pos_list, i, j, Rij);

            /* Find the squared distance between particles. */
            double dist_sq = Rij[0] * Rij[0] + Rij[1] * Rij[1] + Rij[2] * Rij[2];

            /* Because each dimension is at most L/2, so total squared norm is 
             * at most 3 times of that squared.
             */
            assert(dist_sq < 3 * (in.BOX_SIZE/2) * (in.BOX_SIZE/2));

            /* If within the interaction zone + buffer, we save that neighbour. */
            if (dist_sq < neigh_dist_sq) {
                neigh_list[list_pointer] = j;
                list_pointer++;
            }
        }

        /* Put a -1 after each "host" object's neighbours. */
        neigh_list[list_pointer] = -1;
        list_pointer++;
    }
}

void update_neigh_list_pw(Inputs in, int t, 
                        double *pos_list_p, double *pos_list_w, 
                        double cutoff_dist, int *neigh_list_pw, 
                        int num_p, int num_w) {

    /* Points to the next empty slot in the neighbour list. */
    int list_pointer = 0;
    int neigh_arr_length = num_p * (num_w + 1);
    double neigh_dist_sq = (cutoff_dist + SKIN) * (cutoff_dist + SKIN);

    /* Set all values in the neigh_list_pw to be a conspicuous vale of -255. */
    for (int i = 0; i < neigh_arr_length; i++) {
        neigh_list_pw[i] = -255;
    }

    /* Loop over particles */
    for (int i = 0; i < num_p; i++) {
        /* Loop over water */
        for (int j = 0; j < num_w; j++) {

            /* Check that we didn't exceed the array. */
            assert(list_pointer < neigh_arr_length);

            /* Rij stores the relative displacement vector. Points from j to i. */
            double Rij[3];
            get_rel_vector(in, t, pos_list_p, pos_list_w, i, j, Rij);

            /* Find the squared distance between particles. */
            double dist_sq = Rij[0] * Rij[0] + Rij[1] * Rij[1] + Rij[2] * Rij[2];

            /* Because each dimension is at most L/2, so total squared norm is 
             * at most 3 times of that squared.
             */
            assert(dist_sq < 3 * (in.BOX_SIZE/2) * (in.BOX_SIZE/2));

            /* If within the interaction zone + buffer, we save that neighbour. */
            if (dist_sq < neigh_dist_sq) {
                neigh_list_pw[list_pointer] = j;
                list_pointer++;
            }
        }

        /* Put a -1 after each "host" object's neighbours. */
        neigh_list_pw[list_pointer] = -1;
        list_pointer++;
    }
}

int check_update_req(double *disp_list_w, double *disp_list_p, int n_w, int n_p) {
    double largest     = 0;
    double sec_largest = 0;

    /* Go through both lists, find the largest and second largest components of 
     * displacements.
     */
    for (int i = 0; i < 3*n_w; i++) {
        if (disp_list_w[i] > largest) {
            sec_largest = largest;
            largest = disp_list_w[i];
        }
    }

    for (int i = 0; i < 3*n_p; i++) {
        if (disp_list_p[i] > largest) {
            sec_largest = largest;
            largest = disp_list_p[i];
        }
    } 

    #if DEBUG
        printf("Largest disp = %f, Second Largest = %f\n", largest, sec_largest);
        if ((largest + sec_largest) > SKIN) {
            printf("EXCEED!!!\n\n");
        }
    #endif

    /* Need to update the list if the 2 largest displacements are greater than 
     * SKIN, the buffer distance in our neighbour list construction.
     */
    if ((largest + sec_largest) > SKIN) {
        return 1;
    }
    else {
        return 0;
    }
}

/* Calculates temperature by using E = 3/2 NkT = sum(1/2 * m * v^2). */
double calc_temp(Dyn_Vars *dyn_vars, Inputs in) {

    double KE = 0;

    /* Sum the KE from all particles. No mass_water term because it is taken to
     * be exactly 1 by our definition of units.
     */
    for (int i = 0; i < in.N_WATER; i++) {

        int x = 3 * i;
        int y = x + 1;
        int z = y + 1;

        /* Subtract away the x-velocity from the shear flow. */
        double vx_corr = dyn_vars->watvel[x] - in.V_SHEAR * (dyn_vars->watpos[y] / in.BOX_SIZE - 0.5);

        KE += 0.5 * 
                (vx_corr * vx_corr + 
                 dyn_vars->watvel[y] * dyn_vars->watvel[y] +
                 dyn_vars->watvel[z] * dyn_vars->watvel[z]);
    }

    double temp = KE / (3.0/2.0 * in.N_WATER); /* E ~ 3/2 NkT */

    return temp;
}