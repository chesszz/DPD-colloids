#include "utils.h"
#include "sim.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define SKIN 0.8

void refold_positions(Dyn_Vars *dyn_vars, Inputs in);
void calculate_acc(Dyn_Vars *dyn_vars, Inputs in, int *neigh_list);
void get_rel_vector(Dyn_Vars *dyn_vars, Inputs in, int i, int j, double *Rij);
void update_neigh_list(Dyn_Vars *dyn_vars, Inputs in, int *neigh_list, double *disp_list, int num_obj);
int check_update_req(double *disp_list, int num_obj);
double calc_temp(Dyn_Vars *dyn_vars, Inputs in);

/* Evolves the system from the given initial conditions and input parameters.
 * Will print the variables' time evolution into a file if PRINT_STEPS is 
 * defined to be 1 in utils.h.
 */ 
void evolve_system(Dyn_Vars *dyn_vars, Inputs in) {

    #if DEBUG
        int list_recomp_count = 0;
    #endif

    /* Used to accumulate the displacements each time step. */
    double disp;
    double dt = in.TIME_STEP;
    int checkpoint = in.N_STEPS / 10;
    FILE *dyn_out = fopen("dyn.out", "w");
    FILE *temp_out = fopen("temp.out", "w");

    int update_req_w = 1;
    int update_req_p = 1;

    /* List of displacements for each water since the last update. */
    double *disp_list_w = (double *) calloc(3 * in.N_WATER, sizeof(double));
    /* List of displacements for each particle since the last update. */
    double *disp_list_p = (double *) calloc(3 * in.N_PARTICLES, sizeof(double));

    /* List of neighbours for each water. Boundaries are marked by -1. Note
     * that because we only consider pairs of waters i, j where i < j, thus
     * we have at most N(N-1)/2 pairs. However, we also need ~N dividers, hence
     * we have a total of N(N+1)/2 entries.*/
    int *neigh_list_w = (int *) calloc(in.N_WATER * (in.N_WATER + 1) / 2, sizeof(int));
    /* List of neighbours for each particle. */
    int *neigh_list_p = (int *) calloc(in.N_PARTICLES * (in.N_PARTICLES + 1) / 2, sizeof(int));

    if (disp_list_w == NULL || disp_list_p == NULL || neigh_list_w == NULL || neigh_list_p == NULL) {
        fprintf(stderr, "Error in allocating memory.\n");
        exit(1);
    }

    #if PRINT_STEPS
        fprintf(dyn_out, "Time Step, Particle Number, Particle Type, Position, Velocity, Acceleration\n");
    #endif

    /* MAIN TIME LOOP */
    for (int t = 0; t < in.N_STEPS; t++) {

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
            disp_list_w[i]      += disp;
        }
        for (int i = 0; i < 3*in.N_PARTICLES; i++) {
            disp = dt * dyn_vars->partvel[i];
            /* Add displacement to both the position and displacement list. */
            dyn_vars->partpos[i] += disp;
            disp_list_p[i]       += disp;
        }

        /* Wrap the objects back to the box if they escaped. Handles both water
         * and particles. 
         */
        refold_positions(dyn_vars, in);

        /* Update neighbour list if required. */
        if (update_req_w) {
            update_neigh_list(dyn_vars, in, neigh_list_w, disp_list_w, in.N_WATER);
            update_req_w = 0;

            #if DEBUG
                printf("List recomputed: %d\n", list_recomp_count);
                list_recomp_count++;
            #endif          
        }
        if (update_req_p) {
            update_neigh_list(dyn_vars, in, neigh_list_p, disp_list_p, in.N_PARTICLES);
            update_req_p = 0;
        }

        /* Velocity Verlet 3: F(t+dt) */
        calculate_acc(dyn_vars, in, neigh_list_w);
        //calculate_acc(dyn_vars, in, neigh_list_p);

        /* Velocity Verlet 4: V(t+dt) */
        for (int i = 0; i < 3*in.N_WATER; i++) {
            dyn_vars->watvel[i] += 0.5 * dt * dyn_vars->watacc[i];
        }
        for (int i = 0; i < 3*in.N_PARTICLES; i++) {
            dyn_vars->partvel[i] += 0.5 * dt * dyn_vars->partacc[i];
        }

        /* Check if the particles have moved too much since the last update
         * that requires a new neighbour list to be created.
         */
        update_req_w = check_update_req(disp_list_w, in.N_WATER);
        //update_req_p = check_update_req(disp_list_p, in.N_PARTICLES);

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

        #if PRINT_STEPS
            /* At each time step, print the positions, velocities, and 
             * accelerations for each particle.
             */
            for (int i = 0; i < in.N_WATER; i++) {
                fprintf(dyn_out, "%d, %d, %d, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", 
                    t, i, 0, /* 0 means water */
                    dyn_vars->watpos[3*i], dyn_vars->watpos[3*i+1], dyn_vars->watpos[3*i+2], 
                    dyn_vars->watvel[3*i], dyn_vars->watvel[3*i+1], dyn_vars->watvel[3*i+2],
                    dyn_vars->watacc[3*i], dyn_vars->watacc[3*i+1], dyn_vars->watacc[3*i+2]);
            }
            for (int i = 0; i < in.N_PARTICLES; i++) {
                fprintf(dyn_out, "%d, %d, %d, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", 
                    t, i, 1, /* 1 means particle */
                    dyn_vars->partpos[3*i], dyn_vars->partpos[3*i+1], dyn_vars->partpos[3*i+2], 
                    dyn_vars->partvel[3*i], dyn_vars->partvel[3*i+1], dyn_vars->partvel[3*i+2],
                    dyn_vars->partacc[3*i], dyn_vars->partacc[3*i+1], dyn_vars->partacc[3*i+2]);
            }

            /* At each time step, print the temperature. */
            fprintf(temp_out, "%d, %f\n", t, calc_temp(dyn_vars, in));
        #endif

        /* Print a message every 10% complete. */
        if (t % checkpoint == 0) {
            fprintf(stderr, "Time step %d / %d complete.\n", t, in.N_STEPS);
        }  
    }

    fclose(dyn_out);
    fclose(temp_out);

    free(disp_list_w);
    free(disp_list_p);
    free(neigh_list_w);
    free(neigh_list_p);
}

/* Enforce periodic BCs. TODO: Lees-Edwards. */
void refold_positions(Dyn_Vars *dyn_vars, Inputs in) {

    double size = in.BOX_SIZE;

    /* If water escapes from the left or right side of the box, wrap it back
     * to the inside. Each dimension is handled independently.
     */
    for (int i = 0; i < 3*in.N_WATER; i++) {
        if (dyn_vars->watpos[i] > size) {
            dyn_vars->watpos[i] -= size;
        }
        else if (dyn_vars->watpos[i] < 0) {
            dyn_vars->watpos[i] += size;
        }
    }

    for (int i = 0; i < 3*in.N_PARTICLES; i++) {
        if (dyn_vars->partpos[i] > size) {
            dyn_vars->partpos[i] -= size;
        }
        else if (dyn_vars->partpos[i] < 0) {
            dyn_vars->partpos[i] += size;
        }
    }
}

/* Calculate the forces acting on each particle. */
void calculate_acc(Dyn_Vars *dyn_vars, Inputs in, int *neigh_list) {

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

    /* Water-Water interactions. */

    int i = 0; /* "host" object counter */

    /* Traverse through the array until the "host" object has reached the 
     * last object.
     */
    for (int neigh_counter = 0; i < in.N_WATER; neigh_counter++) {

        int j = neigh_list[neigh_counter];

        /* Check that we don't exceed the array. */
        assert(neigh_counter < in.N_WATER * (in.N_WATER + 1) / 2);
        /* We should never reach the -255 part because we should break out using 
         * the i < in.WATER condition by then. 
         */
        assert (j != -255);
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
        get_rel_vector(dyn_vars, in, i, j, Rij);

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

    /* TODO: Water-Particle interactions and Particle-Particle interactions. */
}

/* Gives a vector pointing from j to i, which is computed by Ri - Rj. 
 * Enforces periodic BC by ensuring that all separations are at most L/2. 
 */
void get_rel_vector(Dyn_Vars *dyn_vars, Inputs in, int i, int j, double *Rij) {

    /* For all 3 dimensions. */
    for (int k = 0; k < 3; k++) {

        /* Takes Ri - Rj. */
        Rij[k] = dyn_vars->watpos[3*i+k] - dyn_vars->watpos[3*j+k];

        /* Enforces periodic BC's nearest image. */
        if (Rij[k] > in.BOX_SIZE/2) {
            Rij[k] -= in.BOX_SIZE;
        }
        else if (Rij[k] < -in.BOX_SIZE/2) {
            Rij[k] += in.BOX_SIZE;
        }
    }
}

void update_neigh_list(Dyn_Vars *dyn_vars, Inputs in, int *neigh_list, double *disp_list, int num_obj) {

    /* Points to the next empty slot in the neighbour list. */
    int list_pointer = 0;
    int neigh_arr_length = num_obj * (num_obj + 1) / 2;
    double neigh_dist_sq = (1 + SKIN) * (1 + SKIN);

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
            get_rel_vector(dyn_vars, in, i, j, Rij);

            /* Find the squared distance between particles. */
            double dist_sq = Rij[0] * Rij[0] + Rij[1] * Rij[1] + Rij[2] * Rij[2];

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

    /* Set the displacement list to be 0. */
    for (int i = 0; i < 3*num_obj; i++) {
        disp_list[i] = 0;
    }
}

int check_update_req(double *disp_list, int num_obj) {
    double largest     = 0;
    double sec_largest = 0;

    /* Go through the list, find the largest and second largest components of 
     * displacements.
     */
    for (int i = 0; i < 3*num_obj; i++) {
        if (disp_list[i] > largest) {
            sec_largest = largest;
            largest = disp_list[i];
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
        KE += 0.5 * 
                (dyn_vars->watvel[3*i+0] * dyn_vars->watvel[3*i+0] + 
                 dyn_vars->watvel[3*i+1] * dyn_vars->watvel[3*i+1] +
                 dyn_vars->watvel[3*i+2] * dyn_vars->watvel[3*i+2]);
    }

    double temp = KE / (3.0/2.0 * in.N_WATER); /* E ~ 3/2 NkT */

    return temp;
}