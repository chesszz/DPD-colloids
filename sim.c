#include "utils.h"
#include "sim.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

void refold_positions(Dyn_Vars *dyn_vars, Inputs in);
void calculate_acc(Dyn_Vars *dyn_vars, Inputs in);
void get_rel_vector(Dyn_Vars *dyn_vars, Inputs in, int i, int j, double *Rij);
double calc_temp(Dyn_Vars *dyn_vars, Inputs in);

/* Evolves the system from the given initial conditions and input parameters.
 * Will print the variables' time evolution into a file if PRINT_STEPS is 
 * defined to be 1 in utils.h.
 */ 
void evolve_system(Dyn_Vars *dyn_vars, Inputs in) {

    double dt = in.TIME_STEP;
    int checkpoint = in.N_STEPS / 10;
    FILE *dyn_out = fopen("dyn.out", "w");
    FILE *temp_out = fopen("temp.out", "w");

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
            dyn_vars->watpos[i] += dt * dyn_vars->watvel[i];
        }
        for (int i = 0; i < 3*in.N_PARTICLES; i++) {
            dyn_vars->partpos[i] += dt * dyn_vars->partvel[i];
        }

        /* Wrap the particles back to the box if they escaped. */
        refold_positions(dyn_vars, in);

        /* Velocity Verlet 3: F(t+dt) */
        calculate_acc(dyn_vars, in);

        /* Velocity Verlet 4: V(t+dt) */
        for (int i = 0; i < 3*in.N_WATER; i++) {
            dyn_vars->watvel[i] += 0.5 * dt * dyn_vars->watacc[i];
        }
        for (int i = 0; i < 3*in.N_PARTICLES; i++) {
            dyn_vars->partvel[i] += 0.5 * dt * dyn_vars->partacc[i];
        }

        /* Print out the positions and velocities of everything for each t. */
        #if DEBUG
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
}

/* Enforce periodic BCs. TODO: Lees-Edwards. */
void refold_positions(Dyn_Vars *dyn_vars, Inputs in) {

    double size = in.BOX_SIZE;

    /* If water escapes from the left or right side of the box, wrap it back
     * to the inside. Each dimension is handled independetly.
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
void calculate_acc(Dyn_Vars *dyn_vars, Inputs in) {

    /* Zero out all accelerations. */
    for (int i = 0; i < 3*in.N_WATER; i++) {
        dyn_vars->watacc[i] = 0;
    }
    for (int i = 0; i < 3*in.N_PARTICLES; i++) {
        dyn_vars->partacc[i] = 0;
    }

    /* TODO: Change to something O(N), eg. fincham ralston loop */

    /* Calculate sqrt of dt for use in F_r. Calculate it just once here. */
    double sqrt_dt= sqrt(in.TIME_STEP);
    /* Also calculate sigma, the Brownian noise strength. */
    double sigma =  sqrt(2 * in.DAMP_CONST * in.TEMP);

    /* Water-Water interactions. */
    for (int i = 0; i < in.N_WATER; i++) {
        for (int j = i+1; j < in.N_WATER; j++) {

            /* Rij stores the relative displacement vector. Points from j to i. */
            double Rij[3];
            get_rel_vector(dyn_vars, in, i, j, Rij);

            /* Find the squared distance between particles. */
            double dist_sq = Rij[0] * Rij[0] + Rij[1] * Rij[1] + Rij[2] * Rij[2];

            #if DEBUG
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

                #if DEBUG
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

                #if DEBUG
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