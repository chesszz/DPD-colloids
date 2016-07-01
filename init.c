#include "utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define TWO_D 0

void fill_gaussian_two_tuple(double *tuple, double stdev);

/* Reads inputs from stdin that correspond to the vaious parameters. Returns
 * an Inputs struct that stores those parameters.
 */
Inputs get_inputs(void) {

    Inputs in;

    scanf("%d %d %d %lf %lf %lf %lf %lf %lf", 
        &in.N_WATER, &in.N_PARTICLES,
        &in.N_STEPS, &in.TIME_STEP, &in.BOX_SIZE, &in.V_SHEAR,
        &in.DAMP_CONST, &in.SPRING_CONST, &in.M_PARTICLE);

    return in;
}

/* Allocates space for the dynamical variables - position, velocity, 
 * acceleration. Fills the positions randomly within the box, and then fills
 * the velocities with a Maxwell-Boltzmann distribution. Returns a Dyn_Vars
 * pointer that stores these arrays.
 */
Dyn_Vars *initialise(Inputs in) {

    /* Space for the struct that holds all the dynamical variables. */
     Dyn_Vars *dyn_vars = (Dyn_Vars *) malloc(sizeof(Dyn_Vars));

    if (dyn_vars == NULL) {
        fprintf(stderr, "Error in memory allocation!\n");
        exit(1);
    }

    dyn_vars->t = 0;

    /* Allocates arrays of length 3 * N_WATER and 3 * N_PARTICLES. */
    dyn_vars->watpos  = (double *) calloc(3 * in.N_WATER,     sizeof(double));
    dyn_vars->watvel  = (double *) calloc(3 * in.N_WATER,     sizeof(double));
    dyn_vars->watacc  = (double *) calloc(3 * in.N_WATER,     sizeof(double));
    dyn_vars->partpos = (double *) calloc(3 * in.N_PARTICLES, sizeof(double));
    dyn_vars->partvel = (double *) calloc(3 * in.N_PARTICLES, sizeof(double));
    dyn_vars->partacc = (double *) calloc(3 * in.N_PARTICLES, sizeof(double));

    if (dyn_vars->watpos == NULL || dyn_vars->watvel == NULL || 
        dyn_vars->watacc == NULL || dyn_vars->partpos == NULL ||
        dyn_vars->partvel == NULL || dyn_vars->partacc == NULL) {
        fprintf(stderr, "Error in memory allocation!\n");
        exit(1);
    }

    // dyn_vars->watpos[0] = 1; // 1
    // dyn_vars->watpos[1] = 0.2;
    // dyn_vars->watpos[2] = 0;
    // dyn_vars->watpos[3] = 0; // 2
    // dyn_vars->watpos[4] = 1;
    // dyn_vars->watpos[5] = 0;
    // dyn_vars->watpos[6] = 4; // 3
    // dyn_vars->watpos[7] = 4.8;
    // dyn_vars->watpos[8] = 0;

    for (int i = 0; i < in.N_WATER; i++) {

        int x = 3 * i;
        int y = x + 1;
        int z = y + 1;

        /* Random positions from 0 to BOX_SIZE. */
        dyn_vars->watpos[x]   = in.BOX_SIZE * (double) rand() / RAND_MAX;
        dyn_vars->watpos[y] = in.BOX_SIZE * (double) rand() / RAND_MAX;
        dyn_vars->watpos[z] = in.BOX_SIZE * (double) rand() / RAND_MAX;

        #if TWO_D
            dyn_vars->watpos[z] = 0; /* TODO: Uncomment for 2D. */
        #endif

        /* http://scicomp.stackexchange.com/questions/19969/how-do-i-generate-maxwell-boltzmann-variates-using-a-uniform-distribution-random */
        /* Maxwell-Boltzmann distribution has all components of velocity
         * randomly drawn from a Gaussian. 
         */
        double gaussian_tuple[2];

        /* Variance is 1 since we take temp and mass to be 1 by definition of 
         * our units.
         */
        fill_gaussian_two_tuple(gaussian_tuple, 1);
        dyn_vars->watvel[x] = gaussian_tuple[0];
        dyn_vars->watvel[y] = gaussian_tuple[1];

        /* Add a x-velocity gradient in the y-direction from the shear. We have
         * a positive x-velocity at y = BOX_SIZE, and a negative x-velocity at
         * y = 0.
         */
        dyn_vars->watvel[x] += in.V_SHEAR * (dyn_vars->watpos[y] / in.BOX_SIZE - 0.5);

        fill_gaussian_two_tuple(gaussian_tuple, 1);
        dyn_vars->watvel[z] = gaussian_tuple[1];

        #if TWO_D
            dyn_vars->watvel[z] = 0; /* TODO: Uncomment for 2D. */
        #endif
    }
    
    return dyn_vars;
}

/* Take a tuple with at least 2 slots and fills them with Gaussian values that
 * have sigma = 1.
 */
void fill_gaussian_two_tuple(double *tuple, double stdev) {

    double x, y, s;

    do {
        x = (2.0 * rand() / RAND_MAX) - 1; 
        y = (2.0 * rand() / RAND_MAX) - 1; 
        s = x * x + y * y;
    } while (s >= 1);

    double z1 = x * sqrt(-2 * log(s) / s);
    double z2 = y * sqrt(-2 * log(s) / s);

    tuple[0] = z1*stdev;
    tuple[1] = z2*stdev;
}