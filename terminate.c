#include <stdlib.h>
#include <stdio.h>
#include "sim.h"

/* Prints some stats on the screen. */
void output_state(Dyn_Vars *dyn_vars, Inputs in) {
    fprintf(stderr, "Number of water particles: %d\n", in.N_WATER);
    fprintf(stderr, "Number of hard particles: %d, Mass: %f\n", in.N_PARTICLES, in.M_PARTICLE);
    fprintf(stderr, "Box Size: %f, Shear Velocity: %f\n", in.BOX_SIZE, in.V_SHEAR);
    fprintf(stderr, "Damping Constant: %f, Spring Constant: %f\n", in.DAMP_CONST, in.SPRING_CONST);
    fprintf(stderr, "Simulated %d time steps of length %f\n", in.N_STEPS, in.TIME_STEP);
}

/* Frees all the dynamical variables lists and the struct itself. */
void terminate(Dyn_Vars *dyn_vars) {

    free(dyn_vars->watpos);
    free(dyn_vars->watvel); 
    free(dyn_vars->watacc);
    free(dyn_vars->partpos);
    free(dyn_vars->partvel);
    free(dyn_vars->partacc);
    free(dyn_vars);

}