#ifndef SIM_H
#define SIM_H

#define DEBUG 0
#define VERBOSE 0

/* SHOULD BE 0 */
#define TWO_D 0

/* SHOULD BE 0 */
#define PRINT_WATER 0
#define PRINT_PARTICLES 0

/* SHOULD BE 1 */
#define OFF_AXIS_FORCE_ON 1

/* SHOULD BE 1. Can be >1 but not 0. */
#define NUM_PROG_STEPS 1

/* Used to contain all the input parameters to be passed around. */
typedef struct {
    int    N_WATER;
    int    N_PARTICLES;

    int    N_STEPS;
    double TIME_STEP;
    double BOX_SIZE;
    double SHEAR_RATE;
    
    double DAMP_CONST;
    double SPRING_CONST;

    double R_PARTICLE_AVG;
    double R_PARTICLE_STDEV;
    double E_SC;
    double E_CC;

    double *PART_RADII;
    double *PART_MASS;
} Inputs;

/* Contains all the dynamical variable arrays + time counter. */
typedef struct {
    int t;
    double *watpos;
    double *watvel;
    double *watacc;
    double *watomg;
    double *watalp;

    double *partpos;
    double *partvel;
    double *partacc;
    double *partomg;
    double *partalp;
} Dyn_Vars;

/* Evolves the system from the given initial conditions and input parameters.
 * Will print the variables' time evolution into a file if PRINT_STEPS is 
 * defined to be 1 in utils.h.
 */ 
void evolve_system(Dyn_Vars *dyn_vars, Inputs in);
void calc_viscosity(Inputs in);
void calc_viscosity_book(Inputs in);

#endif /* SIM_H */
