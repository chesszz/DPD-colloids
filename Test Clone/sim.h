#ifndef SIM_H
#define SIM_H

#define DEBUG 0
#define VERBOSE 0

/* SHOULD BE 0 */
#define TWO_D 1

/* SHOULD BE 0 */
#define PRINT_WATER 1
#define PRINT_PARTICLES 1

/* SHOULD BE 1 */
#define OFF_AXIS_FORCE_ON 0

/* Used to contain all the input parameters to be passed around. */
typedef struct {
    int N_WATER;
    int N_PARTICLES;

    int N_STEPS;
    double TIME_STEP;
    double BOX_SIZE;
    double SHEAR_RATE;
    
    double DAMP_CONST;
    double SPRING_CONST;

    double M_PARTICLE;
    double R_SC;
    double R_CC;
    double E_SC;
    double SIGMA_SC;
    double E_CC;
    double SIGMA_CC;
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
void calc_viscosity(Inputs in, int type);

#endif /* SIM_H */
