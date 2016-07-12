#ifndef UTILS_H
#define UTILS_H

#define DEBUG 0
#define VERBOSE 0

#define TWO_D 1

#define PRINT_WATER 1
#define PRINT_PARTICLES 1

/* Used to contain all the input parameters to be passed around. */
typedef struct {
    int N_WATER;
    int N_PARTICLES;

    int N_STEPS;
    double TIME_STEP;
    double BOX_SIZE;
    double V_SHEAR;
    
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
    double *wattoq;

    double *partpos;
    double *partvel;
    double *partacc;
    double *partomg;
    double *parttoq;
} Dyn_Vars;

#endif /* UTILS_H */