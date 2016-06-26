#ifndef UTILS_H
#define UTILS_H

#define DEBUG 0
#define PRINT_STEPS 1

/* Used to contain all the input parameters to be passed around. */
typedef struct {
    int N_WATER;
    int N_PARTICLES;

    int N_STEPS;
    double TIME_STEP;
    double TEMP;
    double BOX_SIZE;
    double V_SHEAR;
    
    double DAMP_CONST;
    double SPRING_CONST;
    double M_PARTICLE;
} Inputs;

/* Contains all the dynamical variable arrays. */
typedef struct {
    double *watpos;
    double *watvel;
    double *watacc;
    double *partpos;
    double *partvel;
    double *partacc;
} Dyn_Vars;

#endif /* UTILS_H */