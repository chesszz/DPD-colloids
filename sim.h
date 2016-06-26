#ifndef SIM_H
#define SIM_H

#include "utils.h"

/* Evolves the system from the given initial conditions and input parameters.
 * Will print the variables' time evolution into a file if PRINT_STEPS is 
 * defined to be 1 in utils.h.
 */ 
void evolve_system(Dyn_Vars *dyn_vars, Inputs in);

#endif /* SIM_H */
