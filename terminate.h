#ifndef TERMINATE_H
#define TERMINATE_H

#include "sim.h"

/* Prints some stats on the screen. */
void output_state(Dyn_Vars *dyn_vars, Inputs in);

/* Frees all the dynamical variables lists and the struct itself. */
void terminate(Dyn_Vars *dyn_vars);

#endif /* TERMINATE_H */