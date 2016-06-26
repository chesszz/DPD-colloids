#ifndef INIT_H
#define INIT_H

/* Reads inputs from stdin that correspond to the vaious parameters. Returns
 * an Inputs struct that stores those parameters.
 */
Inputs get_inputs(void);

/* Allocates space for the dynamical variables - position, velocity, 
 * acceleration. Fills the positions randomly within the box, and then fills
 * the velocities with a Maxwell-Boltzmann distribution. Returns a Dyn_Vars
 * pointer that stores these arrays.
 */
Dyn_Vars *initialise(Inputs in);

#endif /* INIT_H */