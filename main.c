#include "init.h"
#include "sim.h"
#include "terminate.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

int main (void) {

    srand(time(0));

    Inputs inputs = get_inputs();
    Dyn_Vars *dyn_vars = initialise(inputs);
    output_state(dyn_vars, inputs);

    clock_t start = clock(), diff;
    evolve_system(dyn_vars, inputs);
    diff = clock() - start;

    output_state(dyn_vars, inputs);
    terminate(dyn_vars);

    double msec = diff * 1000.0 / CLOCKS_PER_SEC;
    printf("Time taken %.3f seconds for %d time steps. %.3f ms per step.", msec/1000, inputs.N_STEPS, msec/inputs.N_STEPS);
    
    return 0;
}

