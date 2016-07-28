#include "sim.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void fill_gaussian_two_tuple(double *tuple, double stdev);
void initialise_pos_vel(Inputs in, double *pos_list, double *vel_list, int num, double *mass, int var_mass_mode);
void initialise_part_size_mass(Inputs in);

/* Reads inputs from stdin that correspond to the vaious parameters. Returns
 * an Inputs struct that stores those parameters.
 */
Inputs get_inputs(void) {

    Inputs in;

    scanf("%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
        &in.N_WATER, &in.N_PARTICLES,
        &in.N_STEPS, &in.TIME_STEP, &in.BOX_SIZE, &in.SHEAR_RATE,
        &in.DAMP_CONST, &in.SPRING_CONST, 
        &in.R_PARTICLE_AVG, &in.R_PARTICLE_STDEV, 
        &in.E_SC, &in.E_CC);

    in.PART_RADII = (double *) calloc(in.N_PARTICLES, sizeof(double));
    in.PART_MASS = (double *) calloc(in.N_PARTICLES, sizeof(double));

    if (in.PART_RADII  == NULL || in.PART_MASS == NULL) {
        fprintf(stderr, "Error in memory allocation!\n");
        exit(1);
    }

    return in;
}

/* Allocates space for the dynamical variables - position, velocity, 
 * acceleration. Fills the positions randomly within the box, and then fills
 * the velocities with a Maxwell-Boltzmann distribution. Returns a Dyn_Vars
 * pointer that stores these arrays.
 */
Dyn_Vars *initialise(Inputs in) {

    FILE *part_radii = fopen("part_radii.out", "w");

    #if TWO_D
        fprintf(stderr, "TWO DIMENSIONAL\n\n\n");
    #endif

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
    dyn_vars->watomg  = (double *) calloc(3 * in.N_WATER,     sizeof(double));
    dyn_vars->watalp  = (double *) calloc(3 * in.N_WATER,     sizeof(double));

    dyn_vars->partpos = (double *) calloc(3 * in.N_PARTICLES, sizeof(double));
    dyn_vars->partvel = (double *) calloc(3 * in.N_PARTICLES, sizeof(double));
    dyn_vars->partacc = (double *) calloc(3 * in.N_PARTICLES, sizeof(double));
    dyn_vars->partomg = (double *) calloc(3 * in.N_PARTICLES, sizeof(double));
    dyn_vars->partalp = (double *) calloc(3 * in.N_PARTICLES, sizeof(double));

    if (dyn_vars->watpos  == NULL || dyn_vars->watvel  == NULL || 
        dyn_vars->watacc  == NULL || dyn_vars->watomg  == NULL || 
        dyn_vars->watalp  == NULL ||
        dyn_vars->partpos == NULL || dyn_vars->partvel == NULL || 
        dyn_vars->partacc == NULL || dyn_vars->partomg == NULL || 
        dyn_vars->partalp == NULL) {
        fprintf(stderr, "Error in memory allocation!\n");
        exit(1);
    }
    
    initialise_part_size_mass(in);

    for (int i = 0; i < in.N_PARTICLES; i++) {
        fprintf(part_radii, "%f\n", in.PART_RADII[i]);
        fprintf(stderr, "Radius: %f, Mass:%f\n", in.PART_RADII[i], in.PART_MASS[i]);
    }
    fclose(part_radii);
    
    /* Creating a new variaable is needed to pass mass of water to the function 
     * since it takes a double pointer. */
    double M_WATER = 1.0;
    initialise_pos_vel(in, dyn_vars->watpos, dyn_vars->watvel, in.N_WATER, &M_WATER, 0);
    initialise_pos_vel(in, dyn_vars->partpos, dyn_vars->partvel, in.N_PARTICLES, in.PART_MASS, 1);

    // dyn_vars->watpos[0] = -2; // 1
    // dyn_vars->watpos[1] = 2.5;
    // dyn_vars->watvel[0] = 10.0;
    // dyn_vars->watpos[3] = 2; // 2
    // dyn_vars->watpos[4] = 2;
    // dyn_vars->watpos[6] = 0; // 3
    // dyn_vars->watpos[7] = 0;
    // dyn_vars->watpos[9] = 4; // 4
    // dyn_vars->watpos[10] = -0.5;
    // dyn_vars->watvel[9] = -10;

    // dyn_vars->partvel[1] = 2;
    // dyn_vars->partpos[1] = 4;
    
    return dyn_vars;
}

/* Take a tuple with at least 2 slots and fills them with Gaussian values that
 * have sigma = 1. Uses the Box-Muller transform.
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

/* Initialise the position and velocities using a random uniform distribution
 * and  Maxwell-Boltzmann distribution respectively. Requires the mass since
 * the M-B distribution depends on the mass. If var_mass_mode is true, then we 
 * interpret the mass as a list of masses. Otherwise it is 
 * interpreted as a normal pointer o a double.
 */
void initialise_pos_vel(Inputs in, double *pos_list, double *vel_list, int num, 
                        double *mass, int var_mass_mode) {
    for (int i = 0; i < num; i++) {

        int x = 3 * i;
        int y = x + 1;
        int z = y + 1;

        /* Random positions from -BOX_SIZE/2 to BOX_SIZE/2. */
        pos_list[x] = in.BOX_SIZE * ((double) rand() / RAND_MAX - 0.5);
        pos_list[y] = in.BOX_SIZE * ((double) rand() / RAND_MAX - 0.5);
        pos_list[z] = in.BOX_SIZE * ((double) rand() / RAND_MAX - 0.5);

        #if TWO_D
            pos_list[z] = 0; /* TODO: Uncomment for 2D. */
        #endif

        /* http://scicomp.stackexchange.com/questions/19969/how-do-i-generate-maxwell-boltzmann-variates-using-a-uniform-distribution-random */
        /* Maxwell-Boltzmann distribution has all components of velocity
         * randomly drawn from a Gaussian with mean 0 and stdev sqrt(kT/m). 
         */
        double gaussian_tuple[2];
        double mass_i;

        /* Interpret "mass" as a list of particle masses. */
        if (var_mass_mode) {
            mass_i = mass[i];
        }
        /* Interpret mass as a simple pointer to a double. Dereference the pointer. */
        else {
            mass_i = *mass;
        }

        /* Stdev is 1/sqrt(m) since we take kT to be 1 by definition of 
         * our units.
         */
        double vel_stdev = 1.0 / sqrt(mass_i);

        fill_gaussian_two_tuple(gaussian_tuple, vel_stdev);
        vel_list[x] = gaussian_tuple[0];
        vel_list[y] = gaussian_tuple[1];

        /* Add a x-velocity gradient in the y-direction from the shear. */
        /* Commented out because it is not included in the actual velocity, this
         * is instead added from the dr/dt term. */
        //vel_list[x] += in.SHEAR_RATE * pos_list[y];

        fill_gaussian_two_tuple(gaussian_tuple, vel_stdev);
        vel_list[z] = gaussian_tuple[1];

        #if TWO_D
            vel_list[z] = 0; /* TODO: Uncomment for 2D. */
        #endif
    }
}

/* Fills in the particle sizes by sampling from a Gaussian distribution with 
 * mean = R_PARTICLE_AVG and standard deviation = R_PARTICLE_STDEV */
void initialise_part_size_mass(Inputs in) {

    double box_vol = in.BOX_SIZE * in.BOX_SIZE * in.BOX_SIZE;

    double particle_vol = 0.0;
    for (int i = 0; i < in.N_PARTICLES; i++) {
        particle_vol += 4.0 / 3 * M_PI * in.PART_RADII[i] * in.PART_RADII[i] * in.PART_RADII[i];
    }

    /* Number of particles per unit volume, excluding the particle-occupied 
     * volumes. This is the number density and also the mass density since each
     * water has mass 1. We assume the particles hve this same density as 
     * well. */
    double rho_particle = in.N_WATER / (box_vol - particle_vol);
    rho_particle = 3; /*TODO */

    double gaussian_tuple[2];

    /* Fill in the particle sizes with the Gaussian distribution in pairs. */
    for (int i = 0; i < in.N_PARTICLES / 2; i++) {

        /* Ensure that the resulting radii are non-negative. */
        do {
            fill_gaussian_two_tuple(gaussian_tuple, in.R_PARTICLE_STDEV);
        } while (gaussian_tuple[0] <= -in.R_PARTICLE_AVG || gaussian_tuple[1] <= -in.R_PARTICLE_AVG);
        
        in.PART_RADII[2*i  ] = gaussian_tuple[0] + in.R_PARTICLE_AVG;
        in.PART_RADII[2*i+1] = gaussian_tuple[1] + in.R_PARTICLE_AVG;

        in.PART_MASS[2*i]   = rho_particle * 4.0 / 3 * M_PI * 
                              in.PART_RADII[2*i] * in.PART_RADII[2*i] * in.PART_RADII[2*i];
        in.PART_MASS[2*i+1] = rho_particle * 4.0 / 3 * M_PI * 
                              in.PART_RADII[2*i+1] * in.PART_RADII[2*i+1] * in.PART_RADII[2*i+1];
    }

    /* If we have an odd number of particles, we manually fill in the last one
     * by itself. */
    if (in.N_PARTICLES % 2 == 1) {

        /* Ensure that the resulting radius is non-negative. */
        do {
            fill_gaussian_two_tuple(gaussian_tuple, in.R_PARTICLE_STDEV);
        } while (gaussian_tuple[0] <= -in.R_PARTICLE_AVG);

        in.PART_RADII[in.N_PARTICLES-1] = gaussian_tuple[0] + in.R_PARTICLE_AVG;

        in.PART_MASS[in.N_PARTICLES-1] = rho_particle * 4.0 / 3 * M_PI * 
                                          in.PART_RADII[in.N_PARTICLES-1] * 
                                          in.PART_RADII[in.N_PARTICLES-1] * 
                                          in.PART_RADII[in.N_PARTICLES-1];
    }
}