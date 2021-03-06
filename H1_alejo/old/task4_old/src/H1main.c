/*
 H1main.c

 Created by Anders Lindman on 2013-10-31.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "helper.h"
#include "H1lattice.h"
#include "H1potential.h"
#include "verlet.h"
#include "H1equilibration.h"

/* Main program */
int main(int argc, char *argv[])
{
    /*
     Descriptions of the different functions in the files H1lattice.c and
     H1potential.c are listed below.
    */

    /*
     Function that generates a fcc lattice in units of [Å]. Nc is the number of
     primitive cells in each direction and a0 is the lattice parameter. The
     positions of all the atoms are stored in pos which should be a matrix of the
     size N x 3, where N is the number of atoms. The first, second and third column
     correspond to the x,y and z coordinate respectively.
    */

    //Check the correct number of timesteps were given
    if(argc < 6){
        printf("Incorrect number of parameters, exiting\n");
        return(1);
    }

    //Read values of parameters and calculate number of needed timesteps
    unsigned int total_time = strtoul(argv[1], NULL, 10);
    double dt = atof(argv[2]);             //in picoseconds, e.g. 0.001 = 1 femtosecond
    unsigned int n_timesteps = total_time/dt;
    double inter_temp = atof(argv[3]);
    double t_eq = atof(argv[4]);
    double p_eq = atof(argv[5]);


    //system("clear"); //Clear screen before printing results
    printf("Total Time:            %d ps\n", total_time);
    printf("Time step size:        %f\n", dt);
    printf("Number of time steps:  %d\n", n_timesteps);
    printf("Desired intermediate T:  %f\n", inter_temp);
    printf("Desired equilibration T:  %f\n", t_eq);
    printf("Desired equilibration P:  %f\n", p_eq);

    //Lattice parameters
    const unsigned int Nc = 4;
    const unsigned int N = 4*Nc*Nc*Nc; // 4 total atoms in an FCC unit cell
    const double a0 = 4.030037915453137; //Å

    double pos[N][3];
    double v_0[N][3];
    double m[N];
    double time_array[n_timesteps];

    double *T = calloc(n_timesteps+1, sizeof(double));
    double *V = calloc(n_timesteps+1, sizeof(double));
    double *E = calloc(n_timesteps+1, sizeof(double));
    double *Temp = calloc(n_timesteps+1, sizeof(double));
    double *Pressure = calloc(n_timesteps+1, sizeof(double));
    double *a0_ev = calloc(n_timesteps+1, sizeof(double));

    /* Initial conditions */
    /* Displacements in Ångstroms */
    printf("Initializing FCC lattice coordinates\n");
    init_fcc(pos, Nc, a0);
    deviate_fcc(pos, N, a0);

    arange(time_array, 0, n_timesteps, dt);

    T[0] = 0;
    V[0] = get_energy_AL(pos, a0*Nc, N);
    E[0] = T[0] + V[0];

    for(int i = 0; i < N; i++){
        m[i] = 27/12*1.244e-3; //Cross multiplication from mass of Carbon
        v_0[i][0] = 0.0;
        v_0[i][1] = 0.0;
        v_0[i][2] = 0.0;
    }

    printf("Long routine: Simulating Time Evolution for the Kinetic, \n");
    printf("Potential, Total Energy and virial term using Verlet. Scaling \n");
    printf("of velocities and positions are done at each time step.\n");
    verlet_inter_melting(n_timesteps, a0, Nc, N, m, v_0, pos,
      T, V, E, dt, inter_temp, t_eq, p_eq, Temp, Pressure, a0_ev, time_array);

    //Calculating time averages for Pressure and Temperature
    //calc_time_average(n_timesteps, Temp, Temp_exp);
    //calc_time_average(n_timesteps, Pressure, Pressure_exp);

    // Calculating averages after equilibration for Pressure and Temperature
    double Temp_equilibrium = calc_eq_average(n_timesteps, Temp, n_timesteps/2);
    double Pressure_equilibrium = calc_eq_average(n_timesteps, Pressure, 200);
    double a0_equilibrium = calc_eq_average(n_timesteps, a0_ev, n_timesteps/2);

    printf("Writing Results to Disk\n");
    write_energies_file("./output/energy.csv", time_array, n_timesteps, T, V, E);
    write_temperatures_file("./output/temperature.csv", time_array, n_timesteps, Temp);
    write_temperatures_file("./output/pressure.csv", time_array, n_timesteps, Pressure);
    write_temperatures_file("./output/a0.csv", time_array, n_timesteps, a0_ev);

    printf("Final average values after equilibration:\n");
    printf("T: %f\n", Temp_equilibrium);
    printf("P: %f\n", Pressure_equilibrium);
    printf("a0: %f\n", a0_equilibrium);
}
