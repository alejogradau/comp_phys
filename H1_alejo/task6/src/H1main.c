/*
 H1main.c

 Created by Anders Lindman on 2013-10-31.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "H1helper.h"
#include "H1lattice.h"
#include "H1potential.h"
#include "H1verlet.h"
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
    if(argc < 4){
        printf("Incorrect number of parameters, exiting\n");
        return(1);
    }

    //Read values of parameters and calculate number of needed timesteps
    double total_time = atof(argv[1]);
    double dt = atof(argv[2]);             //in picoseconds, e.g. 0.001 = 1 femtosecond
    unsigned int enable_scaling = strtoul(argv[3], NULL, 10);
    unsigned int n_timesteps = total_time/dt;
    double temp_inter = 0, delta_temp = 0, temp_eq = 0, p_eq = 0, t_eq = 0;
    if(enable_scaling){
      if(argc < 9){
        printf("Incorrect number of parameters, exiting\n");
        return(1);
      }
      temp_inter = atof(argv[4]);
      temp_eq = atof(argv[5]);
      delta_temp = atof(argv[6]);
      p_eq = atof(argv[7]);
      t_eq = atof(argv[8]);
    }
    unsigned int steps_equilibration = t_eq/dt;

    //Clear screen before printing results
    //system("clear");
    printf("Total Time Simulation: %f ps\n", total_time);
    printf("Time step size:        %f ps\n", dt);
    printf("Number of time steps:  %d\n", n_timesteps);
    printf("Scaling:               %d\n", enable_scaling);
    printf("Intermediate T:        %f K\n", temp_inter);
    printf("T_Eq:                  %f K\n", temp_eq);
    printf("deltaT:                  %f K\n", delta_temp);
    printf("P_Eq:                  %f bar\n", p_eq);
    printf("Equilibration Time:    %f ps\n\n", t_eq);

    //Lattice parameters
    const unsigned int Nc = 4;
    const unsigned int N = 4*Nc*Nc*Nc; // 4 total atoms in an FCC unit cell
    double a0 = 4.030283615073347; //Å

    double pos[N][3];
    double v_0[N][3];
    double m[N];
    double time_array[n_timesteps+1];

    // Memory allocation for T=700ºC
    double *T_700 = calloc(n_timesteps+1, sizeof(double));
    double *V_700 = calloc(n_timesteps+1, sizeof(double));
    double *E_700 = calloc(n_timesteps+1, sizeof(double));
    double *Temp_700 = calloc(n_timesteps+1, sizeof(double));
    double *Pressure_700 = calloc(n_timesteps+1, sizeof(double));
    double *a0_ev_700 = calloc(n_timesteps+1, sizeof(double));

    // Memory allocation for T=(700+dt)ºC
    double temp_up = temp_eq + delta_temp;
    double *T_700up = calloc(n_timesteps+1, sizeof(double));
    double *V_700up = calloc(n_timesteps+1, sizeof(double));
    double *E_700up = calloc(n_timesteps+1, sizeof(double));
    double *Temp_700up = calloc(n_timesteps+1, sizeof(double));
    double *Pressure_700up = calloc(n_timesteps+1, sizeof(double));
    double *a0_ev_700up = calloc(n_timesteps+1, sizeof(double));

    // Memory allocation for T=(700-dt)ºC
    double temp_down = temp_eq - delta_temp;
    double *T_700down = calloc(n_timesteps+1, sizeof(double));
    double *V_700down = calloc(n_timesteps+1, sizeof(double));
    double *E_700down = calloc(n_timesteps+1, sizeof(double));
    double *Temp_700down = calloc(n_timesteps+1, sizeof(double));
    double *Pressure_700down = calloc(n_timesteps+1, sizeof(double));
    double *a0_ev_700down = calloc(n_timesteps+1, sizeof(double));

    /* Initial conditions */
    /* Displacements in Ångstroms */
    printf("Initializing FCC lattice coordinates\n");
    init_fcc(pos, Nc, a0);
    deviate_fcc(pos, N, a0);

    arange(time_array, 0, n_timesteps+1, dt);

    T_700[0] = 0;
    V_700[0] = get_energy_AL(pos, a0*Nc, N);
    E_700[0] = T_700[0] + V_700[0];

    for(int i = 0; i < N; i++){
        m[i] = 27/12*1.244e-3; //Cross multiplication from mass of Carbon
        v_0[i][0] = 0.0;
        v_0[i][1] = 0.0;
        v_0[i][2] = 0.0;
    }

    printf("Long routine: Simulating Time Evolution for the Kinetic, \n");
    printf("Potential, Total Energy and virial term using Verlet. Scaling \n");
    printf("of velocities and positions are done at each time step.\n");
    // Function changes the value of a0, v_0, and pos, which are called by ref.
    verlet_equilibration(n_timesteps, &a0, Nc, N, m, v_0, pos,
      T_700, V_700, E_700, dt, enable_scaling, temp_inter, temp_eq, p_eq, t_eq,
      Temp_700, Pressure_700, a0_ev_700, time_array);

    T_700up[0] = 0;  // No valuable information here.
    V_700up[0] = get_energy_AL(pos, a0*Nc, N);
    E_700up[0] = T_700up[0] + V_700up[0];
    verlet_thermostat(n_timesteps, a0, Nc, N, m, v_0, pos,
      T_700up, V_700up, E_700up, dt, enable_scaling, temp_inter, temp_up, t_eq,
      Temp_700up, Pressure_700up, a0_ev_700up, time_array);

    T_700down[0] = 0;  // No valuable information here.
    V_700down[0] = get_energy_AL(pos, a0*Nc, N);
    E_700down[0] = T_700down[0] + V_700down[0];
    verlet_thermostat(n_timesteps, a0, Nc, N, m, v_0, pos,
      T_700down, V_700down, E_700down, dt, enable_scaling, temp_inter, temp_down, t_eq,
      Temp_700down, Pressure_700down, a0_ev_700down, time_array);

      /* Algorithm:
       * Equilibrate system at T = 973.17
       * Calculate the total energy E of the system
       * Deviate from T = 973.17 by +dT (small)
       * Keep V and N constant in the mean time
       * Calculate the total energy of the system
       * Deviate from T = 973.17 by -dT (small)
       * Keep V and N constant in the mean time
       * Calculate the total energy of the system
       * Perform a finite difference approximation of the derivative
       *
       * Variables needed:
       * Energy of the system at the equilibration temperatures T+dT and T-dT
       * The REAL equilibration temperatures T+dT and T-dT
       */

    // Calculating averages at T=700ºC
    double Temp_equilibrium_700 = calc_eq_average(n_timesteps, Temp_700,
      steps_equilibration);
    double Pressure_equilibrium_700 = calc_eq_average(n_timesteps, Pressure_700,
      steps_equilibration);
    double E_equilibrium_700 = calc_eq_average(n_timesteps, E_700,
      steps_equilibration);

    // Calculating averages at T=700+dt ºC
    double Temp_equilibrium_700up = calc_eq_average(n_timesteps, Temp_700up,
      steps_equilibration);
    double Pressure_equilibrium_700up = calc_eq_average(n_timesteps, Pressure_700up,
      steps_equilibration);
    double E_equilibrium_700up = calc_eq_average(n_timesteps, E_700up,
      steps_equilibration);

    // Calculating averages at T=700-dt ºC
    double Temp_equilibrium_700down = calc_eq_average(n_timesteps, Temp_700down,
      steps_equilibration);
    double Pressure_equilibrium_700down = calc_eq_average(n_timesteps, Pressure_700down,
      steps_equilibration);
    double E_equilibrium_700down = calc_eq_average(n_timesteps, E_700down,
      steps_equilibration);

    // Calculate finite differences
    double diff_E = E_equilibrium_700up - E_equilibrium_700down;
    double diff_T = Temp_equilibrium_700up - Temp_equilibrium_700down;

    double cv = calc_cv_NV(diff_E, diff_T);


    //Shift Potential and Total Energy so E = 0 after equilibration at T=700ºC
    const double E_shift = E_700[n_timesteps];
    for(int i = 0; i < n_timesteps; i++){
        V_700[i] -= E_shift;
        E_700[i] -= E_shift;
        V_700up[i] -= E_shift;
        E_700up[i] -= E_shift;
        V_700down[i] -= E_shift;
        E_700down[i] -= E_shift;
    }


    printf("Writing Results to Disk\n\n");
    write_energies_file("./out/energy_700.csv",
      time_array, n_timesteps, T_700, V_700, E_700);
    write_observable_file("./out/temperature_700.csv", "temperature_700",
      time_array, n_timesteps, Temp_700);
    write_observable_file("./out/pressure_700.csv", "pressure_700",
      time_array, n_timesteps, Pressure_700);
    write_observable_file("./out/a0_700.csv", "a0_700",
      time_array, n_timesteps, a0_ev_700);

    write_energies_file("./out/energy_700up.csv",
      time_array, n_timesteps, T_700up, V_700up, E_700up);
    write_observable_file("./out/temperature_700up.csv", "temperature_700up",
      time_array, n_timesteps, Temp_700up);
    write_observable_file("./out/pressure_700up.csv", "pressure_700up",
      time_array, n_timesteps, Pressure_700up);
    write_observable_file("./out/a0_700up.csv", "a0_700up",
      time_array, n_timesteps, a0_ev_700up);

    write_energies_file("./out/energy_700down.csv",
      time_array, n_timesteps, T_700down, V_700down, E_700down);
    write_observable_file("./out/temperature_700down.csv", "temperature_700down",
      time_array, n_timesteps, Temp_700down);
    write_observable_file("./out/pressure_700down.csv", "pressure_700down",
      time_array, n_timesteps, Pressure_700down);
    write_observable_file("./out/a0_700down.csv", "a0_700down",
      time_array, n_timesteps, a0_ev_700down);


    printf("Final average values after equilibration T=973.15 K:\n");
    printf("T_700: %f K\n", Temp_equilibrium_700);
    printf("E_700: %f eV\n", E_equilibrium_700);
    printf("P_700: %f bar\n\n", Pressure_equilibrium_700);

    printf("Final average values after equilibration:\n");
    printf("T_700up: %f K\n", Temp_equilibrium_700up);
    printf("E_700up: %f eV\n", E_equilibrium_700up);
    printf("P_700up: %f bar\n\n", Pressure_equilibrium_700up);

    printf("Final average values after equilibration:\n");
    printf("T_700down: %f K\n", Temp_equilibrium_700down);
    printf("E_700down: %f eV\n", E_equilibrium_700down);
    printf("P_700down: %f bar\n\n", Pressure_equilibrium_700down);

    printf("Heat capacity:\n");
    printf("Cv: %f\n", cv);
}
