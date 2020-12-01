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
    double temp_inter= 0, temp_eq = 0, p_eq = 0, t_eq = 0;

    if(enable_scaling){
      if(argc < 8){
        printf("Incorrect number of parameters, exiting\n");
        return(1);
      }
      temp_inter = atof(argv[4]);
      temp_eq = atof(argv[5]);
      p_eq = atof(argv[6]);
      t_eq = atof(argv[7]);
    }

    //Clear screen before printing results
    //system("clear");
    printf("Total Time Simulation: %f ps\n", total_time);
    printf("Time step size:        %f ps\n", dt);
    printf("Number of time steps:  %d\n", n_timesteps);
    printf("Scaling:               %d\n", enable_scaling);
    printf("T_intermediate:        %f K\n", temp_inter);
    printf("T_Eq:                  %f K\n", temp_eq);
    printf("P_Eq:                  %f bar\n", p_eq);
    printf("Equilibration Time:    %f ps\n\n", t_eq);

    //Lattice parameters
    const unsigned int Nc = 4;
    const unsigned int N = 4*Nc*Nc*Nc; // 4 total atoms in an FCC unit cell
    const double a0 = 4.030283615073347; //Å

    double pos[N][3];
    double v_0[N][3];
    double m[N];
    double time_array[n_timesteps+1];
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

    arange(time_array, 0, n_timesteps+1, dt);

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
    lattice_velocity_verlet_scaled(n_timesteps, a0, Nc, N, m, v_0, pos,
      T, V, E, dt, enable_scaling,
      temp_inter, temp_eq, p_eq, t_eq, Temp, Pressure, a0_ev, time_array);

    //Calculating time averages for Pressure and Temperature
    //calc_time_average(n_timesteps, Temp, Temp_exp);
    //calc_time_average(n_timesteps, Pressure, Pressure_exp);

    // Calculating averages after equilibration
    double Temp_equilibrium = calc_eq_average(n_timesteps, Temp,
      n_timesteps-t_eq);
    double Pressure_equilibrium = calc_eq_average(n_timesteps, Pressure,
      n_timesteps-t_eq);
    double a0_equilibrium = calc_eq_average(n_timesteps, a0_ev,
      n_timesteps-t_eq);
    double K_equilibrium = calc_eq_average(n_timesteps, T, n_timesteps/2);
    double V_equilibrium = calc_eq_average(n_timesteps, V, n_timesteps/2);

    // Calculating variances after Equilibration
    double K_var = calc_eq_var(n_timesteps, T, K_equilibrium,
      n_timesteps-t_eq);
    double V_var = calc_eq_var(n_timesteps, V, V_equilibrium,
      n_timesteps-t_eq);

    /* Calculating heat capacity using either kinetic
     * or potential energy variance (eq. 57 and 58)
     */

    double cv_K = calc_cv_NVE(N, Temp_equilibrium, K_var);
    double cv_V = calc_cv_NVE(N, Temp_equilibrium, V_var);


    //Shift Potential and Total Energy so E = 0 after equilibration
    const double E_shift = E[n_timesteps];
    for(int i = 0; i < n_timesteps; i++){
        V[i] -= E_shift;
        E[i] -= E_shift;
    }

    //Calculating time averages for Pressure and Temperature
    //calc_time_average(n_timesteps, Temp, Temp_exp);
    //calc_time_average(n_timesteps, Pressure, Pressure_exp);

    printf("Writing Results to Disk\n\n");
    write_energies_file("./out/energy.csv", time_array, n_timesteps, T, V, E);
    write_observable_file("./out/temperature.csv", "temperature", time_array, n_timesteps, Temp);
    write_observable_file("./out/pressure.csv", "pressure", time_array, n_timesteps, Pressure);
    write_observable_file("./out/a0.csv", "a0", time_array, n_timesteps, a0_ev);

    printf("Final average values after equilibration:\n");
    printf("T: %f\n", Temp_equilibrium);
    printf("P: %f\n", Pressure_equilibrium);
    printf("a0: %f\n\n", a0_equilibrium);
    //printf("Kinetic: %f\n", K_equilibrium);
    //printf("Potential: %f\n", V_equilibrium);

    printf("Variance from the mean values after equilibration:\n");
    printf("var Kinetic: %f\n", K_var);
    printf("var Potential: %f\n\n", V_var);

    printf("Heat capacity using either kinetic or potential energy variances:\n");
    printf("Cv_K: %f\n", cv_K);
    printf("Cv_V: %f\n", cv_V);
}
