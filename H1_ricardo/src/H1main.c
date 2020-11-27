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
    if(argc != 4){
        printf("Incorrect number of parameters, exiting\n");
        return(1);
    }
    
    //Read values of parameters and calculate number of needed timesteps
    unsigned int total_time = strtoul(argv[1], NULL, 10);
    double dt = atof(argv[2]);             //in picoseconds, e.g. 0.001 = 1 femtosecond
    unsigned int enable_scaling = strtoul(argv[3], NULL, 10);
    unsigned int n_timesteps = total_time/dt;
    
    //Clear screen before printing results
    system("clear");
    printf("Total Time:            %d ps\n", total_time);
    printf("Time step size:        %f\n", dt);
    printf("Number of time steps:  %d\n", n_timesteps);
    printf("Scaling:               %d\n", enable_scaling);
    
    unsigned int Nc = 4;
    unsigned int N = 4*Nc*Nc*Nc; // 4 total atoms in an FCC unit cell
    double pos[N][3];
    double v_0[N][3];
    double m[N];
    double a0 = 4.030283615073347; //Units?
    
    printf("Initializing FCC lattice coordinates\n");
    init_fcc(pos, Nc, a0);
    deviate_fcc(pos, N, a0);
    
    double L = N*a0;
//    double energy;
    double time_array[n_timesteps];
    double *T = calloc(n_timesteps+1, sizeof(double));
    double *V = calloc(n_timesteps+1, sizeof(double));
    double *E = calloc(n_timesteps+1, sizeof(double));
    double *Temp = calloc(n_timesteps+1, sizeof(double));
    double *Pressure = calloc(n_timesteps+1, sizeof(double));

    /* Initial conditions */
    /* Displacements in Ångstroms */
    arange(time_array, 0, n_timesteps, dt);
    
    T[0] = 0;
    V[0] = get_energy_AL(pos, L, N);
    E[0] = T[0] + V[0];
    
    for(int i = 0; i < N; i++){
        m[i] = 27/12*1.244e-3; //Cross multiplication from mass of Carbon
        v_0[i][0] = 0.0;
        v_0[i][1] = 0.0;
        v_0[i][2] = 0.0;
    }
    
    printf("Simulating Time Evolution for the Kinetic, Potential, and Total Energy\n");
    lattice_velocity_verlet_scaled(n_timesteps, L, N, m, v_0, pos, T, V, E, dt, enable_scaling);

    printf("Calculating Temperature at every time step\n");
    double kb = 8.617333262145e-5;
    double factor = 2/(3*kb);
    for(int i = 0; i < n_timesteps; i++){
        Temp[i] = (T[i]/N) * factor;
    }
    
//    printf("Calculating Pressure at every time step\n");
//    double num_cells = Nc*Nc*Nc;
//    double volume = num_cells*pow(a0, 3.0);
//    for(int i = 0; i < n_timesteps; i++){
//        Pressure[i] = 1.0/(3.0*volume) * 2.0 *(T[i]-V[i]);
//        Pressure[i] /= 624e-7; //Conversion to bars
//    }
    
    
    printf("Writing Results to Disk\n");
    write_energies_file("./output/energy.csv", time_array, n_timesteps, T, V, E);
    write_temperatures_file("./output/temperature.csv", time_array, n_timesteps, Temp);
//    write_temperatures_file("./output/pressure.csv", time_array, n_timesteps, Pressure);
    printf("Last average values:\n");
    printf("T: %f\n", Temp[n_timesteps-1]);
//    printf("P: %f\n", Pressure[n_timesteps-1]);
}
