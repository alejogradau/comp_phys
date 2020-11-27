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
int main()
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

    unsigned int Nc = 4;
    unsigned int N = 4*Nc*Nc*Nc; // 4 total atoms in an FCC unit cell
    double pos[N][3];
    double v_0[N][3];
    double m[N];
    double a0 = 4.030283615073347; //Units?

    printf("Initializing FCC lattice coordinates\n");
    init_fcc(pos, Nc, a0);
    deviate_fcc(pos, N, a0);

    /*
     Function that calculates the potential energy in units of [eV]. pos should be
     a matrix containing the positions of all the atoms, L is the length of the
     supercell and N is the number of atoms.
    */

    double L = N*a0;
//    double energy;
    int verlet_timesteps = 1; double dt = 0.001;
    double time_array[n_timesteps];
    double *T = calloc(verlet_timesteps+1, sizeof(double));
    double *V = calloc(verlet_timesteps+1, sizeof(double));
    double *E = calloc(verlet_timesteps+1, sizeof(double));

    /* Initial conditions */
    /* Displacements in Ångstroms */
    arange(time_array, 0, n_timesteps, dt);

//    T[0] = get_virial_AL(pos, L, N);
    T[0] = 0;
    V[0] = get_energy_AL(pos, L, N);
    E[0] = T[0] + V[0];

    for(int i = 0; i < N; i++){
        m[i] = 27/12*1.244e-3;
        v_0[i][0] = 0.0;
        v_0[i][1] = 0.0;
        v_0[i][2] = 0.0;
    }

    printf("Simulating Time Evolution for the Kinetic, Potential, and Total Energy\n");
    lattice_velocity_verlet(n_timesteps, L, N, m, v_0, pos, T, V, E, dt);

    printf("Writing Results to Disk\n");
    write_energies_file("./output/energy.csv", time_array, n_timesteps, T, V, E);

//    printf("Calculating potential energy in eV\n");
//    energy = get_energy_AL(pos, L, N);

    /*
     Function that calculates the virial in units of [eV]. pos should be a matrix
     containing the positions of all the atoms, L is the length of the supercell
     and N is the number of atoms.
    */

//    double virial;
    //    printf("Calculating virial in eV\n");
    //    virial = get_virial_AL(pos, L, N)/256;

    double kb = 8.617333262145e-5;
    double temperature = (T[n_timesteps-1]/256) * 2/(3*kb);
    printf("T: %f\n", temperature);

    /*
     Function that calculates the forces on all atoms in units of [eV/Å]. the
     forces are stored in f which should be a matrix of size N x 3, where N is the
     number of atoms and column 1,2 and 3 correspond to the x,y and z component of
     the force resepctively . pos should be a matrix containing the positions of
     all the atoms, L is the length of the supercell and N is the number of atoms.
    */

//    double f[N][3];
//    printf("Calculating forces\n");
//    get_forces_AL(f,pos, L, N);

//    printf("A_0:    %f\n", a0);
//    printf("Energy: %f\n", energy);
//    printf("Virial: %f\n", virial);
}
