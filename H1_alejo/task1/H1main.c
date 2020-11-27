/*
 H1main.c

 Created by Anders Lindman on 2013-10-31.
 */

 /*
 *
 * Compile me as:
 * make
 * clang -o H1main H1lattice.o H1potential.o H1main.o
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "H1lattice.h"
#include "H1potential.h"


/******************************************************************************
 * Helper functions
 *****************************************************************************/
 /*
  * constructs a0 array
  * @array - array to be filled with values
  * @start - start value
  * @n_points - number of stamps in array
  * @step - step between two consecutive values
 */
 void arange(double *array, double start, int n_points, double step){
     for(int i = 0; i < n_points; i++){
 	array[i] = start + i*step;
     }
 }
 /*
  * Writes file
  * @fname - File name
  * @a0 - array of lattice constant values
  * @E_pot - array with potential energy values
  * @n_points - number of points
 */
void write_to_file(char *fname, double *a0,
		   double *E_pot_cell, int n_points)
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "a0, energy\n");
    for(int i = 0; i < n_points; ++i){
	    fprintf(fp, "%f,%f\n", a0[i], E_pot_cell[i]);
    }
    fclose(fp);
}

/* Main program */
int main()
{
    /* Declaration of variables */
    int N = 4;
    int natoms = 4 * pow(N,3);
    int dim = 3;
    int n_points = 20;
    double step = 0.005;
    double E_pot;
    double a0[n_points];  // units = Å
    arange(a0, 4.0, n_points, step);

    /* Dynamic memory allocation of arrays */
    double *E_pot_cell = calloc(n_points, sizeof(double));
    double  (*X)[dim] = malloc(sizeof (double [natoms][dim]));
    double *cell_vol = calloc(n_points, sizeof(double));

    /* Loop for potential energy vs lattice constant relation */
    for(int ii = 0; ii < n_points; ++ii){
      printf("%d\n", ii);
      /* (Over)write matrix X: atoms positions on a fcc lattice. */
      init_fcc(X, N, a0[ii]);
      printf("Writing atoms positions\n");

      /* Calculate potential energy per unit cell */
      E_pot = get_energy_AL(X, N * a0[ii], natoms);
      printf("Calculate potential energy\n");
      E_pot_cell[ii] = E_pot / pow(N,3);
      cell_vol[ii] = pow(a0[ii],3);

      printf("Calculate potential energy per cell\n");
    }

    write_to_file("csv/energy_volume.csv", cell_vol, E_pot_cell, n_points);
    /* Print a0 and potential energy */
    // printf("%f,%f,%f\n", a0, E_pot, E_pot_cell);




    /*
     Code for generating a uniform random number between 0 and 1. srand should only
     be called once.
    */
    /*
     srand(time(NULL));
     double random_value;
     random_value = (double) rand() / (double) RAND_MAX;
    */

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
    /*
     init_fcc(pos, Nc, a0);
    */

    /*
     Function that calculates the potential energy in units of [eV]. pos should be
     a matrix containing the positions of all the atoms, L is the length of the
     supercell and N is the number of atoms.
    */
    /*
     double energy;
     energy = get_energy_AL(pos, L, N);
    */

    /*
     Function that calculates the virial in units of [eV]. pos should be a matrix
     containing the positions of all the atoms, L is the length of the supercell
     and N is the number of atoms.
    */
    /*
     double virial;
     virial = get_virial_AL(pos, L, N);
    */

    /*
     Function that calculates the forces on all atoms in units of [eV/Å]. the
     forces are stored in f which should be a matrix of size N x 3, where N is the
     number of atoms and column 1,2 and 3 correspond to the x,y and z component of
     the force resepctively . pos should be a matrix containing the positions of
     all the atoms, L is the length of the supercell and N is the number of atoms.
    */
    /*
     get_forces_AL(f,pos, L, N);
    */



}
