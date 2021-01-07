/******************************************************************************
 * H3
 ******************************************************************************
 * Routine for ...
 *
 *
 */

 /******************************************************************************
  * Macro defines
  *****************************************************************************/
 #define PI 3.14159
 #define hbar 0.6582  // eV * fs

 /******************************************************************************
  * Includes
  *****************************************************************************/
 #include <time.h>    //time
 #include <stdlib.h>  //srand, rand, strtol
 #include <math.h>    //pow, sin, exp
 #include <stdio.h>   //printf
 #include "helper.h"

/******************************************************************************
 * MAIN
 *****************************************************************************/
int main(int argc, char *argv[])
{
    // Initialize Variables
    double d = 0.5;  // Angstroms
    double x0 = 0.0;  // Angstroms
    double xmax = 3 * d;  // Maximum grid reach is 3 times the width of Gaussian
    int N = 101;  // Number of grid points
    double dx = (2*xmax) / (N-1);  // Space grid steps
    double m = 104.455;  // eV * (fs)^2 / (Å^2)



    // Memory allocation
    double *x_arr = calloc(N, sizeof(double)); // Space grid
    double *gaussian = calloc(N, sizeof(double));  // Configuration m coordinates

    space_grid(x_arr, N, x0, dx);
    array_to_file("grid.csv", N, x_arr);

    // File names
    //char fhistogram[20] = "./out/histogram.csv";


}
