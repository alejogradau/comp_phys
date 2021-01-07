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
    double x0 = 0.0;  // Angstroms; expectation value of position
    double E0 = 0.1;  // eV
    double m = 104.455;  // eV * (fs)^2 / (Å^2)
    double p0 = sqrt(2 * m * E0);  // expectation value of momentum

    double xmax = 3 * d;  // Maximum grid reach is 3 times the width of Gaussian
    int N = 101;  // Number of grid points
    double dx = (2*xmax) / (N-1);  // Space grid steps



    // Memory allocation
    double *x_arr = calloc(N, sizeof(double)); // Space grid
    double *p_arr = calloc(N, sizeof(double)); // Momentum array
    double *gaussian = calloc(N, sizeof(double));  // Gaussian wave packet
    double *prob_density = calloc(N, sizeof(double));  // Prob. density dist.

    // Build space grid
    space_grid(x_arr, N, x0, dx);
    array_to_file("./out/grid.csv", N, x_arr, "index, grid point\n");

    // Build gaussian wave packet and probability density distribution
    gaussian_packet(gaussian, x_arr, N, x0, p0, d);
    probability_density(prob_density, gaussian, N);
    array_to_file("./out/prob_density.csv", N, prob_density,
                  "index, probability density\n");



    // File names
    //char fhistogram[20] = "./out/histogram.csv";


}
