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
 #include "fft.h"
 #include <complex.h>  // Complex numbers

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

    printf("p0 = %f\n", p0);
    printf("i^2 = %f + i%f\n", creal(I*I), cimag(I*I));
    double xmax = 3 * d;  // Maximum grid reach is 3 times the width of Gaussian
    int N = 101;  // Number of grid points
    double dx = (2*xmax) / (N-1);  // Space grid steps



    // Memory allocation
    double *x_arr = calloc(N, sizeof(double)); // Space grid
    double *p_arr = calloc(N, sizeof(double)); // Momentum array
    double complex *gaussian_x = calloc(N, sizeof(double complex));  // Gaussian wave packet position
    double complex *gaussian_p = calloc(N, sizeof(double complex));  // Gaussian wave packet momentum
    double *prob_density_x = calloc(N, sizeof(double));  // Prob. density dist. position
    double *prob_density_p = calloc(N, sizeof(double));  // Prob. density dist. momentum

    // Build space grid
    space_grid(x_arr, N, x0, dx);
    array_to_file("./out/grid_x.csv", N, x_arr,
                  "index, grid point position space\n");

    // Build gaussian wave packet and probability density distribution
    gaussian_packet(gaussian_x, x_arr, N, x0, p0, d);
    probability_density(prob_density_x, gaussian_x, N);
    array_to_file("./out/prob_density_x.csv", N, prob_density_x,
                  "index, probability density position space\n");

    // Fill in momentum array
    fft_momentum_shift(p_arr, dx, N);
    array_to_file("./out/grid_p.csv", N, p_arr,
                  "index, grid point momentum space\n");

    // Perform fft, and calculate prob. dens. in p space
    powerspectrum(gaussian_x, gaussian_p, N);
    powerspectrum_shift(gaussian_p, N);
    probability_density(prob_density_p, gaussian_p, N);
    array_to_file("./out/prob_density_p.csv", N, prob_density_p,
                  "index, probability density momentum space\n");






    // File names
    //char fhistogram[20] = "./out/histogram.csv";


}
