/******************************************************************************
 * H2 task2
 ******************************************************************************
 * Routine for Variational Monte Carlo
 *
 *
 */

 /*************************************************************
  * Macro defines
  *************************************************************/
 #define rand_num ((double) rand() / (double) RAND_MAX)

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
  double alpha = 0.1;
  unsigned int N = 1e6;
  unsigned int n_accepted;
  int n_bins = 20;
  double bin_size;  // Returned by map_to_int() function
  double burn_factor = 0.0;
  double burn_period = burn_factor*N;  //"Burn-in" is burn_factor% of total run
  double d = 1.5;  // symmetric displacement generated in gen_trial_change
  double Z = 27/16;

  // Memory allocation
  double *conf_m = calloc(6, sizeof(double));  // Configuration m coordinates
  double *pos_large = calloc(N/2.0, sizeof(double)); // Electron1's rad position. TOO LARGE (N/2)

  // File names
  //char fhistogram[20] = "./out/histogram.csv";

  conf_m[0] = 2000.0;
  conf_m[1] = 2000.0;
  conf_m[2] = 2000.0;
  conf_m[3] = -2000.0;
  conf_m[4] = -2000.0;
  conf_m[5] = -2000.0;
  printf("Electron 1 initial coordinates (%f,%f,%f)\n",
          conf_m[0], conf_m[1], conf_m[2]);
  printf("Electron 2 initial coordinates (%f,%f,%f)\n",
          conf_m[3], conf_m[4], conf_m[5]);

  n_accepted = mc_integration_metropolis(N, alpha, burn_factor, d,
                                         conf_m, pos_large);
  printf("Metropolis done\n");
  printf("n_accepted = %u\n", n_accepted);

  double *pos = calloc(n_accepted, sizeof(double)); // E1's rad position.
  int *int_pos = calloc(n_accepted, sizeof(int)); // E1's rad int position.
  //printf("Memory allocated\n");

  //file_to_array("./out/pos_large.csv", n_accepted, pos);
  //printf("file_to_array OK\n");

  //printf("Slicing array\n");
  slice_array(N, pos_large, n_accepted, pos);
  free(pos_large);
  //printf("Array sliced\n");

  //printf("array_to_file ...\n");
  //array_to_file("./out/pos.csv", n_accepted, pos);
  //printf("array_to_file OK\n");

  bin_size = map_to_int(n_accepted, pos, int_pos, n_bins);
  //printf("Array mapped to integers\n");

  build_histogram(n_accepted, int_pos);
  radial_density_file("./out/radial_density.csv", bin_size, n_bins, Z);

  printf("Metropolis integration for N=%u\n", N);
  printf("Accepted steps (excluding burn-in period): %d\n", n_accepted);
  printf("Acceptance rate:            %f%%\n", n_accepted/(N-burn_period)*100);
  printf("Look above for E expectation value, sigma_E, and sigma_n\n");
}
