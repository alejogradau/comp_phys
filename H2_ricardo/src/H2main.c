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
  double alpha;
  double alpha_range = 0.24-0.06;
  int n_alphas = 20;  // Without considering the first alpha --> one additional
  int n_runs = 1;  // Number of independent runs
  double dalpha = alpha_range/(n_alphas);
  unsigned int N = 1e6;
  unsigned int n_accepted;
  int n_bins = 20;
  double bin_size;  // Returned by map_to_int() function
  double burn_factor = 0.01;  // With extreme initial cond. Neq = 100
  double burn_period = burn_factor*N;  //"Burn-in" is burn_factor% of total run
  unsigned int start = round(burn_period);
  unsigned int n_production = N-start;  // Steps in production run
  double d = 1.0;  // symmetric displacement generated in gen_trial_change
  double Z = 27/16;

  // Memory allocation
  double *conf_m = calloc(6, sizeof(double));  // Configuration m coordinates
  double *pos = calloc(n_production, sizeof(double)); // Electron1's rad position.

  // File names
  //char fhistogram[20] = "./out/histogram.csv";


//  conf_m[0] = 15;
//  conf_m[1] = 0;
//  conf_m[2] = 0;
//  conf_m[3] = -15;
//  conf_m[4] = 0;
//  conf_m[5] = 0;
    

    for(int i = 0; i < n_alphas; i++){
        alpha = 0.1 + (i * dalpha);
        for(int run = 0; run < n_runs; run++){
            // Generate random initial configurations within the range (-2.0, 2.0)
            conf_m[0] = (rand_num-0.5)*4.0;
            conf_m[1] = (rand_num-0.5)*4.0;
            conf_m[2] = (rand_num-0.5)*4.0;
            conf_m[3] = (rand_num-0.5)*4.0;
            conf_m[4] = (rand_num-0.5)*4.0;
            conf_m[5] = (rand_num-0.5)*4.0;
            printf("Electron 1 initial coordinates (%f,%f,%f)\n",
                    conf_m[0], conf_m[1], conf_m[2]);
            printf("Electron 2 initial coordinates (%f,%f,%f)\n",
                    conf_m[3], conf_m[4], conf_m[5]);
            
            n_accepted = mc_integration_metropolis(N, alpha, burn_factor, d, conf_m, pos, run);

            printf("--------------------------------\n\n");
        }
    }
  //printf("Metropolis done\n");
  //printf("n_accepted = %u\n", n_accepted);

//  int *int_pos = calloc(n_production, sizeof(int)); // E1's rad int position.
  //printf("Memory allocated\n");
  //printf("Slicing array\n");
  //slice_array(N, pos_large, n_accepted, pos);
  //free(pos_large);
  //printf("Array sliced\n");

  //file_to_array("./out/pos_large.csv", n_accepted, pos); // DONT DO THIS
  //printf("file_to_array OK\n");


  //printf("array_to_file ...\n");
  //array_to_file("./out/pos.csv", n_accepted, pos);
  //printf("array_to_file OK\n");

//  bin_size = map_to_int(n_production, pos, int_pos, n_bins);
  //printf("Array mapped to integers\n");
//
//  build_histogram(n_accepted, int_pos);
//  radial_density_file("./out/radial_density.csv", bin_size+1, n_bins, Z);
}
