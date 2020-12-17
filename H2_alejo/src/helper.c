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

double gen_trial_change(double accepted, double rn, double d)
{
    //printf("The random number is: %f\n", rn);
    return accepted + (d*(rn-0.5));
}

/* Calculate the magnitude of a vector, given its cartesian coordinates */
double vector_magnitude(double x, double y, double z)
{
  double x_sq = pow(x, 2.0);
  double y_sq = pow(y, 2.0);
  double z_sq = pow(z, 2.0);

  double r_sq = x_sq + y_sq + z_sq;
  double r = sqrt(r_sq);

  return r;
}

/* Calculate the relative probability (q) of the trial configuration (t)
 * relative to the current configuration (m) */
double relative_prob(double x1_m, double y1_m, double z1_m,
                     double x2_m, double y2_m, double z2_m,
                     double x1_t, double y1_t, double z1_t,
                     double x2_t, double y2_t, double z2_t,
                     double alpha)
{

  double r1_t = vector_magnitude(x1_t, y1_t, z1_t);
  double r1_m = vector_magnitude(x1_m, y1_m, z1_m);
  double r2_t = vector_magnitude(x2_t, y2_t, z2_t);
  double r2_m = vector_magnitude(x2_m, y2_m, z2_m);

  // r12_t cartesian coordinates configuration (t)
  double x12_t = x1_t - x2_t;
  double y12_t = y1_t - y2_t;
  double z12_t = z1_t - z2_t;

  // r12_m cartesian coordinates configuration (m)
  double x12_m = x1_m - x2_m;
  double y12_m = y1_m - y2_m;
  double z12_m = z1_m - z2_m;

  double r12_t = vector_magnitude(x12_t, y12_t, z12_t);
  double r12_m = vector_magnitude(x12_m, y12_m, z12_m);

  double r1_tm = r1_t - r1_m;
  double r2_tm = r2_t - r2_m;
  //double r12_tm = r12_t - r12_m;

  double numerator_t = 1 + (alpha*r12_t);
  double numerator_m = 1 + (alpha*r12_m);

  /* CORRELATED ELECTRONS */
  double q = exp(-4*(r1_tm + r2_tm)) *
              exp((r12_t/numerator_t) - (r12_m/numerator_m));

  /* UNCORRELATED ELECTRONS */
  //double q = exp(-4*(r1_tm + r2_tm));
  return q;
}

double local_energy(double x1, double x2,
  double y1, double y2,
  double z1, double z2, double alpha)
{

    double r1 = vector_magnitude(x1, y1, z1);
    double r2 = vector_magnitude(x2, y2, z2);

    // Unit vector components
    double sx1 = x1 / r1;
    double sx2 = x2 / r2;
    double sy1 = y1 / r1;
    double sy2 = y2 / r2;
    double sz1 = z1 / r1;
    double sz2 = z2 / r2;

    // r12 cartesian coordinates
    double x12 = x1 - x2;
    double y12 = y1 - y2;
    double z12 = z1 - z2;

    // r12_unitary cartesian coordinates
    double sx12 = sx1 - sx2;
    double sy12 = sy1 - sy2;
    double sz12 = sz1 - sz2;

    // Calculate the dot product on the second term of the local energy (eq.7)
    double r12_dot = (sx12*x12) + (sy12*y12) + (sz12*z12);

    double r12 = vector_magnitude(x12, y12, z12);

    double numerator = 1 + (alpha*r12);
    double numerator_sq = pow(numerator, 2.0);
    double numerator_cub = pow(numerator, 3.0);
    double numerator_tetra = pow(numerator, 4.0);

    double energy = -4 + (r12_dot/(r12*numerator_sq)) - (1/(r12*numerator_cub))
            - (1/(4*numerator_tetra)) + (1/r12);

    return energy;
}

/* Performs Monte Carlo integration using the Metropolis algorithm
 * @ N - Number of time steps
 * @ alpha - Trial wave function parameter
 * @ conf_m - Coordinates of the configuration m (2 electrons)
 * @ pos - Array of electron 1's radial position.
 */
 unsigned int mc_integration_metropolis(unsigned int N,
              double alpha, double burn_factor, double d,
              double *conf_m, double *pos, int run)
{
    unsigned int n_accepted = 0;
    double acceptance_ratio;
    unsigned int n_production;  // Steps in production run
    double burn_period = burn_factor*N;  //"Burn-in" is burn_factor% of total run
    unsigned int start = round(burn_period);
    double E_i = 0;
    double E_mean = 0;
    double E2_mean = 0;
    double sigma_E;
    double sigma_n;

    char fname1[80];
    sprintf(fname1, "./out/local_energy_alpha%.2f_run%d.csv", alpha, run);
    //FILE *fp = fopen(fname1, "w");
    //fprintf(fp, "time, local energy\n");
    //FILE *gp = fopen("./out/pos.csv", "w");
    //fprintf(gp, "time, position\n");
    FILE *hp = fopen("./out/configurations.csv", "w");
    fprintf(hp, "x1_m, y1_m, z1_m, x2_m, y2_m, z2_m\n");

    // Declare Variables
    double x1_m, y1_m, z1_m, x2_m, y2_m, z2_m;  // current configurations (m)
    double x1_t, y1_t, z1_t, x2_t, y2_t, z2_t;  // trial configurations (t)

    // Electron 1 initial positions
    x1_m = conf_m[0];
    y1_m = conf_m[1];
    z1_m = conf_m[2];
    // Electron 2 initial positions
    x2_m = conf_m[3];
    y2_m = conf_m[4];
    z2_m = conf_m[5];

    srand(time(NULL));
    for (unsigned long i = 0; i < N; i++)
    {
        //Generates trial changes based on the current accepted values
        // Electron 1
        x1_t = gen_trial_change(x1_m, rand_num, d);
        y1_t = gen_trial_change(y1_m, rand_num, d);
        z1_t = gen_trial_change(z1_m, rand_num, d);

        // Electron 2
        x2_t = gen_trial_change(x2_m, rand_num, d);
        y2_t = gen_trial_change(y2_m, rand_num, d);
        z2_t = gen_trial_change(z2_m, rand_num, d);

        //Calculates the relative probability q = p_t/p_m
        double q = relative_prob(x1_m, y1_m, z1_m, x2_m, y2_m, z2_m,
                                 x1_t, y1_t, z1_t, x2_t, y2_t, z2_t,
                                 alpha);
        printf("The relative probability is %f\n", q);

        /* Decide if trial change is accepted based on q
         * If not, the configuration is not updated */
        double r = rand_num;
        printf("The random number is: %f\n", r);
        printf("If %f > %f, update configuration.\n", q, r);
        if (q >= r)
        {
            // Electron 1
            x1_m = x1_t;
            y1_m = y1_t;
            z1_m = z1_t;
            printf("Electron 1 updated coordinates (%f,%f,%f)\n",
                    x1_m, y1_m, z1_m);

            // Electron 2
            x2_m = x2_t;
            y2_m = y2_t;
            z2_m = z2_t;
            printf("Electron 2 updated coordinates (%f,%f,%f)\n",
                    x2_m, y2_m, z2_m);
            n_accepted++;
        }

        if (i >= burn_period)
        {
          unsigned int indx = i-start;
          printf("Burn period is over, production run:\n");
          pos[indx] = vector_magnitude(x1_m,y1_m,z1_m);
          printf("pos[%u] = %f\n", indx, pos[indx]);
          E_i = local_energy(x1_m, x2_m, y1_m, y2_m, z1_m, z2_m, alpha);
          //fprintf(fp, "%u, %f\n", indx, E_i);
          //fprintf(gp, "%d, %f\n", indx, pos[indx]);
          fprintf(hp, "%f, %f, %f, %f, %f, %f\n",
                      x1_m, y1_m, z1_m, x2_m, y2_m, z2_m);
          E_mean += E_i;
          E2_mean += pow(E_i, 2.0);
        }

    }
    //fclose(fp);
    //fclose(gp);
    fclose(hp);

    acceptance_ratio = n_accepted*100/N;
    n_production = N-start;
    E_mean /= n_production;
    E2_mean /= n_production;
    sigma_E = sqrt(E2_mean - pow(E_mean, 2.0));
    sigma_n = sigma_E/sqrt(n_production);

    printf("Metropolis integration for N=%u\n", N);
    printf("Symmetric displacement parameter d=%f\n", d);
    printf("Accepted steps (including burn-in period): %d\n", n_accepted);
    printf("Acceptance-rejection ratio:                %f\n", acceptance_ratio);
    printf("E expectation value:                       %f\n", E_mean);
    printf("sigma_E:                                   %f\n", sigma_E);
    printf("sigma_n:                                   %f\n\n", sigma_n);

    return n_accepted;
}

/* Calculates the gradient of the logarithm of the trial wave function with
 * respect to the parameter alpha.
 */
double nabla_wavefunction(double alpha, double r12)
{
  double r12_inv = 1/r12;
  double numerator = r12_inv + alpha;
  double numerator_sq = pow(numerator, 2.0);
  double nabla = - 1/numerator_sq;
  return nabla;
}


/* Performs variational Monte Carlo using the Metropolis algorithm
 * and steepest descent.
 * @ N - Number of time steps
 * @ alpha - Trial wave function parameter
 * @ burn_factor - Fraction of steps to be discarded as "burn-in"
 * @ d - Symmetric displacement parameter
 * @ pos - Array of electron 1's radial position.
 * @ n_p - number of iterations of the damped steepest descent
 * @ beta - exponent of the steepest descent scaling factor
 */
 void variational_mc(unsigned int N, double alpha, double burn_factor,
                     double d, int n_p, double beta)
{
    int p;
    double gamma_p;
    double acceptance_ratio;
    unsigned int n_production;  // Steps in production run
    double burn_period = burn_factor*N;  //"Burn-in" is burn_factor% of total run
    unsigned int start = round(burn_period);
    double r1;
    double r2;
    double r12;
    double nabla_wave_i;
    double nabla_E_p;
    double sigma_E;
    double sigma_n;

    //char fname1[80];
    //sprintf(fname1, "./out/local_energy_alpha%.2f_run%d.csv", alpha, run);
    char fname2[80];
    sprintf(fname2, "./out/steepest_descent_np%d_beta%.2f.csv", n_p, beta);
    //FILE *fp = fopen(fname1, "w");
    //fprintf(fp, "time, local energy\n");
    //FILE *hp = fopen("./out/configurations.csv", "w");
    //fprintf(hp, "x1_m, y1_m, z1_m, x2_m, y2_m, z2_m\n");
    FILE *jp = fopen(fname2, "w");
    fprintf(jp, "p, alpha_p, E, sigma_n\n");

    // Declare Variables
    double x1_m, y1_m, z1_m, x2_m, y2_m, z2_m;  // current configurations (m)
    double x1_t, y1_t, z1_t, x2_t, y2_t, z2_t;  // trial configurations (t)

    for(p = 1; p < n_p; p++)
    {
      unsigned int n_accepted = 0;
      double E_i = 0;
      double E_mean = 0;
      double nabla_wave_mean = 0;
      double E_nabla_wave_mean = 0;
      double E2_mean = 0;
      
      // Electron 1 initial positions
      x1_m = (rand_num-0.5)*4.0;
      y1_m = (rand_num-0.5)*4.0;
      z1_m = (rand_num-0.5)*4.0;
      // Electron 2 initial positions
      x2_m = (rand_num-0.5)*4.0;
      y2_m = (rand_num-0.5)*4.0;
      z2_m = (rand_num-0.5)*4.0;

      srand(time(NULL));
      for (unsigned long i = 0; i < N; i++)
      {
          //Generates trial changes based on the current accepted values
          // Electron 1
          x1_t = gen_trial_change(x1_m, rand_num, d);
          y1_t = gen_trial_change(y1_m, rand_num, d);
          z1_t = gen_trial_change(z1_m, rand_num, d);

          // Electron 2
          x2_t = gen_trial_change(x2_m, rand_num, d);
          y2_t = gen_trial_change(y2_m, rand_num, d);
          z2_t = gen_trial_change(z2_m, rand_num, d);

          //Calculates the relative probability q = p_t/p_m
          double q = relative_prob(x1_m, y1_m, z1_m, x2_m, y2_m, z2_m,
                                   x1_t, y1_t, z1_t, x2_t, y2_t, z2_t,
                                   alpha);
          //printf("The relative probability is %f\n", q);

          /* Decide if trial change is accepted based on q
           * If not, the configuration is not updated */
          double r = rand_num;
          //printf("The random number is: %f\n", r);
          //printf("If %f > %f, update configuration.\n", q, r);
          if (q >= r)
          {
              // Electron 1
              x1_m = x1_t;
              y1_m = y1_t;
              z1_m = z1_t;
              //printf("Electron 1 updated coordinates (%f,%f,%f)\n",
              //        x1_m, y1_m, z1_m);

              // Electron 2
              x2_m = x2_t;
              y2_m = y2_t;
              z2_m = z2_t;
              //printf("Electron 2 updated coordinates (%f,%f,%f)\n",
              //        x2_m, y2_m, z2_m);
              n_accepted++;
          }

          if (i >= burn_period)
          {
            unsigned int indx = i-start;
            //printf("Burn period is over, production run:\n");
            r1 = vector_magnitude(x1_m,y1_m,z1_m);
            r2 = vector_magnitude(x2_m,y2_m,z2_m);
            r12 = r1 - r2;
            E_i = local_energy(x1_m, x2_m, y1_m, y2_m, z1_m, z2_m, alpha);
            nabla_wave_i = nabla_wavefunction(alpha, r12);
            //fprintf(fp, "%u, %f\n", indx, E_i);
            //fprintf(hp, "%f, %f, %f, %f, %f, %f\n",
            //            x1_m, y1_m, z1_m, x2_m, y2_m, z2_m);
            E_mean += E_i;
            nabla_wave_mean += nabla_wave_i;
            E_nabla_wave_mean += (E_i * nabla_wave_i);
            E2_mean += pow(E_i, 2.0);
          }

      }
      //fclose(fp);
      //fclose(hp);

      acceptance_ratio = n_accepted*100/N;
      n_production = N-start;

      E_mean /= n_production;
      nabla_wave_mean /= n_production;
      E_nabla_wave_mean /= n_production;
      nabla_E_p = 2 * (E_nabla_wave_mean - (E_mean * nabla_wave_mean));
      gamma_p = pow(p, -beta);

      E2_mean /= n_production;
      sigma_E = sqrt(E2_mean - pow(E_mean, 2.0));
      sigma_n = sigma_E/sqrt(n_production);

      fprintf(jp, "%d, %f, %f, %f\n", p, alpha, E_mean, sigma_n);
      alpha -= gamma_p * nabla_E_p;

      printf("Metropolis integration for N =             %u\n", N);
      printf("Symmetric displacement parameter d =       %f\n", d);
      printf("Accepted steps (including burn-in period): %d\n", n_accepted);
      printf("Acceptance-rejection ratio:                %f\n", acceptance_ratio);
      printf("E expectation value:                       %f\n", E_mean);
      printf("sigma_E:                                   %f\n", sigma_E);
      printf("sigma_n:                                   %f\n", sigma_n);
      printf("Steepest descent iteration p =             %d\n", p);
      printf("gamma in iteration p =                     %f\n", gamma_p);
      printf("nabla_E in iteration p =                     %f\n", nabla_E_p);
      printf("Alpha in iteration p =                     %f\n\n\n", alpha);
    }
}

void calc_autocorrelation(unsigned long length, double time_series[length], double autocorrelation[length/2])
{
    long double average_f = 0, average_f2 = 0, average_fifik = 0, var_f, autocsum = 0, k_max = 10000;
    FILE *fp = fopen("./out/autocorrelation.csv", "w");

    printf("Calculating expectation values\n");
    for (unsigned long i = 0; i < length; i++) {
        average_f += time_series[i];
        average_f2 += pow(time_series[i], 2.0);
    }
    average_f /= length;
    average_f2 /= length;

    var_f = average_f2 - pow(average_f, 2.0);

    printf("Calculating autocorrelations\n");
    autocorrelation[0] = (average_f2 - pow(average_f, 2.0))/var_f; //Avoids going through an N-sized loop

    for (unsigned long k = 1; k < k_max; k++){
        for (unsigned long i = 0; i < length-k; i++){
            average_fifik += (time_series[i]-average_f) * (time_series[i + k]-average_f);
        }
        average_fifik /= (length-k);
        autocorrelation[k] = (average_fifik)/var_f;
        autocsum += 2*autocorrelation[k];
    }
    printf("s: %Lf\n", autocsum);
    printf("phi_k=s: %f\n", autocorrelation[(int)(autocsum)]);

    printf("Printing results to file\n");
    fprintf(fp, "time, autocorrelation \n");
    for (unsigned long k = 0; k < k_max; k++){
        fprintf(fp, "%ld, %f\n", k, autocorrelation[k]);
    }
    fclose(fp);
}

void file_to_array(char *fname, unsigned int length, double data[length])
{
    FILE *fp;
    fp = fopen(fname, "r");

    for (int i = 0; i < length; i++)
    {
        fscanf(fp, "%lf", &data[i]);
    }
    printf("Reading file with: %d lines\n", length);
    fclose(fp);
}

void calc_block_average(unsigned long length, double time_series[length], double block_average[length/2])
{
    unsigned long num_blocks, block_index;
    double average_f = 0, average_f2 = 0, var_f, bloc_avg, bloc_size_avg, bloc_size_avg_2, var_F, s;
    FILE *fp = fopen("./out/block_averaging.csv", "w");
    fprintf(fp, "time, s \n");

    for (unsigned long i = 0; i < length; i++) {
        average_f += time_series[i];
        average_f2 += pow(time_series[i], 2.0);
    }
    average_f /= length;
    average_f2 /= length;

    var_f = average_f2 - pow(average_f, 2.0);

    for (int block_size = length; block_size > 0; block_size--) {
        s = 0;
        if(length%block_size == 0){
            num_blocks = length/block_size;
            bloc_size_avg = 0;
            bloc_size_avg_2 = 0;

            for (int block = 0; block < num_blocks; block++) {
                bloc_avg = 0;

                for (int block_elem = 0; block_elem < block_size; block_elem++)
                {
                    block_index = (block*block_size)+block_elem;
                    bloc_avg += time_series[block_index];
                }
                bloc_avg /= block_size;
                bloc_size_avg += bloc_avg;
                bloc_size_avg_2 += pow(bloc_avg, 2.0);
            }
            bloc_size_avg /= num_blocks;
            bloc_size_avg_2 /= num_blocks;
            var_F = bloc_size_avg_2 - pow(bloc_size_avg, 2.0);
            s = (block_size*var_F)/var_f;
            fprintf(fp, "%d, %f\n", block_size, s);
        }
    }
    fclose(fp);
}

/* Copies an slice of an array to a shorter array, to throw away empty elements
 * @ N - 2 times length of pos_large arrays
 * @ n_accepted - number of accepted configurations on the Metropolis alg.
 * @ pos_large[] - positions array, TOO LARGE
 * @ pos[] - positions array, shorter copy made
 */
 void slice_array(unsigned int N, double pos_large[],
                  unsigned int n_accepted, double *pos)
{
  //FILE *fp = fopen("./out/radial.csv", "w");
  //fprintf(fp, "time, radial\n");
  for (int i = 0; i < n_accepted; i++)
  {
    pos[i] = pos_large[i];
    //fprintf(fp, "%d, %f, %f\n", i, pos[i], pos_large[i]);
  }
  //fclose(fp);
}

/* Write an array (2nd column) and its index+1 (1st column) to file
 * @ fname - File's name
 * @ length - length of the array
 * @ *array - array of integer values
 */
void int_array_to_file(char *fname, int length, int *array)
{
  //printf("opening...\n");
  FILE *fp = fopen(fname, "w");
  //printf("opened...\n");
  fprintf(fp, "bin, frequency\n");
  //printf("header...\n");
  for (int i = 0; i < length; i++)
  {
    //printf("%d, %d\n", i+1, array[i]);
    fprintf(fp, "%d, %d\n", i+1, array[i]);
  }
  fclose(fp);
}

/* Write an array (2nd column) and its index+1 (1st column) to file
 * @ fname - File's name
 * @ length - length of the array
 * @ *array - array of double values
 */
void array_to_file(char *fname, int length, double *array)
{
  printf("opening...\n");
  FILE *fp = fopen(fname, "w");
  printf("opened...\n");
  fprintf(fp, "bin, frequency\n");
  printf("header...\n");
  for (int i = 0; i < length; i++)
  {
    printf("%d, %f\n", i+1, array[i]);
    fprintf(fp, "%d, %f\n", i+1, array[i]);
  }
  fclose(fp);
}

/* Maps a positions array (pos) to an integer positions array (int_pos), to be
 * used to build a histogram.
 * @ n_accepted - length of pos and int_pos arrays
 * @ pos[] - positions array
 * @ int_pos[] - integer positions arrays
 * @ n_bins - number of bins on histogram to build
 * RETURNS: bin_size (to be used in radial_density_file)
 */
double map_to_int(unsigned int n_accepted, double pos[n_accepted],
                int *int_pos, int n_bins)
{
  double maxval = 0;

  for (int i = 0; i < n_accepted; i++)
  {
    if (pos[i] > maxval)
    {
      maxval = pos[i];
    }
  }

  //printf("maxval = %f\n", maxval);

  double bin_size = maxval/n_bins;

  for (int i = 0; i < n_accepted; i++)
  {
    int_pos[i] = round((pos[i]/bin_size) + 0.5);
    //printf("%f --> %d\n", pos[i], int_pos[i]);
  }
  return bin_size;
}

/* Build a histogram.
 * @ numvalues - Same value as n_accepted: length of pos and int_pos arrays
 * @ values[] - Same as pos[]: positions array
 * SOURCE: https://www.techwalla.com/articles/how-to-create-a-histogram-using-c-programming-code
 */
void build_histogram(unsigned int numvalues, int *values)
{
  int i = 0;

  // Set maxval
  int maxval = 0;
  for (i=0; i<numvalues; i++)
  {
    if (values[i] > maxval)
    {
      maxval = values[i];
    }
  }

  // Set minval
  int minval = maxval;
  for (i=0; i<numvalues; i++)
  {
    if (values[i] < minval)
    {
      minval = values[i];
    }
  }

  // Size of frequency array
  int freqsize = maxval - minval + 1;

  // Declare and initialize frequency array
  int frequency[freqsize];
  for (i=0; i<freqsize; i++)
  {
    frequency[i] = 0;
  }

  // Count the frequency of each integer
  for (i = 0 ; i < numvalues ; i++)
  {
    int index = values[i] - minval;
    frequency[index]++;
  }
  //printf("minval and maxval: (%d, %d)\n", minval, maxval);
  int_array_to_file("./out/histogram.csv", freqsize, frequency);
}

/* Writes a file of radial probability density with the corresponding integer
 * positions to compare with the sampled histogram.
 * @ bin_size - Size of bin, returned by map_to_int
 * @ n_bins - number of bins (for some reason the real val is n_bins+1)
 * @ Z - Number of electrons
 */
void radial_density_file(char *fname, double bin_size, int n_bins, double Z)
{
  double rad = 0;
  double rad_sq;
  double density = 0;
  double Z_cb = pow(Z,3.0);

  FILE *fp = fopen(fname, "w");
  fprintf(fp, "bin, radial density\n");

  for (int i = 0; i < n_bins+1; i++)
  {
    rad = bin_size * (i + 0.5);
    rad_sq = pow(rad, 2.0);
    density = Z_cb * 4 * rad_sq * exp(-2 * Z * rad);
    fprintf(fp, "%d, %f\n", i+1, density);
  }
  fclose(fp);
}