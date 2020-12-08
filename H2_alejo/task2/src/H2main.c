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

 /******************************************************************************
  * Helper functions
  *****************************************************************************/
 void mc_integration_metropolis(unsigned int N, double alpha);
 double gen_trial_change(double accepted, double r);
 double vector_magnitude(double x, double y, double z);
 double relative_prob(double x1_m, double y1_m, double z1_m,
   double x2_m, double y2_m, double z2_m,
   double x1_t, double y1_t, double z1_t,
   double x2_t, double y2_t, double z2_t,
   double alpha);
 double local_energy(double x1, double x2,
   double y1, double y2,
   double z1, double z2, double alpha);
 void calc_autocorrelation(unsigned long length, double time_series[length],
   double autocorrelation[length/2]);
 void file_to_array(unsigned int length, double data[length]);
 void calc_block_average(unsigned long length, double time_series[length],
   double block_average[length/2]);





int main(int argc, char *argv[])
{
  // Define parameters
  double alpha = 0.1;
  mc_integration_metropolis(1e5, alpha);
}

void mc_integration_metropolis(unsigned int N, double alpha)
{
    unsigned int n_accepted = 0;
    double burn_period = 0.01*N; //"Burn-in" period set to 1% of total run
    // current configurations (m)
    double x1_m, y1_m, z1_m, x2_m, y2_m, z2_m;
    // trial configurations (t)
    double x1_t, y1_t, z1_t, x2_t, y2_t, z2_t;
    double E_i = 0, E_mean = 0, E2_mean = 0, sigma_E, sigma_n;
    FILE *fp = fopen("./out/local_energy.csv", "w");

    srand(time(NULL));

    // Generate values for initial guesses
    // Electron 1
    x1_m = 20.0;
    y1_m = 20.0;
    z1_m = 20.0;
    printf("Electron 1 initial coordinates (%f,%f,%f)\n", x1_m, y1_m, z1_m);

    // Electron 2
    x2_m = 20.1;
    y2_m = 20.1;
    z2_m = 20.1;
    printf("Electron 1 initial coordinates (%f,%f,%f)\n", x2_m, y2_m, z2_m);

    fprintf(fp, "time, local energy\n");
    for (unsigned long i = 0; i < N; i++)
    {
        //Generates trial changes based on the current accepted values
        // Electron 1
        x1_t = gen_trial_change(x1_m, rand_num);
        y1_t = gen_trial_change(y1_m, rand_num);
        z1_t = gen_trial_change(z1_m, rand_num);

        // Electron 2
        x2_t = gen_trial_change(x2_m, rand_num);
        y2_t = gen_trial_change(y2_m, rand_num);
        z2_t = gen_trial_change(z2_m, rand_num);

        //Calculates the relative probability q = p_t/p_m
        double q = relative_prob(x1_m, y1_m, z1_m, x2_m, y2_m, z2_m,
                          x1_t, y1_t, z1_t, x2_t, y2_t, z2_t,
                          alpha);

        // Decide if trial change is accepted or not based on q
        if (q > rand_num)
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

            E_i = local_energy(x1_m, x2_m, y1_m, y2_m, z1_m, z2_m, alpha);
            fprintf(fp, "%d, %f\n", n_accepted, E_i);
            E_mean += E_i;
            E2_mean += pow(E_i, 2.0);
            n_accepted++;
        }
    }
    fclose(fp);

    E_mean /= n_accepted;
    E2_mean /= n_accepted;
    sigma_E = sqrt(E2_mean - pow(E_mean, 2.0));
    sigma_n = sigma_E/sqrt(n_accepted);

    printf("Metropolis integration for N=%d\n", N);
    printf("Accepted steps (excluding burn-in period): %d\n", n_accepted);
    printf("Acceptance rate:                           %f%%\n", n_accepted/(N-burn_period)*100);
    printf("E expectation value:                       %f\n", E_mean);
    printf("sigma_E:                                   %f\n", sigma_E);
    printf("sigma_n:                                   %f\n\n", sigma_n);
}

double gen_trial_change(double accepted, double rn)
{
    double d = 1.5;
    return accepted + d*(rn-0.5);
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

  double q = exp(-4*(r1_tm + r2_tm)) *
              exp((r12_t/numerator_t) - (r12_m/numerator_m));
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

    double sx12 = sx1 - sx2;
    double sy12 = sy1 - sy2;
    double sz12 = sz1 - sz2;

    // Calculate the dot product on the second term of the local energy (eq.7)
    double r12_dot = sx12*x12 + sy12*y12 + sz12*z12;

    double r12 = vector_magnitude(x12, y12, z12);

    double numerator = 1 + alpha*r12;
    double numerator_sq = pow(numerator, 2.0);
    double numerator_cub = pow(numerator, 3.0);
    double numerator_tetra = pow(numerator_sq, 2.0);

    double energy = -4 + r12_dot/(r12*numerator_sq) - 1/(r12*numerator_cub)
            - 1/(4*numerator_tetra) + 1/r12;

    return energy;
}

void calc_autocorrelation(unsigned long length, double time_series[length], double autocorrelation[length/2]){
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

void file_to_array(unsigned int length, double data[length]){
    FILE *fp;
    fp = fopen("./dat/MC.txt", "r");

    for (int i = 0; i < length; i++) {
        fscanf(fp, "%lf", &data[i]);
    }
    printf("Reading file with: %d lines\n", length);
}

void calc_block_average(unsigned long length, double time_series[length], double block_average[length/2]){
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
