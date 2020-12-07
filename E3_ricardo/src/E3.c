/******************************************************************************
 * E3
 ******************************************************************************
 * Routine for MC integration
 *
 *
 */

/*************************************************************
 * Macro defines
 *************************************************************/

//C99 does not require it to be defined in math.h anymore
#define M_PI 3.14159265358979323846264338327

//F(x) and P(x) from lecture
//#define F(x) (4/(1+x*x))
//#define P(x) ((4-2*x)/3)

//F(x) and P(x) for task 1 and task 2
#define F(x) (x*(1-x))
#define P(x) ((M_PI/2)*sin(M_PI*x))

#define rand_num (use_gsl ? gsl_rng_uniform(q) : ((double) rand() / (double) RAND_MAX))


/******************************************************************************
 * Includes
 *****************************************************************************/
#include <time.h>    //time
#include <stdlib.h>  //srand, rand, strtol
#include <math.h>    //pow, sin, exp
#include <stdio.h>   //printf
#include <gsl/gsl_rng.h>


void mc_integration(unsigned int N, unsigned int importance_sampling, unsigned int use_gsl);
void mc_integration_metropolis(unsigned int N, unsigned int use_gsl);
double gen_trial_change(double accepted, double r);
double metropolis_weight(double x, double y, double z);
double metropolis_function(double x, double y, double z);
void calc_autocorrelation(unsigned long length, double time_series[length], double autocorrelation[length/2]);
void file_to_array(unsigned int length, double data[length]);
void calc_block_average(unsigned long length, double time_series[length], double block_average[length/2]);

int main(int argc, char *argv[]){
    //Clear screen before printing results
    system("clear");
    
    unsigned int importance_sampling = argc > 1 ? strtol(argv[1], NULL, 10) : 0;
    unsigned int use_gsl = argc > 2 ? strtol(argv[2], NULL, 10) : 0;

    mc_integration(1e2, 1, 1);
    mc_integration(1e3, 1, 1);
    mc_integration(1e4, 1, 1);
    mc_integration(1e5, 1, 1);
    mc_integration(1e5, 1, 1);

    mc_integration_metropolis(1e7, 0);
    mc_integration_metropolis(1e7, 1);
    
    unsigned long length = 1000000;
    double *data = calloc(length, sizeof(double));
    double *autocorrelation = calloc(length/2, sizeof(double));
    
    file_to_array(length, data);
    calc_autocorrelation(length, data, autocorrelation);
    calc_block_average(length, data, autocorrelation);
}

void mc_integration(unsigned int N, unsigned int importance_sampling, unsigned int use_gsl){
    double x, I_n, f_mean = 0, f2_mean = 0, sigma_f, sigma_n;
    const gsl_rng_type *T; /* static info about rngs */
    gsl_rng *q; /* rng instance */ gsl_rng_env_setup(); /* setup the rngs */

    if(use_gsl){
        T = gsl_rng_default; /* specify default rng */ q = gsl_rng_alloc(T); /* allocate default rng */
        gsl_rng_set(q,time(NULL)); /* Initialize rng */
    } else{
        srand(time(NULL));
    }
    
    FILE *fp1 = fopen("./out/generated_xi.csv", "w");
    FILE *fp2 = fopen("./out/generated_p_xi.csv", "w");
    fprintf(fp1, "x_i\n");
    fprintf(fp2, "x_i\n");
    
    for (unsigned int i = 0; i < N; i++) {
        x = (double) rand() / (double) RAND_MAX;
        I_n = importance_sampling ? F(x)/P(x) : F(x);
        f_mean += I_n;
        f2_mean += pow(I_n, 2.0);
        fprintf(fp1, "%f\n", x);
        fprintf(fp2, "%f\n", P(x));
    }
    
    if(use_gsl){gsl_rng_free(q);}
    fclose(fp1);
    fclose(fp2);
    
    f_mean /= N;
    f2_mean /= N;
    sigma_f = sqrt(f2_mean - pow(f_mean, 2.0));
    sigma_n = sigma_f/sqrt(N);
    
    printf("MC integration for N=%d\n", N);
    printf("GSL Random Number Generator:               %s\n", use_gsl ? "True" : "False");
    printf("Importance Sampling:                       %s\n", importance_sampling ? "Enabled" : "Disabled");
    printf("I_N:                                       %f\n", f_mean);
    printf("sigma_f:                                   %f\n", sigma_f);
    printf("sigma_i:                                   %f\n\n", sigma_n);
}

void mc_integration_metropolis(unsigned int N, unsigned int use_gsl){
    unsigned int n_accepted = 0;
    double burn_period = 0.01*N; //"Burn-in" period set to 1% of total run
    double x_m, y_m, z_m, x_n, y_n, z_n;
    double p_m = 0, p_n, I_n = 0, f_mean = 0, f2_mean = 0, sigma_f, sigma_n;
    const gsl_rng_type *T; /* static info about rngs */
    gsl_rng *q; /* rng instance */ gsl_rng_env_setup(); /* setup the rngs */

    if(use_gsl){
        T = gsl_rng_default; /* specify default rng */ q = gsl_rng_alloc(T); /* allocate default rng */
        gsl_rng_set(q,time(NULL)); /* Initialize rng */
    } else{
        srand(time(NULL));
    }
    
    //Generate values for initial guesses
    x_m = rand_num;
    y_m = rand_num;
    z_m = rand_num;
    
    //Calculates the probability of initial guesses
    p_m = metropolis_weight(x_m, y_m, z_m);
    
    for (unsigned int i = 0; i < N; i++) {
        //Generates trial changes based on the current accepted values
        x_n = gen_trial_change(x_m, rand_num);
        y_n = gen_trial_change(y_m, rand_num);
        z_n = gen_trial_change(z_m, rand_num);
        
        //Calculates the probability of proposed values
        p_n = metropolis_weight(x_n, y_n, z_n);
        
        //If probability of proposed step is higher than accepted step or than random number
        if (p_n > p_m || p_n/p_m > rand_num){
            x_m = x_n;
            y_m = y_n;
            z_m = z_n;
            
            //If steps are taken after the "burn-in" period, add them to the integral
            if (i > burn_period){
                I_n = metropolis_function(x_m, y_m, z_m);
                f_mean += I_n;
                f2_mean += pow(I_n, 2.0);
                p_m = p_n;
                n_accepted++;
            }
        }
    }
    
    if(use_gsl){gsl_rng_free(q);}
    
    f_mean /= n_accepted;
    f2_mean /= n_accepted;
    sigma_f = sqrt(f2_mean - pow(f_mean, 2.0));
    sigma_n = sigma_f/sqrt(n_accepted);
    
    printf("Metropolis integration for N=%d\n", N);
    printf("GSL Random Number Generator:               %s\n", use_gsl ? "True" : "False");
    printf("Accepted steps (excluding burn-in period): %d\n", n_accepted);
    printf("Acceptance rate:                           %f%%\n", n_accepted/(N-burn_period)*100);
    printf("I_N:                                       %f\n", f_mean);
    printf("sigma_f:                                   %f\n", sigma_f);
    printf("sigma_i:                                   %f\n\n", sigma_n);
}

double gen_trial_change(double accepted, double r) {
    return accepted + 2*(r-0.5);
}

double metropolis_function(double x, double y, double z){
    double x_2 = pow(x, 2.0);
    double y_2 = pow(y, 2.0);
    double z_2 = pow(z, 2.0);
    
    return (x_2 + x_2*y_2 + x_2*y_2*z_2);
}

double metropolis_weight(double x, double y, double z){
    return pow(M_PI,-1.5)*exp(-(pow(x,2.0)+pow(y,2.0)+pow(z,2.0)));
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

                for (int block_elem = 0; block_elem < block_size; block_elem++) {
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
