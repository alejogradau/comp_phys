/******************************************************************************
 * E4main
 ******************************************************************************
 * Routine that runs the velocity verlet algorithm
 * Use as template to construct your program!
 *
 * Compile me as:
 * clang E4main.c -o  -lgsl -lgslcblas
 * clang E1code4.c -o code4 -lm
 */

 /*************************************************************
  * Macro defines
  *************************************************************/
 //C99 does not require it to be defined in math.h anymore
 #define M_PI 3.14159265358979323846264338327

/******************************************************************************
 * Includes
 *****************************************************************************/
#include <stdio.h>   //fopen, fprintf
#include <stdlib.h>  //calloc
#include <stdint.h>  //uint64_t
#include <math.h>    //pow, fabs function (absolute values double)
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>  // Gaussian distribution

#include "fft.h"  // interface to fft routine

/******************************************************************************
 * Helper functions
 *****************************************************************************/
double calc_acc(double q, double omega);
void arange(double *array, double start, int len_t, double dt);
void write_trajectories_file(char *fname, double *time_array,
			   double *position, double *velocity, double *acceleration,
		     long int n_timesteps);
void BD3(long int n_timesteps, double dt,
         double *q_arr, double *v_arr, double *a_arr,
         double omega, double c0, double vth, double eta,
         gsl_rng * r1, gsl_rng * r2, double sigma);
/* For FFT */
void read_data(char *fname, double *time_array, double *signal);
void write_spectrum_file(char *fname, double *frequencies,
		                     double *spectrum, int n_points);

/******************************************************************************
 * MAIN
 *****************************************************************************/
int main()
{
  // GENERATE THREE DIFFERENT RANDOM NUMBER GENERATOR
  // Instantiate the random number generator
  const gsl_rng_type * T1;
  gsl_rng * r1;
  const gsl_rng_type * T2;
  gsl_rng * r2;

  /* create a generator chosen by the
     environment variable GSL_RNG_TYPE */
  gsl_rng_env_setup();
  T1 = gsl_rng_default;
  r1 = gsl_rng_alloc (T1);
  T2 = gsl_rng_default;
  r2 = gsl_rng_alloc (T2);

  // Calculate mass of silica particle
  double density = 1.66e-4;  // mass_asu/Å^3
  double radius = (2.79e4)/2;  // Å
  double volume = (4.0/3.0) * M_PI * pow(radius, 3.0);
  double mass = density * volume;

  /* Declaration of variables */
  double T = 297;  // K
  const double kb = 8.617333262145e-5;  // eV/K
  double tau = 147.3e6;  // picoseconds
  double eta = 1/tau;
  long int n_timesteps = 100000;
  double dt = 1e6;  // picoseconds
  double dtau = 50*dt;  // Sampling timestep (picoseconds) for fft
  double omega = 2*M_PI*3.1e-9;  // picoseconds^-1
  double c0 = exp(-eta*dt);
  double vth = sqrt(kb*T/mass);
  double equilibration_time = 1e10;  // 1e9 is 1 ms.
  // timestep after which equilibration is already achieved
  int equilibration_i = equilibration_time/dt;
  long int N_POINTS = (n_timesteps+1)-equilibration_i;  // For fft

  // RNG's sigmas
  double sigma = 1.0;

  /* ARRAYS: Initialize and allocate memory */
  double *q_arr = calloc(n_timesteps+1, sizeof(double));
  double *v_arr = calloc(n_timesteps+1, sizeof(double));
  double *a_arr = calloc(n_timesteps+1, sizeof(double));
  double time_array[n_timesteps];
  arange(time_array, 0, n_timesteps, dt);

  /* Initial conditions */
  q_arr[0] = 10.0;
  v_arr[0] = 10.0;

  printf("Report of Simulation Variables:\n\n");
  printf("Total Simulation Time: %f ps\n", n_timesteps*dt);
  printf("Time step size:        %f ps\n", dt);
  printf("Number of time steps:  %ld\n", n_timesteps);
  printf("Temperature:           %f K\n", T);
  printf("Tau:                   %f ps\n", tau);
  printf("Equilibration Time:    %f ps\n", equilibration_time);
  printf("Equilibration timestep:%d ps\n", equilibration_i);
  printf("Initial position:      %f Å\n", q_arr[0]);
  printf("Initial velocity:      %f Å/ps\n", v_arr[0]);



  //printf("omega times 10^9: %lf\n", omega*(1e9));
  //printf("pi: %f\n", M_PI);

  // Solve Langevin's equation using BD algorithm on page 16.
  BD3(n_timesteps, dt, q_arr, v_arr, a_arr, omega, c0, vth, eta,
      r1, r2, sigma);

  /* Slice arrays to throw away data points before equilibration.
   * equilibration_i is the place you want to start your subset.
   */

  double *sliced_q = calloc((n_timesteps+1)-equilibration_i, sizeof(double));
  double *sliced_v = calloc((n_timesteps+1)-equilibration_i, sizeof(double));
  double *sliced_a = calloc((n_timesteps+1)-equilibration_i, sizeof(double));
  double *sliced_time = calloc((n_timesteps+1)-equilibration_i, sizeof(double));
  printf("Memory allocated for sliced\n");

  for(long int j = equilibration_i; j < (n_timesteps+1); j++)
  {
    sliced_q[j-equilibration_i] = q_arr[j];
    sliced_v[j-equilibration_i] = v_arr[j];
    sliced_a[j-equilibration_i] = a_arr[j];
    sliced_time[j-equilibration_i] = time_array[j] - equilibration_time;
  }
  printf("Sliced arrays created\n");

  /* Write trajectories, velocities, and accelerations to file.
   */
  write_trajectories_file("./out/non_eq_trajectories.csv", time_array,
                q_arr, v_arr, a_arr, n_timesteps);
  printf("Non-equilibrated files created\n");

  /* After equilibration, write trajectories, velocities,
   * and accelerations to file.
   */
  write_trajectories_file("./out/trajectories.csv", sliced_time,
                sliced_q, sliced_v, sliced_a, (n_timesteps-equilibration_i));
  printf("Equilibrated files created\n");
}

/*
 * Perform the Brownian Dynamics 3 algorithm to solve Langevin's equation.
 * RETURNS: void, but fills in @q_arr and @v_arr
 * @n_timesteps - The number of time steps to be performed
 * @dt - timestep
 * @v_arr - velocity array with initial condition
 * @q_arr - position array with initial condition
 * @omega - natural radian frequency of the particle under the harmonic force
 * @c0 - exponent of the friction coefficient and dt
 * @vth - thermal velocity

 */
void BD3(long int n_timesteps, double dt,
         double *q_arr, double *v_arr, double *a_arr,
         double omega, double c0, double vth, double eta,
         gsl_rng * r1, gsl_rng * r2, double sigma)
{
    // Initialize variables
    double q = q_arr[0];
    double v = v_arr[0];
    double c0_sqrt = sqrt(c0);
    double a = calc_acc(q, omega);
    //printf("%f\n", a*(1e12));
    a_arr[0] = a;
    //printf("%f\n", a_arr[0]*(1e12));

    for (int i = 1; i < n_timesteps + 1; i++)
    {
        double gauss1 = gsl_ran_gaussian(r1, sigma);
        double gauss2 = gsl_ran_gaussian(r2, sigma);

        /* v(t+dt/2) */
        v = (dt * 0.5 * a) + (c0_sqrt * v) + (vth * sqrt(1-c0) * gauss1);

        /* q(t+dt) */
        q += dt * v;

        /* a(t+dt) */
        a = calc_acc(q, omega);
        //printf("a_outside: %f\n", a*(1e15));

        v = (0.5 * c0_sqrt * a * dt) + (c0_sqrt * v)
            + (vth * sqrt(1-c0) * gauss2);

        //printf("a: %f v: %f q: %f\n", a, v, q);
	      /* Save position and velocity to arrays */
        q_arr[i] = q;
        v_arr[i] = v;
        a_arr[i] = a;
    }

}

/*
 * Calculate the acceleration on a particle experiencing three different forces
 * harmonic force, friction force, and stochastic force
 * @a - acceleration
 * @q - current position
 * @omega - natural radian frequency of particle under harmonic force
 * @eta - Friction coefficient
 * @v - velocity
 */
double calc_acc(double q, double omega)
{
    // Initialize variables
    double omega_sq = pow(omega, 2.0);
    double a = -omega_sq * q;
    //printf("xi: %f, omega_sq: %f, eta: %f, a: %f\n", xi, omega_sq, eta, a);
    //printf("a_inside: %f\n", a*(1e15));
    return a;
}

/*
 * constructs time array
 * @array - array to be filled with time values
 * @start - start value
 * @len_t - number of times stamps in array
 * @dt - time step between two consecutive times
*/
void arange(double *array, double start, int len_t, double dt){
    for(int i = 0; i < len_t; i++){
        array[i] = start + i*dt;
    }
}

/*
 * Writes trajectories and velocities to file
 * Unit conversion is done to match the plots on the Science paper.
 * position vs time: (nm) vs (ms)
 * velocity vs time: (mm/s) vs (ms)
 *
 * @fname - File name
 * @time_array - array of time values
 * @position - array with particle's positions
 * @velocity - array with particle's velocities
 * @n_timesteps - number of timesteps
*/
void write_trajectories_file(char *fname, double *time_array,
			   double *position, double *velocity, double *acceleration,
		     long int n_timesteps)
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp,
          "time (ms), position (nm), velocity (mm/s), acceleration (Å/ps^2)\n");
    for(int i = 0; i < n_timesteps; ++i)
    {
      fprintf(fp, "%f,%f,%f,%f\n", time_array[i]*(1e-9),
	      position[i]*(1e-1), velocity[i]*(1e5), acceleration[i]);
    }
    fclose(fp);
}

/*
 * reads time_array and signal data from file
 * @fname - File name
 * @time_array - array of time values
 * @signal - array with signal values
*/
void read_data(char *fname, double *time_array, double *signal)
{
    FILE *fp = fopen(fname, "r");

    /* if file no found
     * error out and exit code 1
     */
    if(fp == NULL){
	perror("error:");
	exit(1);
    }

    /* skip header */
    fseek(fp, strlen("time (ms), velocity (mm/s)\n"), SEEK_SET);
    char line[128] = {0};
    char *token;
    int i = 0;
    while(fgets(line, sizeof(line), fp) != NULL){
	token = strtok(line, ",");
	time_array[i] = strtod(token, NULL);
	token = strtok(NULL, ",");
	signal[i] = strtod(token, NULL);
	i++;
	memset(line, 0, sizeof(line));
	token = NULL;
    }
    fclose(fp);
}

/*
 * constructs time array
 * @fname - File name
 * @frequencies - array of frequency values
 * @power - array with power values
 * @n_points - number of points
*/
void write_spectrum_file(char *fname, double *frequencies,
		   double *power, int n_points)
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "frequencies, power\n");
    for(int i = 0; i < n_points; ++i){
	    fprintf(fp, "%f,%f\n", frequencies[i], power[i]);
    }
    fclose(fp);
}
