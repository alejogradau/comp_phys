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


/******************************************************************************
 * Helper functions
 *****************************************************************************/
double calc_acc(double q, double omega);
void arange(double *array, double start, int len_t, double dt);
void write_to_file(char *fname, double *time_array,
			   double *position, double *velocity, double *acceleration,
		     double n_timesteps);
void write_fft_M_files(double *time_array, double *velocity,
         		           double n_timesteps, int M, char *cases);
void BD3(double n_timesteps, double dt,
         double *q_arr, double *v_arr,
         double omega, double c0, double vth, double eta,
         gsl_rng * r1, gsl_rng * r2, double sigma);

int main()
{
  printf("main\n");
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
  char cases[80];
  sprintf(cases, "high");
  double tau = 48.5e6;  // picoseconds
  double eta = 1/tau;
  double n_timesteps = 1e5;  // Updated below
  int M = 100;  // Number of segments of long trajectory
  n_timesteps = M*n_timesteps;
  double dt = 1e6;  // picoseconds
  double omega = 2*M_PI*3.1e-9;  // picoseconds^-1
  double c0 = exp(-eta*dt);
  double vth = sqrt(kb*T/mass);
  double equilibration_time = 1e10;  // 1e9 is 1 ms.
  // timestep after which equilibration is already achieved
  int equilibration_i = equilibration_time/dt;

  // RNG's sigmas
  double sigma = 1.0;

  /* ARRAYS: Initialize and allocate memory */
  double *q_arr = calloc(n_timesteps+1, sizeof(double));
  double *v_arr = calloc(n_timesteps+1, sizeof(double));
  //double *a_arr = calloc(n_timesteps+1, sizeof(double));
  double *time_array = calloc(n_timesteps, sizeof(double));
  arange(time_array, 0, n_timesteps, dt);

  /* Initial conditions */
  q_arr[0] = 10.0;
  v_arr[0] = 10.0;

  printf("Report of Simulation Variables:\n\n");
  printf("Total Simulation Time: %f ps\n", n_timesteps*dt);
  printf("Time step size:        %f ps\n", dt);
  printf("Number of time steps:  %f\n", n_timesteps);
  printf("Number of segments:    %d\n", M);
  printf("Temperature:           %f K\n", T);
  printf("Tau:                   %f ps\n", tau);
  printf("Equilibration Time:    %f ps\n", equilibration_time);
  printf("Equilibration timestep:%d\n", equilibration_i);
  printf("Initial position:      %f Å\n", q_arr[0]);
  printf("Initial velocity:      %f Å/ps\n\n", v_arr[0]);



  //printf("omega times 10^9: %lf\n", omega*(1e9));
  //printf("pi: %f\n", M_PI);

  // Solve Langevin's equation using BD algorithm on page 16.
  BD3(n_timesteps, dt, q_arr, v_arr, omega, c0, vth, eta,
      r1, r2, sigma);

  /* Slice arrays to throw away data points before equilibration.
   * equilibration_i is the place you want to start your subset.
   */
  // WORKS WITH n_timesteps+1, but I believe you don't need that +1, you arent using it
  double *sliced_q = calloc((n_timesteps+1)-equilibration_i, sizeof(double));
  double *sliced_v = calloc((n_timesteps+1)-equilibration_i, sizeof(double));
  double *sliced_time = calloc((n_timesteps+1)-equilibration_i, sizeof(double));
  printf("Memory allocated for sliced\n");

  for(long int j = equilibration_i; j < (n_timesteps+1); j++)
  {
    sliced_q[j-equilibration_i] = q_arr[j];
    sliced_v[j-equilibration_i] = v_arr[j];
    sliced_time[j-equilibration_i] = time_array[j] - equilibration_time;
  }
  printf("Sliced arrays created\n");

  /* Write trajectories, velocities, and accelerations to file.
   */
  //write_to_file("./out/non_eq_trajectories.csv", time_array,
  //              q_arr, v_arr, a_arr, n_timesteps);
  //printf("Non-equilibrated files created\n");

  /* After equilibration, write trajectories, velocities,
   * and accelerations to file.
   */
  //write_to_file("./out/trajectories.csv", sliced_time,
  //              sliced_q, sliced_v, sliced_a, (n_timesteps-equilibration_i));
  //printf("Equilibrated files created\n");

  write_fft_M_files(sliced_time, sliced_v, (n_timesteps-equilibration_i), M, cases);
  printf("%d velocity signal files for FFT created in ./out/\n", M);

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
void BD3(double n_timesteps, double dt,
         double *q_arr, double *v_arr,
         double omega, double c0, double vth, double eta,
         gsl_rng * r1, gsl_rng * r2, double sigma)
{
    // Initialize variables
    double q = q_arr[0];
    double v = v_arr[0];
    double c0_sqrt = sqrt(c0);
    double a = calc_acc(q, omega);
    //printf("%f\n", a*(1e12));
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
void write_to_file(char *fname, double *time_array,
			   double *position, double *velocity, double *acceleration,
		     double n_timesteps)
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
 * Writes velocities to M segmented files for FFT analysis
 * Unit conversion is done to match the plots on the Science paper.
 * velocity vs time: (mm/s) vs (ms)
 *
 * @fname - File name
 * @time_array - array of time values
 * @velocity - array with particle's velocities
 * @n_timesteps - number of timesteps
*/
void write_fft_M_files(double *time_array, double *velocity,
		                double n_timesteps, int M, char *cases)
{
    char fname[80];

    for(int k = 0; k < M; ++k)
    {
      sprintf(fname, "../task2bfft/out/velocities_%s_%d.csv", cases, k);
      FILE *fp = fopen(fname, "w");
      fprintf(fp, "time (ms), velocity (mm/s)\n");

      for(int i = k*(n_timesteps/M); i < (k+1)*(n_timesteps/M); ++i)
      {
        fprintf(fp, "%f,%f\n", time_array[i]*(1e-9), velocity[i]*(1e5));
      }
      fclose(fp);
      printf("velocities_%s_%d.csv file created in ./out/\n", cases, k);
    }
}
