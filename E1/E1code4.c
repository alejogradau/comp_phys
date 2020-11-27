/******************************************************************************
 * E1code4
 ******************************************************************************
 * Routine that runs the velocity verlet algorithm
 * Use as template to construct your program!
 *
 * Compile me as:
 * clang E1code4.c -o code4 -lm
 */

/******************************************************************************
 * Includes
 *****************************************************************************/
#include <stdio.h>   //fopen, fprintf
#include <stdlib.h>  //malloc
#include <stdint.h>  //uint64_t
#include <math.h>    //pow, fabs function (absolute values double)

/*
 * Calculate the acceleration
 * @a - vector that is filled with acceleration
 * @u - vector with the current positions
 * @m - vector with masses
 * @kappa - Spring constant
 * @size_of_u - the size of the position, acceleration and mass array
 */
void calc_acc(double *a, double *u, double *m, double kappa, int size_of_u)
{
    /* Declaration of variables */
    int i;

    /* Calculating the acceleration on the boundaries */
    a[0] = kappa*(- 2*u[0] + u[1])/m[0];
    a[size_of_u - 1] = kappa*(u[size_of_u - 2] - 2*u[size_of_u - 1])/m[size_of_u - 1];

    /* Calculating the acceleration of the inner points */
    for (i = 1; i < size_of_u - 1; i++){
        a[i] = kappa*(u[i - 1] - 2*u[i] + u[i + 1])/m[i];
    }
}

/*
 * Perform the velocity verlet alogrithm
 * @n_timesteps - The number of time steps to be performed
 * @n_particles - number of particles in the system
 * @v - array of velocity (Empty allocated array) : sizeof(v) = n_particles
 * @q_n - position of the n'th atom : sizeof(q_n) = n_timesteps+1
 * @dt - timestep
 * @m - vector with masses of atoms sizeof(n_particles)
 * @kappa - Spring constant
 */
void velocity_verlet(int n_timesteps, int n_particles, double *v,
		     double *q_1, double *q_2, double *q_3, double *abs_q_total,
		     double dt, double *m, double kappa)
{
    double q[n_particles];
    double a[n_particles];
    q[0] = q_1[0];
    q[1] = q_2[0];
    q[2] = q_3[0];
    calc_acc(a, q, m, kappa, n_particles);
    for (int i = 1; i < n_timesteps + 1; i++) {
        /* v(t+dt/2) */
        for (int j = 0; j < n_particles; j++) {
            v[j] += dt * 0.5 * a[j];
        }

        /* q(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            q[j] += dt * v[j];
        }

        /* a(t+dt) */
        calc_acc(a, q, m, kappa, n_particles);

        /* v(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            v[j] += dt * 0.5 * a[j];
        }

        /* Save the displacement of the three atoms */
        q_1[i] = q[0];
        q_2[i] = q[1];
        q_3[i] = q[2];

	/* Save absolute total displacement */
        for (int j = 0; j < n_particles; j++) {
            abs_q_total[i] += fabs(q[j]);
        }
    }
}

/*
 * Perform the verlet alogrithm to calculate energies
 * @n_timesteps - The number of time steps to be performed
 * @n_particles - number of particles in the system
 * @v - array of velocity (Empty allocated array) : sizeof(v) = n_particles
 * @T_i - total kinetic E : sizeof(T) = n_timesteps+1
 * @q_i - position of the ith atom : sizeof(q_i) = n_timesteps+1
 * @dt - timestep
 * @m - vector with masses of atoms sizeof(n_particles)
 * @kappa - Spring constant
 */
void energy_verlet(int n_timesteps, int n_particles, double *v,
		   double *T, double *V, double *E,
		   double *q_1, double *q_2, double *q_3,
		   double dt, double *m, double kappa)
{
    double q[n_particles];
    double a[n_particles];
    q[0] = q_1[0];
    q[1] = q_2[0];
    q[2] = q_3[0];
    calc_acc(a, q, m, kappa, n_particles);
    for (int i = 1; i < n_timesteps + 1; i++) {
        /* v(t+dt/2) */
        for (int j = 0; j < n_particles; j++) {
            v[j] += dt * 0.5 * a[j];
        }

        /* q(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            q[j] += dt * v[j];
        }

        /* a(t+dt) */
        calc_acc(a, q, m, kappa, n_particles);

        /* v(t+dt) */
        for (int j = 0; j < n_particles; j++) {
	  v[j] += dt * 0.5 * a[j];
        }

	/* calculate and save the total kinetic energy at time i; T */
        for (int j = 0; j < n_particles; j++) {
	  T[i] += m[j] * v[j] * v[j] / 2.0;
        }

	/* calculate and save the total potential energy at time i; V */
        V[i] += kappa * (pow(q[0],2) + pow(q[1]-q[0],2)
			 + pow(q[2]-q[1],2) + pow(q[2],2)) / 2.0;

	/* calculate and save the total energy at time i; E */
	  E[i] += T[i] + V[i];

    }
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
 * Writes positions file
 * @fname - File name
 * @time_array - array of time values
 * @position_i - array with i'th atom positions
 * @n_timesteps - number of timesteps
*/
void write_positions_file(char *fname, double *time_array,
			   double *position_1, double *position_2, double *position_3,
		     int n_timesteps)
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "time, q_1, q_2, q_3\n");
    for(int i = 0; i < n_timesteps; ++i){
      fprintf(fp, "%f,%f,%f,%f\n", time_array[i],
	      position_1[i], position_2[i], position_3[i]);
    }
    fclose(fp);
}

/*
 * Writes total_position.csv file exactly as code1 writes signal file.
 * This way code3 just needs to look for another file name.
 * @fname - File name
 * @time_array - array of time values
 * @position_i - array with i'th atom positions
 * @n_timesteps - number of timesteps
*/
void write_total_displacements_file(char *fname, double *time_array,
			  double *abs_q_total, int n_timesteps)
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "time, signal\n");
    for(int i = 0; i < n_timesteps; ++i){
      fprintf(fp, "%f,%f\n", time_array[i], abs_q_total[i]);
    }
    fclose(fp);
}

/*
 * Writes energies file
 * @fname - File name
 * @time_array - array of time values
 * @position_i - array with i'th atom positions
 * @n_timesteps - number of timesteps
*/
void write_energy_file(char *fname, double *time_array,
		       double *T, double *V, double *E,
		       int n_timesteps)
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "time, T, V, E\n");
    for(int i = 0; i < n_timesteps; ++i){
      fprintf(fp, "%f,%f,%f,%f\n", time_array[i], T[i], V[i], E[i]);
    }
    fclose(fp);
}

/*
 * Writes q1 energy oscillation file
 * @fname - File name
 * @time_array - array of time values
 * @position_i - array with i'th atom positions
 * @n_timesteps - number of timesteps
*/
void write_q_file(char *fname, double *time_array,
		       double *q, int n_timesteps)
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "time, signal\n");
    for(int i = 0; i < n_timesteps; ++i){
      fprintf(fp, "%f,%f\n", time_array[i], q[i]);
    }
    fclose(fp);
}

int main()
{
  /* Declaration of variables */
  int n_timesteps = 50000; double dt = 0.00005;
  int n_particles = 3; double kappa = 6.25e1;
  double m[n_particles];
  double v[n_particles];
  double q_1[n_timesteps+1];
  double q_2[n_timesteps+1];
  double q_3[n_timesteps+1];
  /* Allocate memory for energy arrays */
  double *T = calloc(n_timesteps+1, sizeof(double));
  double *V = calloc(n_timesteps+1, sizeof(double));
  double *E = calloc(n_timesteps+1, sizeof(double));
  double *abs_q_total = calloc(n_timesteps+1, sizeof(double));
  double time_array[n_timesteps];
  arange(time_array, 0, n_timesteps, dt);

  /* Instantiate variables */
  /* Mass of carbon atom in asu */
  m[0] = 1.244e-3;
  m[1] = 1.244e-3;
  m[2] = 1.244e-3;

  /* Position initial conditions */
  q_1[0] = 0.01;
  q_2[0] = 0;
  q_3[0] = 0;

  /* Velocity initial conditions */
  v[0] = 0;
  v[1] = 0;
  v[2] = 0;

  /* Kinetic energy based on initial conditions; T[0] */
  for (int j = 0; j < n_particles; j++) {
    T[0] += m[j] * v[j] * v[j] / 2.0;
  }

  /* Potential energy based on initial conditions; V[0] */
  V[0] += kappa * (pow(q_1[0],2) + pow(q_2[0]-q_1[0],2)
   + pow(q_3[0]-q_2[0],2) + pow(q_3[0],2)) / 2.0;

   /* Total energy based on initial conditions; E[0] */
   E[0] += T[0] + V[0];

  velocity_verlet(n_timesteps, n_particles, v,
		  q_1, q_2, q_3, abs_q_total,
		  dt, m, kappa);

  /* Velocity initial conditions (again) */
  v[0] = 0;
  v[1] = 0;
  v[2] = 0;

  energy_verlet(n_timesteps, n_particles, v,
		T, V, E,
		q_1, q_2, q_3,
		dt, m, kappa);

  write_positions_file("positions.csv", time_array, q_1, q_2, q_3,
		n_timesteps);

  write_total_displacements_file("abs_total_displacements.csv", time_array,
    			  abs_q_total, n_timesteps);

  write_energy_file("energy.csv", time_array, T, V, E, n_timesteps);

  write_q_file("q1.csv", time_array, q_1, n_timesteps);
  write_q_file("q2.csv", time_array, q_2, n_timesteps);
  write_q_file("q3.csv", time_array, q_3, n_timesteps);
}
