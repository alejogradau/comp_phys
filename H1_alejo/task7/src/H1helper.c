/******************************************************************************
 * Helper functions
 *****************************************************************************/

#include <stdio.h>
#include <math.h>  // round() function
#include <stdlib.h>

 /*
  * constructs a0 array
  * @array - array to be filled with values
  * @start - start value
  * @n_points - number of stamps in array
  * @step - step between two consecutive values
 */
 void arange(double *array, double start, int n_points, double step){
     for(int i = 0; i < n_points; i++){
     array[i] = start + i*step;
     }
 }

/*
 * Writes energies to file
 * @fname - File name
 * @time_array - array of time values
 * @n_timesteps - number of timesteps
 * @T[n_timesteps] - array containing the kinetic energy at each timestep
 * @V[n_timesteps] - array containing the potential energy at each timestep
 * @E[n_timesteps] - array containing the total energy at each timestep
*/
void write_energies_file(char *fname, double *time_array, int n_timesteps, double T[n_timesteps], double V[n_timesteps], double E[n_timesteps])
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "time, T, V, E\n");

    for(int i = 0; i < n_timesteps; ++i){
        fprintf(fp, "%f, %f, %f, %f\n", time_array[i], T[i], V[i], E[i]);
    }
    fclose(fp);
}

/*
 * Writes observable to file
 * @fname - File name
 * @fname - Observable name
 * @time_array - array of time values
 * @n_timesteps - number of timesteps
 * @O[n_timesteps] - array containing the value of the observable at each timestep
*/
void write_observable_file(char *fname, char* oname, double *time_array, int n_timesteps, double O[n_timesteps])
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "time, %s, \n", oname);

    for(int i = 0; i < n_timesteps; ++i){
        fprintf(fp, "%f, %f\n", time_array[i], O[i]);
    }
    fclose(fp);
}

/*
 * Calculates the time average of an observable
 * @n_timesteps - number of timesteps
 * @O[n_timesteps] - array containing the value of the observable at each timestep
 * @O_exp[n_timesteps] - array where the time average will be saved
*/
void calc_time_average(int n_timesteps, double O[n_timesteps], double O_exp[n_timesteps]){
    double sum_of_O = 0;

    //At each timestep
    //Adds O at current time step to all the previous values
    //Divides by the number of the current timestep
    for(int i = 0; i < n_timesteps; i++){
            sum_of_O += O[i];
            O_exp[i] = sum_of_O/(i+1);
    }
}

/*
 * Calculate the distance between all pairs of particles and return an array
 * with these distances.
 * @N - number of particles
 * @positions[N][3] - matrix of atomic positions
 * @distances[N*(N-1)/2] - array that stores the distance between all pairs of
 * particles
 *
*/
void distances_array(unsigned int N, double *distances[N*(N-1)/2],
  double positions[N][3])
{
  int i,j,k;
  int distances_len = N*(N-1)/2;
  for (k = 0; k < distances_len; k++)
  {
    for (i = 0; i < N; i++)
    {
      for (j = i + 1; j < N; j++)
      {
        double rxij = positions[i][0] - positions[j][0];
        double ryij = positions[i][1] - positions[j][1];
        double rzij = positions[i][2] - positions[j][2];
        double rij_sq = (rxij*rxij + ryij*ryij + rzij*rzij);
        double rij = sqrt(rij_sq);
      }
    }
    distances[k] = rij;
  }
}

/*
 * Calculate the distance between all pairs of particles and return an array
 * with the integer distances k, given by formula modified from H1a.
 * @N - number of particles
 * @positions[N][3] - matrix of atomic positions
 * @dr - bin size
 * @int_distances[N*(N-1)/2] - array that stores the integer distance between
 * all pairs of particles
 *
 */
void int_distances_array(unsigned int N, double *int_distances[N*(N-1)/2],
  double positions[N][3], double dr)
{
  int i,j,l;
  int array_len = N*(N-1)/2;
  for (l = 0; l < array_len; l++)
  {
    for (i = 0; i < N; i++)
    {
      for (j = i + 1; j < N; j++)
      {
        double rxij = positions[i][0] - positions[j][0];
        double ryij = positions[i][1] - positions[j][1];
        double rzij = positions[i][2] - positions[j][2];
        double rij_sq = (rxij*rxij + ryij*ryij + rzij*rzij);
        double rij = sqrt(rij_sq);
      }
    }
    int_distances[l] = round((rij/dr) - 0.5);
  }
}
