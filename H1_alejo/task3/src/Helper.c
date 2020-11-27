/******************************************************************************
 * Helper functions
 *****************************************************************************/

#include <stdio.h>

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
 * Writes to file
 * @fname - File name
 * @time_array - array of time values
 * @n_timesteps - number of timesteps
 * @n_particles - number of particles samples
 * @samples - array of size time_steps x n_particles
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
