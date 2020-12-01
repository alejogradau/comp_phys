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
