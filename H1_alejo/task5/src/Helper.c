/******************************************************************************
 * Helper functions
 *****************************************************************************/

#include <stdio.h>
#include <math.h> // pow

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

void write_temperatures_file(char *fname, double *time_array, int n_timesteps, double Temp[n_timesteps])
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "time, Temp\n");

    for(int i = 0; i < n_timesteps; ++i){
        fprintf(fp, "%f, %f\n", time_array[i], Temp[i]);
    }
    fclose(fp);
}

void calc_time_average(int n_timesteps, double O[], double O_exp[]){
    double sum_of_O = 0;

    for(int i = 0; i < n_timesteps; i++){
            sum_of_O += O[i];
            O_exp[i] = sum_of_O/(i+1);
    }
}

/* Calculate the average value of the quantity A after equilibration */
double calc_eq_average(int n_timesteps, double A[], unsigned int i_0){
    double sum_of_A = 0;
    double A_equilibration = 0;
    //unsigned int tau = 100;
    //unsigned int i_0 = 2*tau;  // Equilibration
    for(int i = i_0; i < n_timesteps; i++){
            sum_of_A += A[i];
    }
    A_equilibration = sum_of_A/(n_timesteps - i_0);
    return A_equilibration;
}

/* Calculate the variance of the average of the quantity A after
 * equilibration
 */
double calc_eq_var(int n_timesteps, double A[], double A_equilibration,
  unsigned int i_0){
    double sum = 0;
    double variance = 0;
    for(int i = i_0; i < n_timesteps; i++){
            sum += pow((A[i] - A_equilibration),2);
    }
    variance = sum/(n_timesteps - i_0);
    return variance;
}

/* Calculate heat capacity for the microcanonical ensemble (constant NVE)
 * formula is taken from eq. 58 in MD_LectureNotes
 */
double calc_cv_NVE(unsigned int n_particles, double temp_eq, double variance)
{
  const double kb = 8.617333262145e-5;
  double factor1 = (3.0*n_particles*kb)/2.0;
  double factor2 = 2.0 / ( 3.0*n_particles*pow(kb,2)*pow(temp_eq,2) );
  double cv = factor1 / (1 - (factor2*variance));
  return cv;
}
