/******************************************************************************
 * E1code3
 ******************************************************************************
 * reads h(t) = a*cos(2*pi*f*t + phi)
 * runs the fft of h(t)
 *
 *
 * Compile me as:
 * clang -c fft.c -o fft.o -lgsl -lgslcblas
 * clang E1code3.c fft.o -o code3 -lgsl -lgslcblas
 * Alejo: primera linea compila librerias. Si no has
 * cambiado nada de las librerias no tienes que hacer esto
 * otra vez
*/

/************************************************************
 * Includes
 ************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h> //bool

#include "fft.h"  // interface to fft routine

/******************************************************************************
 * Helper functions
 *****************************************************************************/
 void read_data(char *fname, double *time_array, double *signal);
 void write_to_file(char *fname, double *frequencies,
 		   double *power, int n_points);

/**************************************************************
 * Main routine
 **************************************************************/
int main(int argc, char **argv)
{
  /*
   * Code from where the signal to analyze comes from
   */

   // Declare Variables
   double dt = 0.001;  // miliseconds
   double dtau = 5.0*dt;  // Sampling timestep (miliseconds) for fft
   long int n_timesteps = 100000;
   double equilibration_time = 10;  // ms.
   int equilibration_i = equilibration_time/dt;
   long int N_POINTS = (n_timesteps+1)-equilibration_i;  // For fft

   double time_array[N_POINTS];
   double signal[N_POINTS];

   printf("Reading velocity signal file for FFT\n");
   read_data("./out/velocities.csv", time_array, signal);

  /*
   * Construct array with frequencies
   */
   double frequencies[N_POINTS];
   printf("Shifting frequencies\n");
   fft_freq_shift(frequencies, dtau, N_POINTS);

  /*
   * Do the fft
   */
   double fftd_data[N_POINTS];
   printf("Performing FFT\n");
   powerspectrum(signal, fftd_data, N_POINTS);
   powerspectrum_shift(fftd_data, N_POINTS);

  /*
   * Write fft and frequencies to file depending on signal analyzed
   */
   printf("Writing file\n");
   write_to_file("./out/powerspectrum.csv", frequencies, fftd_data, N_POINTS);

   return 0;
}

/*
 * reads time_array and signal data from file
 * @fname - File name
 * @time_array - array of time values
 * @signal - array with velocity signal values
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
void write_to_file(char *fname, double *frequencies,
		   double *power, int n_points)
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "frequencies, power\n");
    for(int i = 0; i < n_points; ++i){
	    fprintf(fp, "%f,%f\n", frequencies[i], power[i]);
    }
    fclose(fp);
}
