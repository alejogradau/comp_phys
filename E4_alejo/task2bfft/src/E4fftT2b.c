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

   // Declare Variables from its corresponding E4main.c
   double dt = 0.001;  // miliseconds
   int factor = 50;
   char cases[80];
   sprintf(cases, "high");
   double dtau = factor*dt;  // Sampling timestep (miliseconds) for fft
   double n_timesteps = 1e5;
   int M = 100;  // Number of segments of long trajectory
   //long int total_timesteps = n_timesteps*M;
   double equilibration_time = 10;  // ms.
   int equilibration_i = equilibration_time/dt;
   long int N_POINTS = n_timesteps- (equilibration_i/M);  // For fft
   long int SAMP_POINTS = N_POINTS/factor;
   printf("%ld\n", N_POINTS);

   /* ARRAYS: Initialize and allocate memory */
   double *time_array = calloc(N_POINTS, sizeof(double));
   double *signal = calloc(N_POINTS, sizeof(double));
   double *signal_sampled = calloc(SAMP_POINTS, sizeof(double));

   /*
   double time_array[N_POINTS];
   double time_array_sampled[SAMP_POINTS];
   double signal[N_POINTS];
   double signal_sampled[SAMP_POINTS];
   */

   char fname1[80];
   char fname2[80];

   // Perform a fft on each of the M velocity.csv files
   for(int k = 0; k < M; ++k)
   {
     sprintf(fname1, "./out/velocities_%s_%d.csv", cases, k);
     printf("Reading %s file for FFT\n", fname1);
     read_data(fname1, time_array, signal);
     printf("%s Read\n", fname1);

     printf("Sampling %s using dtau = %d*dt\n", fname1, factor);
     for(int i = 0; i < SAMP_POINTS; i++)
     {
       signal_sampled[i] = signal[i*factor];
     }

    /*
     * Construct array with frequencies
     */
     double frequencies[SAMP_POINTS];
     printf("Shifting frequencies\n");
     fft_freq_shift(frequencies, dtau, SAMP_POINTS);

     /*
      * Do the fft
      */
      double fftd_data[SAMP_POINTS];
      printf("Performing FFT\n");
      powerspectrum(signal_sampled, fftd_data, SAMP_POINTS);
      powerspectrum_shift(fftd_data, SAMP_POINTS);

     /*
      * Write fft and frequencies to file depending on signal analyzed
      */
      sprintf(fname2, "./out/powerspectrum_%s_%d_%ddt.csv", cases, k, factor);
      printf("Writing %s file\n", fname2);
      write_to_file(fname2, frequencies, fftd_data, SAMP_POINTS);
      printf("ok\n\n");
   }

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

/*
void file_to_array(char *fname, long int N_POINTS, double data[arr_len])
{
    FILE *fp;
    fp = fopen(fname, "r");

    for (int i = 0; i < length; i++) {
        fscanf(fp, "%lf", &data[i]);
    }
    printf("Reading file with: %d lines\n", length);
}
*/
