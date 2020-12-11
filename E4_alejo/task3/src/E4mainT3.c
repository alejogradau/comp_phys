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
#include <math.h>

#include "fft.h"  // interface to fft routine

/******************************************************************************
 * Helper functions
 *****************************************************************************/
 void read_data(char *fname, double *time_array, double *signal);
 void write_to_file(char *fname, double *frequencies,
 		   double *power, int n_points);
 void calc_autocorrelation(unsigned long length, double time_series[length],
   double autocorrelation[length/2]);

/**************************************************************
 * Main routine
 **************************************************************/
int main(int argc, char **argv)
{
  /*
   * Code from where the signal to analyze comes from
   */

   // Declare Variables
   char cases[80];
   double dt = 0.001;  // miliseconds
   int factor = 50;
   double dtau = factor*dt;  // Sampling timestep (miliseconds) for fft
   long int n_timesteps = 300000;
   double equilibration_time = 10;  // ms.
   int equilibration_i = equilibration_time/dt;
   long int length = n_timesteps-equilibration_i;  // Modified to match E3 autocorrelation (N_POINTS before)
   long int SAMP_POINTS = length/factor;
   printf("%ld\n", length);

   /* ARRAYS: Initialize and allocate memory */
   double *time_array = calloc(length, sizeof(double));
   double *signal = calloc(length, sizeof(double));
   double *autocorrelation = calloc(length/2, sizeof(double));
   double *signal_sampled = calloc(SAMP_POINTS, sizeof(double));

   char fname1[80];
   char fname2[80];
   char fname3[80];

   /*
    * CASE LOW
    */
   //sprintf(cases, "low");
   //sprintf(fname1, "./out/velocities_%s.csv", cases);
   //printf("Reading %s file for autocorrelation\n", fname1);
   //read_data(fname1, time_array, signal);
   //printf("%s Read\n", fname1);

   //calc_autocorrelation(length, signal, autocorrelation);

   //printf("Sampling autocorrelation_%s.csv using dtau = %d*dt\n", cases, factor);
   //for(int i = 0; i < SAMP_POINTS; i++)
   //{
  //   signal_sampled[i] = signal[i*factor];
   //}

   /*
    * Construct array with frequencies
    */
    double frequencies[SAMP_POINTS];
    //printf("Shifting frequencies\n");
    //fft_freq_shift(frequencies, dtau, SAMP_POINTS);

    /*
     * Do the fft
     */
     double fftd_data[SAMP_POINTS];
     //printf("Performing FFT\n");
     //powerspectrum(signal_sampled, fftd_data, SAMP_POINTS);
     //powerspectrum_shift(fftd_data, SAMP_POINTS);

  /*
   * Write autocorrelation to file
   */
   //sprintf(fname2, "./out/autocorrelation_%s.csv", cases);
   //printf("Writing %s file\n", fname2);
   //write_to_file(fname2, time_array, autocorrelation, length/2);

   /*
    * Write fft and frequencies to file depending on signal analyzed
    */
    //sprintf(fname3, "./out/powerspectrum_%s_%ddt.csv", cases, factor);
    //printf("Writing %s file\n", fname3);
    //write_to_file(fname3, frequencies, fftd_data, SAMP_POINTS);


   /*
    * CASE HIGH
    */

   sprintf(cases, "low");
   sprintf(fname1, "./out/velocities_%s.csv", cases);
   printf("Reading %s file for autocorrelation\n", fname1);
   read_data(fname1, time_array, signal);
   printf("%s Read\n", fname1);

   calc_autocorrelation(length, signal, autocorrelation);

   printf("Sampling autocorrelation_%s.csv using dtau = %d*dt\n", cases, factor);
   for(int i = 0; i < SAMP_POINTS; i++)
   {
     signal_sampled[i] = signal[i*factor];
   }

   /*
    * Reconstruct array with frequencies
    */
    printf("Shifting frequencies\n");
    fft_freq_shift(frequencies, dtau, SAMP_POINTS);

    /*
     * Do the fft
     */
     printf("Performing FFT\n");
     powerspectrum(signal_sampled, fftd_data, SAMP_POINTS);
     powerspectrum_shift(fftd_data, SAMP_POINTS);

  /*
   * Write autocorrelation to file
   */
   sprintf(fname2, "./out/autocorrelation_%s.csv", cases);
   printf("Writing %s file\n", fname2);
   write_to_file(fname2, time_array, autocorrelation, length/2);

   /*
    * Write fft and frequencies to file depending on signal analyzed
    */
    sprintf(fname3, "./out/powerspectrum_%s_%ddt.csv", cases, factor);
    printf("Writing %s file\n", fname3);
    write_to_file(fname3, frequencies, fftd_data, SAMP_POINTS);


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

void calc_autocorrelation(unsigned long length, double time_series[length],
  double autocorrelation[length/2])
  {
    long double average_f = 0, average_f2 = 0, average_fifik = 0, var_f, autocsum = 0, k_max = 30000;
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
