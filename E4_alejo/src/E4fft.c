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

/*************************************************************
 * Macro defines
 *************************************************************/
/* DT is 50 times the dt=1e6 in E4main.c */
#define DT 5e7
#define N_POINTS 50000

/******************************************************************************
 * Helper functions
 *****************************************************************************/
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
    fseek(fp, strlen("time, signal\n"), SEEK_SET);
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
 * @time_array - array of time values
 * @signal - array with signal values
 * @n_points - number of points
*/
void write_to_file(char *fname, double *frequencies,
		   double *spectrum, int n_points)
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "time, signal\n");
    for(int i = 0; i < n_points; ++i){
	    fprintf(fp, "%f,%f\n", frequencies[i], spectrum[i]);
    }
    fclose(fp);
}

/**************************************************************
 * Main routine
 **************************************************************/
int main(int argc, char **argv)
{
  /*
   * Code from where the signal to analyze comes from
   */
    bool code1 = false;
    bool code4 = true; // Check code4's output file's name
    bool kinetic_spectrum = false;

    double time_array[N_POINTS];
    double signal[N_POINTS];
    if(code1)
    {
      read_data("signal.csv", time_array, signal);
    }

    if(code4)
    {
      read_data("abs_total_displacements.csv", time_array, signal);
    }

    if(kinetic_spectrum)
    {
      read_data("kinetic.csv", time_array, signal);
    }

    /*
     * Construct array with frequencies
     */
    double frequencies[N_POINTS];
    fft_freq_shift(frequencies, DT, N_POINTS);

    /*
     * Do the fft
     */
    double fftd_data[N_POINTS];
    powerspectrum(signal, fftd_data, N_POINTS);
    powerspectrum_shift(fftd_data, N_POINTS);
    /*
     * Dump fft and frequencies to file depending on signal analyzed
     */
     if(code1){
       write_to_file("powerspectrum_code1.csv",
       frequencies, fftd_data, N_POINTS);
     }

     if(code4)
     {
       write_to_file("powerspectrum_code4.csv",
       frequencies, fftd_data, N_POINTS);
     }

     if(kinetic_spectrum)
     {
       write_to_file("powerspectrum_kinetic.csv",
       frequencies, fftd_data, N_POINTS);
     }
    return 0;
}