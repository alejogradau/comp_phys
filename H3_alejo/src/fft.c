/*
	fft.c
	Program with powerspectrum from (fast) discrete Fourier transform  using GSL
	Created by Martin Gren on 2014-10-22
*/

/******************************************************************************
 * Macro defines
 *****************************************************************************/
#define PI 3.14159
#define hbar 0.6582  // eV * fs

/******************************************************************************
 * Includes
 *****************************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <complex.h>

/* Declare variables */
double planck = 2 * PI * hbar;

/* Makes fft of data and returns the powerspectrum in powspec_data */
void powerspectrum(double *data, double complex *powspec_data, int n) /* input data, output powspec_data, number of timesteps */
{
	/* Declaration of variables */
	int i;
	double complex_coefficient[2*n]; // array for the complex fft data
	double data_cp[n];

	/*make copy of data to avoid messing with data in the transform*/
	for (i = 0; i < n; i++)	{
		data_cp[i] = data[i];
	}

	/*Declare wavetable and workspace for fft*/
	gsl_fft_real_wavetable *real;
	gsl_fft_real_workspace *work;

	/*Allocate space for wavetable and workspace for fft*/
	work = gsl_fft_real_workspace_alloc(n);
	real = gsl_fft_real_wavetable_alloc(n);

	/*Do the fft*/
	gsl_fft_real_transform(data_cp, 1, n, real, work);

	/*Unpack the output into array with alternating real and imaginary part*/
	gsl_fft_halfcomplex_unpack(data_cp, complex_coefficient,1,n);

	/*fill the output powspec_data with the powerspectrum */
	for (i = 0; i < n; i++)	{
		powspec_data[i] = (complex_coefficient[2*i]*complex_coefficient[2*i]+complex_coefficient[2*i+1]*complex_coefficient[2*i+1])/n;
	}

	/*Free memory of wavetable and workspace*/
	gsl_fft_real_wavetable_free(real);
	gsl_fft_real_workspace_free(work);
}


/* Shifts the input powspec_data to center the 0 momentum */
void powerspectrum_shift(double complex *powspec_data, int n) /* input data, timestep, number of timesteps */
{
	/* Declaration of variables */
	int i;

	/* make copy of fft_data as reference for the shift */
	double powspec_cp[n];
	for (i = 0; i < n; i++)	{
		powspec_cp[i] = powspec_data[i];
	}

	/* make shift */
	for (i = 0; i < n; i++)	{
		if (n % 2) /*if n odd*/	{
			if (i<=(n-2)/2)	{
				powspec_data[i] = powspec_cp[(i+(n+1)/2)];
			}
			else {
				powspec_data[i] = powspec_cp[(i+(n+1)/2)%(n)];
			}
		}
		else {
			if (i<n/2) {
				powspec_data[i] = powspec_cp[i+n/2];
			}
			else {
				powspec_data[i] = powspec_cp[(i+n/2)%(n)];
			}
		}
	}
}

/* Makes a momentum array fft_momentum with momentum interval planck/(dx*n) */
void fft_momentum(double *fft_momentum, double dx, int n) /* output momentum array, timestep, number of timesteps */
{
	/* Declaration of variables */
    	int i;
	/* Fill the output aaray with frequencies */
	for (i = 0; i < n; i++)	{
		fft_momentum[i] = planck/dx/n;
	}
}

/* Makes a momentum array fft_momentum with momentum interval planck/(dx*n) with a centered O momentum */
void fft_momentum_shift(double *fft_momentum, double dx, int n) /* output momentum array, timestep, number of timesteps */
{
	/* Declaration of variables */
    	int i;
	/* Fill the output aaray with shifted frequencies */
	for (i = 0; i < n; i++)	{
		if (n % 2) /*if n odd*/	{
			if (i<=(n-2)/2)	{
				fft_momentum[i] = (-(n-1)/2+i)*planck/dx/n;
			}
			else {
				fft_momentum[i] = (i-(n-1)/2)*planck/dx/n;
			}
		}
		else {
			if (i<n/2) {
				fft_momentum[i] = (-n/2+i)*planck/dx/n;
			}
			else {
				fft_momentum[i] = (i-n/2)*planck/dx/n;
			}
		}
	}
}
