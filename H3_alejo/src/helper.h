#ifndef _helper_h
#define _helper_h

#include <complex.h> 

/******************************************************************************
 * Helper functions header
 *****************************************************************************/
extern double *space_grid(double *x_arr, int len_arr,
                   double x0, double dx);

extern double complex *gaussian_packet(double complex *gaussian, double *x_arr,
                                int len_x, double x0, double p0, double d);

extern double *probability_density(double *prob_density,
                                   double complex *gaussian, int len_arr);

extern double gen_trial_change(double accepted, double rn, double d);

extern double vector_magnitude(double x, double y, double z);

extern void file_to_array(char *fname, unsigned int length,
                          double data[length]);

extern void array_to_file(char *fname, int length, double *array, char *header);
#endif
