/*
Helper.h

Created by Anders Lindman on 2013-03-15.
*/

#ifndef _helper_h
#define _helper_h

extern void arange(double *array, double start, int n_points, double step);

extern void write_energies_file(char *fname,
  double *time_array, int n_timesteps,
  double T[n_timesteps], double V[n_timesteps], double E[n_timesteps]);

extern void write_observable_file(char *fname, char* oname, double *time_array,
  int n_timesteps, double O[n_timesteps]);

extern void calc_time_average(int n_timesteps, double O[n_timesteps],
  double O_exp[n_timesteps]);

extern void distances_array(unsigned int N, double *distances[N*(N-1)/2],
  double positions[N][3]);

extern void int_distances_array(unsigned int N, double *int_distances[N*(N-1)/2],
    double positions[N][3], double dr);

#endif
