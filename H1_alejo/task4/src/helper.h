/*
Helper.h

Created by Anders Lindman on 2013-03-15.
*/

#ifndef _helper_h
#define _helper_h

extern  void arange(double *array, double start, int n_points, double step);

extern void write_energies_file(char *fname, double *time_array, int n_timesteps, double T[n_timesteps], double V[n_timesteps], double E[n_timesteps]);

extern void write_temperatures_file(char *fname, double *time_array, int n_timesteps, double Temp[n_timesteps]);

extern void calc_time_average(int n_timesteps, double O[], double O_exp[]);

extern double calc_eq_average(int n_timesteps, double A[], int i_0);

#endif
