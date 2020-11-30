/*
 H1equilibration.h

Header file for H1equilibration.c
 

*/

#ifndef _h1equilibration_h
#define _h1equilibration_h

extern double calc_volume(double Nc, double a0);

extern double calc_temp(double K, double N);

extern double calc_pressure(double volume, double T, double virial);

extern double calc_alpha_t(double T_t, double T_eq, double tau_t, double dt);

extern double calc_alpha_p(double P_t, double P_eq, double tau_p, double dt, double kappa);

extern double calc_eq_average(int n_timesteps, double A[], unsigned int i_0);

extern double calc_eq_var(int n_timesteps, double A[], double A_equilibration, unsigned int i_0);

extern double calc_cv_NVE(unsigned int n_particles, double temp_eq, double variance);

#endif
