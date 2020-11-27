/*
 H1equilibration.h

Header file for H1equilibration.c
 

*/

#ifndef _h1equilibration_h
#define _h1equilibration_h

extern double calc_temp(double K, double N);

extern double calc_alpha_t(double T_eq, double t, double tau_t, double dt, double K, double N);

extern double calc_alpha_p(double P_eq, double t, double tau_p, double dt, double kappa, double K, double V, double N, double cell_length, double num_cells);

#endif
