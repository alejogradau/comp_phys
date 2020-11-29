/*
 H1equilibration.h

Header file for H1equilibration.c
 

*/

#ifndef _h1equilibration_h
#define _h1equilibration_h

extern double calc_temp(double K, double N);

extern double calc_pressure(double Nc, double a0, double T, double virial);

extern double calc_alpha_t(double T_t, double T_eq, double tau_t, double dt);

extern double calc_alpha_p(double P_eq, double tau_p, double dt, double kappa, double T, double virial, double a0, double Nc);

#endif
