
#include <math.h>
#include <stdio.h>

double calc_temp(double K, double N){
    double kb = 8.617333262145e-5;

    double T_t = (K/N)*(2.0/(3.0*kb));
    return T_t;
}

double calc_alpha_t(double T_eq, double t, double tau_t, double dt, double K, double N){
    double kb = 8.617333262145e-5;
    
    double T_t = (K/N)*(2.0/(3.0*kb));
    double alpha_t = 1 + ((2*dt)/tau_t) * (T_eq - T_t)/T_t;

    return alpha_t;
}

double calc_alpha_p(double P_eq, double t, double tau_p, double dt, double kappa, double K, double V, double N, double cell_length, double num_cells){

    double volume = num_cells*pow(cell_length, 3.0);
    double P_t = 1.0/(3.0*volume) * 2.0*(K+V);
    double alpha_p = 1 -kappa*(dt/tau_p)*(P_eq - P_t);

    return alpha_p;
}
