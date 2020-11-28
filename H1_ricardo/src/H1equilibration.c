
#include <math.h>

double calc_temp(double K, double N){
    const double kb = 8.617333262145e-5;
    double factor = 2.0/(3.0*N*kb);
    
    return K * factor;
}

double calc_pressure(double Nc, double a0, double T, double virial){
    const double num_cells = pow(Nc, 3.0);
    const double volume = num_cells*pow(a0, 3.0);
    const double volume_inv = 1.0/volume;
    
    double pressure = volume_inv * ((2.0/3.0)*T - virial);
    return pressure/624e-7; //Conversion to bars
    
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
