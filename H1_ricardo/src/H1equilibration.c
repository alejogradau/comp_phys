
#include <math.h>

double calc_volume(double Nc, double a0){
    const double num_cells = pow(Nc, 3.0); //Nc cells in each direction (x, y, z)
    return num_cells*pow(a0, 3.0);
}

double calc_temp(double K, double N){
    const double kb = 8.617333262145e-5;
    double factor = 2.0/(3.0*N*kb);
    
    return K * factor;
}

double calc_pressure(double volume, double T, double virial){
    double pressure = (1.0/volume) * (2.0/3.0) * (T-virial);
    return pressure/624e-7; //Conversion to bars
    
}

double calc_alpha_t(double T_t, double T_eq, double tau_t, double dt){
    return 1 + ((2*dt)/tau_t) * (T_eq - T_t)/T_t;
}

double calc_alpha_p(double P_t, double P_eq, double tau_p, double dt, double kappa){
    return 1 - kappa*(dt/tau_p)*(P_eq - P_t);
}
