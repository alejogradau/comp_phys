
#include <math.h> //pow

/*
 * Calculates the volume of a super cell
 * @Nc - Number of unit cells in each direction (x, y, z)
 * @a0 - lattice parameter (length of unit cell)
*/
double calc_volume(double Nc, double a0){
    const double num_cells = pow(Nc, 3.0); //Nc cells in each direction (x, y, z)
    return num_cells*pow(a0, 3.0);
}

/*
 * Calculates the instant temperature of a system
 * @K - Total Kinetic Energy
 * @N - Number of particles (atoms) in the system
*/
double calc_temp(double K, double N){
    const double kb = 8.617333262145e-5;
    double factor = 2.0/(3.0*N*kb);

    return K * factor;
}

/*
 * Calculates the instant pressure of a system
 * @volume - Volume of the system
 * @T - Total Kinetic Energy
 * @virial - Virial of the system
*/
double calc_pressure(double volume, double T, double virial){
    double pressure = (1.0/volume) * ((2.0/3.0)*T + virial);
    return pressure/(624e-7); //Conversion to bars
}

/*
 * Calculates the scaling coefficient to equlibrate the temperature of a system to a T_eq
 * @T_t - Instant temperature of the system at time t
 * @T_eq - Desired temperature (Equilibration Temperature)
 * @tau_t - Time constant (Defines the exponential decay rate of the temperature)
 * @dt - size of the time step, also known as Δt
*/
double calc_alpha_t(double T_t, double T_eq, double tau_t, double dt){
    return 1 + ((2*dt)/tau_t) * (T_eq - T_t)/T_t;
}

/*
 * Calculates the scaling coefficient to equlibrate the pressure of a system to a P_eq
 * @P_t - Instant pressure of the system at time t
 * @P_eq - Desired pressure (Equilibration pressure)
 * @tau_p - Time constant (Defines the exponential decay rate of the pressure)
 * @dt - size of the time step, also known as Δt
 * @kappa - isothermal compressibility
*/
double calc_alpha_p(double P_t, double P_eq, double tau_p, double dt, double kappa){
    return 1 - kappa*(dt/tau_p)*(P_eq - P_t);
}

/* Calculate the average value of the quantity A after equilibration */
double calc_eq_average(int n_timesteps, double A[], unsigned int i_0){
    double sum_of_A = 0;
    double A_equilibration = 0;
    //unsigned int tau = 100;
    //unsigned int i_0 = 2*tau;  // Equilibration
    for(int i = i_0; i < n_timesteps; i++){
            sum_of_A += A[i];
    }
    A_equilibration = sum_of_A/(n_timesteps - i_0);
    return A_equilibration;
}

 /* Calculate the variance of the average of the quantity A after
  * equilibration
  */
 double calc_eq_var(int n_timesteps, double A[], double A_equilibration,
   unsigned int i_0){
     double sum = 0;
     double variance = 0;
     for(int i = i_0; i < n_timesteps; i++){
             sum += pow((A[i] - A_equilibration),2);
     }
     variance = sum/(n_timesteps - i_0);
     return variance;
 }

 /* Calculate heat capacity for the microcanonical ensemble (constant NVE)
  * formula is taken from eq. 58 in MD_LectureNotes
  */
 double calc_cv_NVE(unsigned int n_particles, double temp_eq, double variance)
 {
   const double kb = 8.617333262145e-5;
   double factor1 = (3.0*n_particles*kb)/2.0;
   double factor2 = 2.0 / ( 3.0*n_particles*pow(kb,2)*pow(temp_eq,2) );
   double cv = factor1 / (1 - (factor2*variance));
   return cv;
 }
