/*
verlet.h

Header file for verlet.c


*/

#ifndef _verlet_h
#define _verlet_h

extern void calc_acc(double a, double *u, double *m, double kappa, double alpha, int size_of_u);

extern void velocity_verlet(int n_timesteps, int n_particles, double m[n_particles], double v[][n_particles],
                            double q[][n_particles], double dt, double kappa, double alpha);

extern void lattice_velocity_verlet(int n_timesteps, double cell_length, int n_particles, double m[n_particles], double v[n_particles][3],
                                    double q[n_particles][3], double T[n_timesteps], double V[n_timesteps], double E[n_timesteps], double dt);

extern void lattice_velocity_verlet_scaled(int n_timesteps, double a0, double Nc,
  int n_particles,
  double m[n_particles], double v[n_particles][3], double q[n_particles][3],
  double T[n_timesteps], double V[n_timesteps], double E[n_timesteps],
  double dt, unsigned int enable_scaling,
  double temp_inter, double temp_eq, double pressure_eq,
  double t_eq, double Temp[n_timesteps], double Pressure[n_timesteps],
  double a0_ev[n_timesteps], double time_array[n_timesteps+1]);

// Different from lattice_velocity_verlet_scaled, this one has no temp_inter.
extern void verlet_equilibration(int n_timesteps, double *a0, double Nc,
  int n_particles,
  double m[n_particles], double v[n_particles][3], double q[n_particles][3],
  double T[n_timesteps], double V[n_timesteps], double E[n_timesteps],
  double dt, unsigned int enable_scaling,
  double temp_inter, double temp_eq, double pressure_eq,
  double t_eq, double Temp[n_timesteps], double Pressure[n_timesteps],
  double a0_ev[n_timesteps], double time_array[n_timesteps+1]);

// Different from verlet_equilibration, this one has no pressure equilibration!
extern void verlet_thermostat(int n_timesteps, double a0, double Nc,
  int n_particles,
  double m[n_particles], double v[n_particles][3], double q[n_particles][3],
  double T[n_timesteps], double V[n_timesteps], double E[n_timesteps],
  double dt, unsigned int enable_scaling,
  double temp_inter, double temp_eq,
  double t_eq, double Temp[n_timesteps], double Pressure[n_timesteps],
  double a0_ev[n_timesteps], double time_array[n_timesteps+1]);

#endif
