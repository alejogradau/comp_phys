/******************************************************************************
 * verlet
 ******************************************************************************
 * Routine that runs the velocity verlet algorithm
 * Use as template to construct your program!
 *
 * Compile me as:
 * clang E1code4.c -o code4 -lm
 */

/******************************************************************************
 * Includes
 *****************************************************************************/

#include <stdint.h>  //uint64_t
#include <stdio.h>  //uint64_t
#include <math.h>    //pow function
#include "H1potential.h"
#include "H1equilibration.h"

/*
 * Calculate the acceleration
 * @a - vector that is filled with acceleration
 * @u - vector with the current positions
 * @m - vector with masses
 * @kappa - Spring constant
 * @alpha - Anharmonic coupling constant
 * @size_of_u - the size of the position, acceleration and mass array
 */
void calc_acc(double *a, double *u, double *m, double kappa, double alpha, int size_of_u)
{
    /* Declaration of variables */
    int i;

    /* Calculating the acceleration on the boundaries */
    a[0] = (kappa*(- 2*u[0] + u[1]) + alpha*(pow(u[1] - u[0], 2)))/m[0];
    a[size_of_u - 1] = (kappa*(u[size_of_u - 2] - 2*u[size_of_u - 1]) + alpha*(-pow(u[size_of_u - 1] - u[size_of_u - 2], 2)))/m[size_of_u - 1];

    /* Calculating the acceleration of the inner points */
    for (i = 1; i < size_of_u - 1; i++){
        a[i] = (kappa*(u[i - 1] - 2*u[i] + u[i + 1]) + alpha*(pow(u[i + 1] - u[i], 2) - pow(u[i] - u[i-1], 2)))/m[i];
    }
}

/*
 * Calculate the acceleration
 * @a - vector that is filled with acceleration
 * @u - vector with the current positions
 * @m - vector with masses
 * @kappa - Spring constant
 * @alpha - Anharmonic coupling constant
 * @size_of_u - the size of the position, acceleration and mass array
 */
//void calc_lattice_acc(double *a, double *u, double *m, int size_of_u)
//{
//    get_forces_AL(double a, double u, double cell_length, int nbr_atoms);
//}

/*
 * Perform the velocity verlet alogrithm
 * @n_timesteps - The number of time steps to be performed
 * @n_particles - number of particles in the system
 * @v - array of velocity (Empty allocated array) : sizeof(v) = n_particles
 * @q_n - position of the n'th atom : sizeof(q_n) = n_timesteps+1
 * @dt - timestep
 * @m - vector with masses of atoms sizeof(n_particles)
 * @kappa - Spring constant
 */
void velocity_verlet(int n_timesteps, int n_particles, double m[n_particles], double v[][n_particles],
		     double q[][n_particles], double dt, double kappa, double alpha)
{
    double a[n_particles];

    calc_acc(a, q[0], m, kappa, alpha, n_particles);
    for (int i = 1; i < n_timesteps + 1; i++) {
//        /* v(t+dt/2) */
        for (int j = 0; j < n_particles; j++) {
            v[i][j] = v[i-1][j] + dt * 0.5 * a[j];
        }

        /* q(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            q[i][j] = q[i-1][j] + dt * v[i][j];
        }

        /* a(t+dt) */
        calc_acc(a, q[i], m, kappa, alpha, n_particles);


        /* v(t+dt) and T(t*dt)*/
        for (int j = 0; j < n_particles; j++) {
            v[i][j] += dt * 0.5 * a[j];
        }
    }
}

/*
 * Perform the velocity verlet alogrithm with an equilibration period
 * @n_timesteps - The number of time steps to be performed
 * @a0 - Lattice parameter/length of unit cell
 * @Nc - Number of unit cells the super cell has in each direction
 * @n_particles - number of particles in the system
 * @m[n_particles] - vector with masses of atoms sizeof(n_particles)
 * @v - array of velocity (Empty allocated array) : sizeof(v) = n_particles
 * @q -
 * @T[n_particles] -
 * @V[n_particles] -
 * @E[n_particles] -
 * @dt - timestep
 * @enable_scaling -
 * @temp_eq -
 * @pressure_eq -
 * @t_eq - Duration of equilibration routine. Given in picoseconds.
 * @Temp[n_timesteps] -
 * @Pressure[n_timesteps] -
 * @a0_ev[n_timesteps] -
 * @time_array[n_timesteps+1] -
 * @kappa - Spring constant
 */

void lattice_velocity_verlet_scaled(int n_timesteps, double a0, double Nc,
  int n_particles,
  double m[n_particles], double v[n_particles][3], double q[n_particles][3],
  double T[n_timesteps], double V[n_timesteps], double E[n_timesteps],
  double dt, unsigned int enable_scaling,
  double temp_inter, double temp_eq, double pressure_eq,
  double t_eq, double Temp[n_timesteps], double Pressure[n_timesteps],
  double a0_ev[n_timesteps], double time_array[n_timesteps+1])
{
    double a[n_particles][3];
    double alpha_t, alpha_p, volume, virial;
    double L = a0*Nc;
    //at 300K, http://www.knowledgedoor.com/2/elements_handbook/aluminum.html, converted to bars
    double kappa = 1.385e-6;
    double final_eq_timestep = (t_eq/dt)+1;

    FILE *fp = fopen("./out/positions.csv", "w");
    fprintf(fp, "time, x1, y1, x1, x2, y2, z2, x3, y3, z3\n");

    //Fills array with forces every particle experiences
    get_forces_AL(a, q, L, n_particles);

    //Divides forces by each particle's mass to obtain their respective accelerations
    for(int i = 0; i < n_particles; i++){
        a[i][0] /= m[i];
        a[i][1] /= m[i];
        a[i][2] /= m[i];
    }

    //Save values for the initial conditions
    a0_ev[0] = a0;
    V[0] = get_energy_AL(q, L, n_particles);
    E[0] = T[0] + V[0];
    virial = get_virial_AL(q, L, n_particles);
    volume = calc_volume(Nc, a0);
    Temp[0] = calc_temp(T[0], n_particles);
    Pressure[0] = calc_pressure(volume, T[0], virial);

    fprintf(fp, "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", time_array[0],
    q[0][0], q[0][1], q[0][2],
    q[127][0], q[127][1], q[127][2],
    q[255][0], q[255][1], q[255][2]);

    //One verlet step
    for (int i = 1; i < n_timesteps + 1; i++) {
        /* v(t+dt/2) */
        for (int j = 0; j < n_particles; j++) {
            v[j][0] += dt * 0.5 * a[j][0];
            v[j][1] += dt * 0.5 * a[j][1];
            v[j][2] += dt * 0.5 * a[j][2];
        }

        /* q(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            q[j][0] += dt * v[j][0];
            q[j][1] += dt * v[j][1];
            q[j][2] += dt * v[j][2];
        }

        /* a(t+dt) */
        get_forces_AL(a, q, L, n_particles);

        for(int j = 0; j < n_particles; j++){
            a[j][0] /= m[j];
            a[j][1] /= m[j];
            a[j][2] /= m[j];
        }

        /* v(t+dt) and T(t*dt)*/
        for (int j = 0; j < n_particles; j++) {
            v[j][0] += dt * 0.5 * a[j][0];
            v[j][1] += dt * 0.5 * a[j][1];
            v[j][2] += dt * 0.5 * a[j][2];

            T[i] += 0.5*(pow(v[j][0], 2) + pow(v[j][1], 2) + pow(v[j][2], 2))*m[j];
        }

        V[i] = get_energy_AL(q, L, n_particles);
        E[i] = T[i] + V[i];
        virial = get_virial_AL(q, L, n_particles);
        volume = calc_volume(Nc, a0);
        Temp[i] = calc_temp(T[i], n_particles);
        Pressure[i] = calc_pressure(volume, T[i], virial);

        a0_ev[i] = a0;

        if((enable_scaling) && (i < final_eq_timestep)){
          alpha_p = cbrt(calc_alpha_p(Pressure[i], pressure_eq,
            dt*100, dt, kappa));
          if(i < final_eq_timestep/2){
              alpha_t = sqrt(calc_alpha_t(Temp[i], temp_inter,  dt*100, dt));
          }
          if( (i >= final_eq_timestep/2)){
              alpha_t = sqrt(calc_alpha_t(Temp[i], temp_eq,  dt*100, dt));
          }

        for (int j = 0; j < n_particles; j++) {
            for (int k = 0; k < 3; k++) {
                v[j][k] *= alpha_t;
                q[j][k] *= alpha_p;

                a0 *= alpha_p;
                L = a0*Nc;
            }
        }
    }
        fprintf(fp, "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", time_array[i], q[0][0], q[0][1], q[0][2], q[107][0], q[107][1], q[107][2], q[255][0], q[255][1], q[255][2]);
    }
    fclose(fp);
}
