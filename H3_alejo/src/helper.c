/******************************************************************************
 * Macro defines
 *****************************************************************************/
#define PI 3.14159
#define hbar 0.6582  // eV * fs

/******************************************************************************
 * Includes
 *****************************************************************************/
#include <time.h>    //time
#include <stdlib.h>  //srand, rand, strtol
#include <math.h>    //pow, sin, exp
#include <stdio.h>   //printf
#include <complex.h>  // Complex numbers

/* Build space grid
 * @x_arr - array to be filled in with the space grid points
 * @len_arr - the length of the array
 */
double *space_grid(double *x_arr, int len_arr,
                   double x0, double dx)
{

    int x0_index = (len_arr-1) / 2;  // Center index
    int indx_plus;
    int indx_minus;
    x_arr[x0_index] = x0;
    for(int i = 1; i <= x0_index; i++)
    {
        indx_plus = x0_index + i;
        indx_minus = x0_index - i;
        x_arr[indx_plus] = x0 + (i * dx);
        x_arr[indx_minus] = x0 - (i * dx);
    }
    return x_arr;
}

/* Gaussian wave packet
 * @gaussian - array to be filled with wave packet values
 * @x_arr - array filled with discretized space steps
 * @len_x - the length of the space array
 * @x0 - Center of the wave packet
 * @p0 - Magnitude of the mean momentum
 * @d - Width of the Gaussian
 */
double *gaussian_packet(double *gaussian, double *x_arr, int len_x,
                        double x0, double p0, double d)
{
    double d_sq = pow(d, 2.0);
    double factor = 1 / pow(PI * d_sq, 0.25);
    for(int i = 0; i < len_x; i++)
    {
        double x_i0 = x_arr[i] - x0;
        double x_i0_sq = pow(x_i0, 2.0);
        gaussian[i] = factor * exp( - x_i0_sq / (2 * d_sq) )
                             * exp( I * p0 * x_i0 / hbar );
    }
    return gaussian;
}

/* Probability density distribution
 * @prob_density - array to be filled in with the probability density distr.
 * @gaussian - array with the gaussian wave packet
 * @len_arr - the length of the wavefunction array
 */
double *probability_density(double *prob_density, double *gaussian,
                            int len_arr)
{
    double conj_gaussian;
    for(int i = 0; i < len_arr; i++)
    {
        conj_gaussian = conj(gaussian[i]);
        prob_density[i] = gaussian[i] * conj_gaussian;
    }
    return prob_density;
}

/* Calculate the magnitude of a vector, given its cartesian coordinates */
double vector_magnitude(double x, double y, double z)
{
  double x_sq = pow(x, 2.0);
  double y_sq = pow(y, 2.0);
  double z_sq = pow(z, 2.0);

  double r_sq = x_sq + y_sq + z_sq;
  double r = sqrt(r_sq);

  return r;
}

void file_to_array(char *fname, unsigned int length, double data[length])
{
    FILE *fp;
    fp = fopen(fname, "r");

    for (int i = 0; i < length; i++)
    {
        fscanf(fp, "%lf", &data[i]);
    }
    printf("Reading file with: %d lines\n", length);
    fclose(fp);
}

/* Write an array (2nd column) and its index+1 (1st column) to file
 * @ fname - File's name
 * @ length - length of the array
 * @ *array - array of double values
 */
void array_to_file(char *fname, int length, double *array, char *header)
{
  //printf("opening...\n");
  FILE *fp = fopen(fname, "w");
  //printf("opened...\n");
  fprintf(fp, "%s", header);
  //printf("header...\n");
  for (int i = 0; i < length; i++)
  {
    //printf("%d, %f\n", i, array[i]);
    fprintf(fp, "%d, %f\n", i, array[i]);
  }
  fclose(fp);
}
