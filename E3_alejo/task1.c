/******************************************************************************
 * task1
 ******************************************************************************
 * Routine that performs a Monte Carlo integration
 *
 * Compile me as:
 * clang task1.c -o task1 -lm
 */

/******************************************************************************
 * Includes
 *****************************************************************************/
#include <stdio.h>   //fopen, fprintf
#include <stdlib.h>  //malloc
#include <stdint.h>  //uint64_t
#include <math.h>    //pow, rand, fabs function (absolute values double)
#include <time.h>  // time in random

/******************************************************************************
 * Helper functions
 *****************************************************************************/

/* Random number generator between [0,1] */
double randm()
{
  double random_value;
  random_value = ((double) rand() / (double) RAND_MAX);
  return random_value;
}

/* Function evaluation */
double f_eval(double xi)
{
  double eval = 0;
  eval = xi*(1-xi);
  return eval;
}

/* Calculate function average */
double f_mean(unsigned int N)
{
  double average = 0;

  for (int i = 0; i < N; i++)
  {
    double xi = randm();
    average += f_eval(xi);
  }
  average /= N;
  return average;
}

/* Calculate function^2 average */
double f2_mean(unsigned int N)
{
  double average = 0;

  for (int i = 0; i < N; i++)
  {
    double xi = randm();
    average += pow(f_eval(xi),2);
  }
  average /= N;
  return average;
}

/* Variance of the function random evaluation */
double f_var(unsigned int N)
{
  double var = 0;
  var = f2_mean(N) - pow(f_mean(N),2);
  return var;
}

/******************************************************************************
 * Main program
 *****************************************************************************/
int main()
{
  srand(time(NULL));  // Initialize seed

  /* Variables declaration */
  unsigned int N1 = 1e1;  // Number of evaluations
  double I1 = 0;
  double var1 = 0;
  unsigned int N2 = 1e2;  // Number of evaluations
  double I2 = 0;
  double var2 = 0;
  unsigned int N3 = 1e3;  // Number of evaluations
  double I3 = 0;
  double var3 = 0;
  unsigned int N4 = 1e4;  // Number of evaluations
  double I4 = 0;
  double var4 = 0;

  I1 = f_mean(N1);
  var1 = f_var(N1);
  I2 = f_mean(N2);
  var2 = f_var(N2);
  I3 = f_mean(N3);
  var3 = f_var(N3);
  I4 = f_mean(N4);
  var4 = f_var(N4);
  printf("The value of I for N = %d is: %f +/- %f\n", N1, I1, var1);
  printf("The value of I for N = %d is: %f +/- %f\n", N2, I2, var2);
  printf("The value of I for N = %d is: %f +/- %f\n", N3, I3, var3);
  printf("The value of I for N = %d is: %f +/- %f\n", N4, I4, var4);
  printf("You still need to calculate the error instead of the variance!\n");
}
