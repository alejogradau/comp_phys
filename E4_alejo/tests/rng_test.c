/*
 * Testing how the gsl Gaussian random number generator works.
 * KEY: this code generates ONLY ONE rn, and prints it 5 times.
 * Compile me as:
 * clang rng.c -o rng -lgsl -lgslcblas
 */

/******************************************************************************
 * Includes
 *****************************************************************************/
#include <stdio.h>   //fopen, fprintf
#include <stdlib.h>  //malloc
#include <stdint.h>  //uint64_t
#include <math.h>    //pow, fabs function (absolute values double)
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>  // Gaussian distribution

int main()
{
  // Instantiate the random number generator
  const gsl_rng_type * T;
  gsl_rng * r;
  double sigma = 1.0;

  /* create a generator chosen by the
     environment variable GSL_RNG_TYPE */
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  double rand_num = gsl_ran_gaussian(r, sigma);

  for (int i = 0; i < 5; i++)
  {
    printf("%f \n", rand_num);
  }
}
