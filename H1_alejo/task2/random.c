#include <stdlib.h>
#include <stdio.h>
#include <time.h>

double deviate_fcc(double lattice_param)
{
    /*
     * Generates a uniform random displacement between -6.5% and 6.5% the passed
     * lattice parameter for each coordinate in the array.
     */
    double random_value;

    random_value = ((2 * ((double) rand() / (double) RAND_MAX)) - 1)
                        * 0.065 * lattice_param;
    return random_value;
}

int main()
{
  srand(time(NULL));
  double lattice_param = 1.0;
  for (int i = 0; i < 10; i++)
  {
    double randm = deviate_fcc(lattice_param);
    printf("randomly: %f\n", randm);
  }
}
