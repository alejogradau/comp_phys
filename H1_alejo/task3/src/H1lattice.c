/*
H1lattice.c
Program that arranges atoms on a fcc lattice. 
Created by Anders Lindman on 2013-03-15.
*/

#include <stdlib.h>
#include <time.h>

/* Function takes a matrix of size [4*N*N*N][3] as input and stores a fcc lattice in it. N is the number of unit cells in each dimension and lattice_param is the lattice parameter. */
void init_fcc(double positions[][3], int N, double lattice_param)
{
    int i, j, k;
    int xor_value;
    
    for (i = 0; i < 2 * N; i++){
        for (j = 0; j < 2 * N; j++){
            for (k = 0; k < N; k++){
                if (j % 2 == i % 2){
                    xor_value = 0;
                }
                else {
                    xor_value = 1;
                }
                positions[i * N * 2 * N + j * N + k][0] = lattice_param * (0.5 * xor_value + k);
                positions[i * N * 2 * N + j * N + k][1] = lattice_param * (j * 0.5);
                positions[i * N * 2 * N + j * N + k][2] = lattice_param * (i * 0.5);
            }
        }
    }
}


void deviate_fcc(double positions[][3], int N, double lattice_param)
{
    /*
     Generates a uniform random displacement between -6.5% and 6.5% the passed lattice parameter for each coordinate in the array.
    */
     srand(time(NULL));
     double random_value;
    
    for (int i = 0; i < N; i++){
        for (int j = 0; j < 3; j++){
            random_value = ((double) rand() / (double) RAND_MAX  * 2 - 1) * 0.065 * lattice_param;
            positions[i][j] += random_value;
        }
    }
}
