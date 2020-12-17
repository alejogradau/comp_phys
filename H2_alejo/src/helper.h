#ifndef _helper_h
#define _helper_h

/******************************************************************************
 * Helper functions header
 *****************************************************************************/
extern unsigned int mc_integration_metropolis(unsigned int N,
              double alpha, double burn_factor, double d,
              double *conf_m, double *pos_large, int run);

extern void variational_mc(unsigned int N, double alpha, double burn_factor,
                           double d, int n_p, double beta);

extern double gen_trial_change(double accepted, double rn, double d);

extern double vector_magnitude(double x, double y, double z);

extern double relative_prob(double x1_m, double y1_m, double z1_m,
                            double x2_m, double y2_m, double z2_m,
                            double x1_t, double y1_t, double z1_t,
                            double x2_t, double y2_t, double z2_t,
                            double alpha);

extern double local_energy(double x1, double x2,
                           double y1, double y2,
                           double z1, double z2, double alpha);

extern void calc_autocorrelation(unsigned long length,
                                 double time_series[length],
                                 double autocorrelation[length/2]);

extern void file_to_array(char *fname, unsigned int length,
                          double data[length]);

extern void calc_block_average(unsigned long length, double time_series[length],
                               double block_average[length/2]);

extern void slice_array(unsigned int N, double pos_large[],
                        unsigned int n_accepted, double *pos);

extern void int_array_to_file(char *fname, int length, int *array);

extern void array_to_file(char *fname, int length, double *array);

extern double map_to_int(unsigned int n_accepted, double pos[n_accepted],
                         int *int_pos, int n_bins);

extern void build_histogram(unsigned int numvalues, int *values);

extern void radial_density_file(char *fname, double bin_size, int n_bins,
                                double Z);

#endif
