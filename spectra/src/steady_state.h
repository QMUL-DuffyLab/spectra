#include "fluorescence.h"
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

unsigned short int pop_converge(double *y, double *yprev, 
                                unsigned int N, double thresh);
int pop_steady_f   (const gsl_vector *x, void *params, gsl_vector *f);
int pop_steady_df  (const gsl_vector *x, void *params, gsl_matrix *J);
int pop_steady_fdf (const gsl_vector *x, void *params,
                    gsl_vector *f, gsl_matrix *J);
