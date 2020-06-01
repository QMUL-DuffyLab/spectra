#include "fluorescence.h"
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

int pop_steady_f (gsl_vector *x, void *params, gsl_vector *f);
int pop_steady_df (gsl_vector *x, void *params, gsl_matrix *J);
int pop_steady_fdf (gsl_vector *x, void *params,
                    gsl_vector *f, gsl_matrix *J);
