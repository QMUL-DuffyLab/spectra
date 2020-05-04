#include "../../lineshape/src/parameters.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

typedef struct {
  unsigned int N;
  double **kij;
  double *gamma;
} ode_params;

double** rate_calc(unsigned int N, double **eig, double** wij, Parameters *p);
double** jacobian_calc(unsigned int N, double **kij, double *gamma);
int odefunc(double x, const double *y, double *f, void *params);
