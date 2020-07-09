#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "input.h"
#include "../../lineshape/src/parameters.h"

typedef struct {
  unsigned int N;
  double **kij;
  double **Tij;
  double *rates;
  double *chiw;
} ode_params;

double** rate_calc (unsigned int N, double **eig, double *eigvals,
                    double** wij, Parameters *p);
double* relaxation_rates (unsigned int N, double* gamma, double** kij);
double **transfer_matrix (unsigned int N, double* relax, double** kij);
double **final_matrix    (unsigned int N, double* relax, double** Tij);
double **jacmat (ode_params p);
int odefunc (double x, const double *y, double *f, void *params);
int jacobian (double t, const double y[], double *dfdy,
              double dfdt[], void *params);
double* bcs (unsigned int N, double* eigvals);
double trapezoid(double *f, unsigned int n);
