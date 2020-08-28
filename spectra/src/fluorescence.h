#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "input.h"
#include "../../lineshape/src/parameters.h"

#ifndef __FLUORESCENCE_H__
#define __FLUORESCENCE_H__

typedef struct {
  unsigned int N;
  double **kij;
  double **Tij;
  double *rates;
  double *chiw;
} ode_params;

gsl_matrix* array_to_gsl_matrix(unsigned int n1,
            unsigned int n2, double** mat);
double** rate_calc (unsigned int N, double **eig,
                    double** wij, Parameters *p);
double* relaxation_rates (unsigned int N, double* gamma);
double **transfer_matrix (unsigned int N, double* relax, double** kij);
double **final_matrix    (unsigned int N, double* relax, double** Tij);
int odefunc (double x, const double *y, double *f, void *params);
int jacobian (double t, const double y[], double *dfdy,
              double dfdt[], void *params);
double* bcs (unsigned const int N, const double* eigvals,
             const double T);
double trapezoid(double *f, unsigned int n);

#endif
