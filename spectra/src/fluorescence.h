#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "input.h"
#include "invert_matrix.h"
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
void check_detailed_balance(unsigned n, double t, double thresh,
                            double **kij, double **wij);
double** rate_calc (unsigned int N, double **eig,
                    double** wij, Parameters *p);
double* relaxation_rates (unsigned int N, double* gamma);
double **transfer_matrix (unsigned int N, double* relax, double** kij);
int odefunc (double x, const double *y, double *f, void *params);
int jacobian (double t, const double y[], double *dfdy,
              double dfdt[], void *params);
double* bcs (unsigned const int N, const double* eigvals,
             const double T);
double trapezoid(double *f, unsigned int n);
void decompose_transfer_matrix(unsigned n, double **Tij,
                               double **Tij_vr, double **Tij_vr_inv,
                               double **Tij_wr);
double mean_excitation_lifetime(unsigned n, double **Tij_vr,
                                double **Tij_vr_inv,
                                double **Tij_wr, double *p0);
void population(unsigned n, double t, double *pt, double **Tij_vr,
                double **Tij_vr_inv, double **Tij_wr, double *p0);

#endif
