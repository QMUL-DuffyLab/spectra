#include "steady_state.h"

unsigned short int
pop_converge(double *y, double *yprev, unsigned int N, double thresh)
{
  unsigned short int result = 0;
  double *diff = calloc(N, sizeof(double));
  unsigned int conv = 0;
  for (unsigned int i = 0; i < N; i++) {
    diff[i] = fabs(yprev[i] - y[i]);
    if (diff[i] < thresh) {
      conv++;
    }
  }
  if (conv == N) {
    result = 1;
  }
  return result;
}

int
pop_steady_f
(const gsl_vector *x, void *params, gsl_vector *f)
{
  ode_params *p = (ode_params *)params;
  for (unsigned int i = 0; i < p->N; i++) {
    double elem = 0.0;
    for (unsigned int j = 0; j < p->N; j++) {
      if (i == j) {
        elem += -1. * p->rates[i] * gsl_vector_get(x, i)
          + p->chiw[i];
      }
      elem += p->kij[i][j] * gsl_vector_get(x, j);
    }
    gsl_vector_set(f, i, elem);
  }
  return GSL_SUCCESS;
}

int
pop_steady_df
(const gsl_vector *x, void *params, gsl_matrix *J)
{
  ode_params *p = (ode_params *)params;
  for (unsigned int i = 0; i < p->N; i++) {
    for (unsigned int j = 0; j < p->N; j++) {
      if (i == j) {
        gsl_matrix_set(J, i, j, - p->rates[i]);
      }
      gsl_matrix_set(J, i, j, p->kij[i][j]);
    }
  }
  return GSL_SUCCESS;
}

int
pop_steady_fdf
(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J)
{
  ode_params *p = (ode_params *)params;
  for (unsigned int i = 0; i < p->N; i++) {
    double elem = 0.0;
    for (unsigned int j = 0; j < p->N; j++) {
      if (i == j) {
        elem += -1. * p->rates[i] * gsl_vector_get(x, i)
        + p->chiw[i];
        /* elem -= (1. / (1000 * p->gamma[i])) * gsl_vector_get(x, i); */
        gsl_matrix_set(J, i, j, -1. * p->rates[i]);
      }
      elem += p->kij[i][j] * gsl_vector_get(x, j);
      gsl_matrix_set(J, i, j, p->kij[i][j]);
    }
    gsl_vector_set(f, i, elem);
  }
  return GSL_SUCCESS;
}
