#include "steady_state.h"

int
pop_steady_f
(gsl_vector *x, void *params, gsl_vector *f)
{
  ode_params *p = (ode_params *)params;
  for (unsigned int i = 0; i < p->N; i++) {
    double elem = 0.0;
    for (unsigned int j = 0; j < p->N; j++) {
      if (i == j) {
        elem += (-1. / (1000 * p->gamma[i])) * gsl_vector_get(x, i)
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
(gsl_vector *x, void *params, gsl_matrix *J)
{
  ode_params *p = (ode_params *)params;
  for (unsigned int i = 0; i < p->N; i++) {
    for (unsigned int j = 0; j < p->N; j++) {
      if (i == j) {
        gsl_matrix_set(J, i, j, - (1. / (1000 * p->gamma[i])));
      }
      gsl_matrix_set(J, i, j, p->kij[i][j]);
    }
  }
  return GSL_SUCCESS;
}

int
pop_steady_fdf
(gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J)
{
  ode_params *p = (ode_params *)params;
  for (unsigned int i = 0; i < p->N; i++) {
    double elem = 0.0;
    for (unsigned int j = 0; j < p->N; j++) {
      if (i == j) {
        elem += (-1. / (1000 * p->gamma[i])) * gsl_vector_get(x, i)
        + p->chiw[i];
        /* elem -= (1. / (1000 * p->gamma[i])) * gsl_vector_get(x, i); */
        gsl_matrix_set(J, i, j, - (1. / (1000 * p->gamma[i])));
      }
      elem += p->kij[i][j] * gsl_vector_get(x, j);
      gsl_matrix_set(J, i, j, p->kij[i][j]);
    }
    gsl_vector_set(f, i, elem);
  }
  return GSL_SUCCESS;
}
