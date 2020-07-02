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
  unsigned int i, j;
  ode_params *p = (ode_params *)params;
  gsl_matrix *Tij_gsl = gsl_matrix_alloc(p->N, p->N);
  for (i = 0; i < p->N; i++) {
    for (j = 0; j < p->N; j++) {
      gsl_matrix_set(Tij_gsl, i, j, p->Tij[i][j]);
    }
  }
  gsl_blas_dgemv(CblasNoTrans, 1., Tij_gsl, x, 0., f);
  gsl_vector* chiw_gsl = gsl_vector_alloc(p->N);
  for (i = 0; i < p->N; i++) {
    gsl_vector_set(chiw_gsl, i, p->chiw[i]);
  }
  gsl_vector_add(f, chiw_gsl);

  gsl_vector_free(chiw_gsl);
  gsl_matrix_free(Tij_gsl);

  return GSL_SUCCESS;
}

int
pop_steady_df
(const gsl_vector *x, void *params, gsl_matrix *J)
{
  (void)(x); /* suppress warning */
  ode_params *p = (ode_params *)params;
  for (unsigned int i = 0; i < p->N; i++) {
    for (unsigned int j = 0; j < p->N; j++) {
      gsl_matrix_set(J, i, j, p->Tij[i][j]);
    }
  }
  return GSL_SUCCESS;
}

int
pop_steady_fdf
(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J)
{
  ode_params *p = (ode_params *)params;
  pop_steady_f(x, p, f);
  pop_steady_df(x, p, J);

  return GSL_SUCCESS;
}
