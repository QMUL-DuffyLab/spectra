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

gsl_vector*
guess(const ss_init p, const double* boltz, const double* musq,
      unsigned const int max, unsigned const int N)
{
  gsl_vector *x = gsl_vector_alloc(N);
  fprintf(stdout, "Initial population guess:\n");
  for (unsigned i = 0; i < N; i++){
    if (p == BOLTZ) {
      gsl_vector_set(x, i, boltz[i]);
      fprintf(stdout, "%d %12.8e ", i, boltz[i]);
    } else if (p == BOLTZ_MUSQ) {
      gsl_vector_set(x, i, boltz[i] * musq[i]);
      fprintf(stdout, "%d %12.8e ", i, boltz[i] * musq[i]);
    } else if (p == CONST) {
      gsl_vector_set(x, i, 1. / N);
      fprintf(stdout, "%d %8.6f ", i, 1. / N);
    } else if (p == MAX) {
      if (i == max) {
        gsl_vector_set(x, i, 1.);
        fprintf(stdout, "%d %12.8e ", i, gsl_vector_get(x, i));
      } else {
        gsl_vector_set(x, i, 0.);
        fprintf(stdout, "%d %12.8e ", i, gsl_vector_get(x, i));
      }
    } else {
      fprintf(stdout, "ss_init given false value\n");
      exit(EXIT_FAILURE);
    }
  }
  fprintf(stdout, "\n");
  return x;
}

double*
incident(pulse p, unsigned int tau)
{
  /* returns W(\omega), the spectrum of incident light */
  double wn, diff, min;
  unsigned int min_arg;
  double *ww = calloc(tau, sizeof(double));
  double sum = 0.0; /* lorentzian / gaussian won't be normalised */
  for (unsigned int i = 0; i < tau; i++) {
    wn = (i * 2. * M_PI / tau) * (1. / TOFS);
    if (p.type == FLAT) {
      ww[i] = 1. / tau;
    }
    else if (p.type == LORENTZIAN) {
      ww[i] = (1. / (M_PI * p.width)) * pow(p.width, 2.) / 
              (pow(wn - p.centre, 2.) + pow(p.width, 2.));
      sum += ww[i];
    }
    else if (p.type == GAUSSIAN) {
      ww[i] = (1. / (p.width * sqrt(2. * M_PI))) *
              exp(-0.5 * pow(fabs(wn - p.centre) / p.width, 2.));
      sum += ww[i];
    }
    else if (p.type == DELTA) {
      /* the centre of the delta function probably
       * won't coincide exactly with any wn; hence put
       * the pulse at min(centre - wn) */
      diff = fabs(p.centre - wn);
      if (i == 0) {
        min = diff;
        min_arg = i;
      }
      if (diff <= min) {
        min = diff;
        min_arg = i;
      }
    }
  }
  if (p.type == DELTA) {
    /* if p.type is DELTA we haven't changed any element of
     * ww, so they're still all 0. hence we can just set the
     * closest bin to the centre of the delta function to 1
     * and be sure it's still normalised. */
    ww[min_arg] = 1.;
  }
  if (sum != 0.0) {
    for (unsigned int i = 0; i < tau; i++) {
      ww[i] /= sum;
    }
  }
  return ww;
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
  /* gsl_vector_add(f, chiw_gsl); */

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
