#include "steady_state.h"

unsigned short int
pop_converge(double *y, double *yprev, unsigned int N, double thresh)
{
  unsigned short int result = 0;
  double *diff = (double *)calloc(N, sizeof(double));
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
  free(diff);
  return result;
}

double*
guess(const ss_init p, const double* boltz, const double* musq,
      unsigned const int max, unsigned const int N)
{
  short print_guess = 0;
  double sum = 0.;
  double *guess = (double *)calloc(N, sizeof(double));
  if (print_guess) {
    fprintf(stdout, "Initial population guess:\n");
  }
  for (unsigned i = 0; i < N; i++){
    if (p == BOLTZ) {
      guess[i] = boltz[i];
      if (print_guess) {
        fprintf(stdout, "%d %12.8e ", i, boltz[i]);
      }
    } else if (p == MUSQ) {
      /* sum this because it won't be normalised otherwise */
      sum += musq[i];
      guess[i] = musq[i];
      if (print_guess) {
        fprintf(stdout, "%d %12.8e ", i, musq[i]);
      }
    } else if (p == BOLTZ_MUSQ) {
      /* likewise this needs to be normalised here */
      sum += boltz[i] * musq[i];
      guess[i] = boltz[i] * musq[i];
      if (print_guess) {
        fprintf(stdout, "%d %12.8e ", i, boltz[i] * musq[i]);
      }
    } else if (p == CONST) {
      guess[i] = 1. / N;
      if (print_guess) {
        fprintf(stdout, "%d %8.6f ", i, 1. / N);
      }
    } else if (p == MAX) {
      if (i == max) {
        guess[i] = 1.;
      } else {
        guess[i] = 0.;
      }
      if (print_guess) {
        fprintf(stdout, "%d %12.8e ", i, guess[i]);
      }
    } else {
      fprintf(stdout, "ss_init given false value\n");
      exit(EXIT_FAILURE);
    }
  }
  if (p == MUSQ || p == BOLTZ_MUSQ) {
    for (unsigned i = 0; i < N; i++){
      guess[i] /= sum;
    }
  }
  if (print_guess) {
    fprintf(stdout, "\n");
  }
  return guess;
}

int
pop_steady_f
(const gsl_vector *x, void *params, gsl_vector *f)
{
  unsigned int i, j;
  ode_params *p = (ode_params *)params;
  gsl_matrix *Tij_gsl = gsl_matrix_alloc(p->N, p->N);
  for (i = 0; i < p->N; i++) {
    gsl_vector_set(f, i, 0.); /* test */
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
  /* NB: steady state populations actually blow up before
   * reaching a root - something to do with the normalisation
   * of the chiw vector?? or should there be a zeroing of f? */

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

/** Return normalised steady-state populations.
 *
 */
double*
steady_state_populations(double *guess, void *params, unsigned n)
{
  double *populations = (double *)calloc(n, sizeof(double));
  const gsl_multiroot_fdfsolver_type *T;
  int status;
  unsigned i, iter = 0;
  gsl_vector *g = gsl_vector_alloc(n);
  gsl_multiroot_fdfsolver *s;
  gsl_multiroot_function_fdf FDF;
  FDF.f = &pop_steady_f;
  FDF.df = &pop_steady_df;
  FDF.fdf = &pop_steady_fdf;
  FDF.n = n;
  FDF.params = params;
  for (i = 0; i < n; i++) {
    gsl_vector_set(g, i, guess[i]);
  }

  T = gsl_multiroot_fdfsolver_hybridsj;
  s = gsl_multiroot_fdfsolver_alloc(T, n);
  gsl_multiroot_fdfsolver_set(s, &FDF, g);

  do
    {
      iter++;
      status = gsl_multiroot_fdfsolver_iterate (s);

      fprintf(stdout, "iter %d: ", iter);
      for (i = 0; i < n; i++) {
        fprintf(stdout, "x(%2d) = %8.6f ", i, gsl_vector_get(s->x, i));
      }
      fprintf(stdout, "\n");

      if (status)   /* check if solver is stuck */
        break;

      status =
        gsl_multiroot_test_residual (s->f, 1e-8);
    }
  while (status == GSL_CONTINUE && iter < 1000);

  printf ("status = %s\n", gsl_strerror (status));

  double p_i_sum = 0.0;
  /* amazingly there does not seem to be a GSL function for this */
  for (i = 0; i < n; i++) {
    p_i_sum += gsl_vector_get(s->x, i);
  }
  if (p_i_sum <= 1E-10) {
    fprintf(stdout, "Sum of steady-state populations is zero!!!\n");
    /* this stops it from e.g. normalising a vector (0, 1e-23)
     * to (0, 1) and making it look normal */
    p_i_sum = 1.;
  } else {
    for (i = 0; i < n; i++) {
      populations[i] = gsl_vector_get(s->x, i) / p_i_sum;
    }
    fprintf(stdout, "Σ_i p_i = %8.6f\n", p_i_sum);
  }

  gsl_multiroot_fdfsolver_free(s);
  gsl_vector_free(g);

  return populations;
}
