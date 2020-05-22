#include "input.h"
#include "fluorescence.h"

double**
rate_calc(unsigned int N, double **eig, double** wij, Parameters *p)
{
  unsigned int i, j, k;
  double **kij;
  void *vptr;
  kij = calloc(N, sizeof(double));
  for (i = 0; i < N; i++) {
    kij[i] = calloc(N, sizeof(double));
  }

  double elem = 0.0;
  unsigned short print_kij = 0;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      for (k = 0; k < N; k++) {
        vptr = &p[k];
        /* 100 cm^-1 = 53 fs^-1 = 53000 ps^-1 */
        /* this is the wrong way up i think */
        elem = (100.0/53000.0) * (pow(eig[i][k], 2.) * pow(eig[j][k], 2.) *
          p[k].cw(fabs(wij[i][j]), vptr));
        kij[i][j] += elem;
        if (print_kij) {
          /* fprintf(stdout, "%d %d %d %10.6e %10.6e %10.6e %10.6e %10.6e ", i, j, k, */
          /*     wij[i][j], p[k].cw(wij[i][j], vptr), eig[i][k], eig[j][k], elem); */
	}
      }
      if (print_kij) {
        fprintf(stdout, "\n%d %d %13.10e %13.10e", i, j, wij[i][j],
          kij[i][j]);
        fprintf(stdout, "\n");
      }
    }
  }
  return kij;
}

int
jacobian (double t, const double y[], double *dfdy,
          double dfdt[], void *params)
{
  (void)(t);
  double elem;
  /* get parameter(s) from params_ptr; here, just a double */
  ode_params *p = (ode_params *) params;

  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, p->N, p->N);

  gsl_matrix *m_ptr = &dfdy_mat.matrix;	/* m_ptr points to the matrix */

  /* fill the Jacobian matrix as shown */
  unsigned short print_elem = 0;
  for (unsigned int i = 0; i < p->N; i++) {
    for (unsigned int j = 0; j < p->N; j++) {
      if (i == j) {
        /* gsl_matrix_set (m_ptr, i, j, */
        /*     (p->kij[i][j] * y[j]) - (p->gamma[i] * y[i]) + p->chiw[i]); */
        elem = (p->kij[i][j]) - (1. / (p->gamma[i] * 1000));
        if (print_elem) {
          fprintf(stdout, "%d %d %8.6f ", i, j, elem);
        }
        gsl_matrix_set (m_ptr, i, j, elem);
      } else {
        /* gsl_matrix_set (m_ptr, i, j, (p->kij[i][j] * y[j])); */
        elem = (p->kij[i][j]);
        if (print_elem) {
          fprintf(stdout, "%d %d %8.6f ", i, j, elem);
        }
        gsl_matrix_set (m_ptr, i, j, elem);
      }
    }
    if (print_elem) {
      fprintf(stdout, "\n");
    }
    dfdt[i] = 0.0; /* set explicit t dependence of f[i] */
  }

  return GSL_SUCCESS;		/* GSL_SUCCESS defined in gsl/errno.h as 0 */
}

int
odefunc(double x, const double *y, double *f, void *params)
{
  (void)(x); /* suppress unused parameter warning */
  unsigned int i, j;
  ode_params *p = (ode_params *) params;
  for (i = 0; i < p->N; i++) {
    for (j = 0; j < p->N; j++) {
      if (i == j) {
        /* f[i] += p->kij[i][i] * y[i] - (p->gamma[i] * y[i]) + p->chiw[i]; */
        f[i] += p->kij[i][i] * y[i] 
             - (y[i] * (1. / (1000 * p->gamma[i])))
             + p->chiw[i];
      } else {
        f[i] += p->kij[i][j] * y[j];
      }
    }
    /* fprintf(stdout, "i = %d; f[i] = %+f; y[i] = %+f ", i, f[i], y[i]); */
    /* fprintf(stdout, "- y[i] * (1. / 1000 gamma_i) = %+f\n", */
        /* - y[i] * (1. / (1000 * p->gamma[i]))); */
  }
  return GSL_SUCCESS;
}

double*
bcs (unsigned int N, double* eigvals)
{
  /* calculate t = 0 boltzmannised exciton populations */
  double temp = 300.0; /* probably shouldn't hardcode this lmao */
  double beta = 1./temp;
  double sum = 0.0;
  double *populations = calloc(N, sizeof(double));
  for (unsigned int i = 0; i < N; i++) {
    sum += exp(-1. * beta * eigvals[i]);
  }
  for (unsigned int i = 0; i < N; i++) {
    populations[i] = exp(-1. * beta * eigvals[i]) / sum;
  }
  return populations;
}

double**
jacmat (ode_params *p)
{
  double **Jij = calloc(p->N, sizeof(double));
  for (unsigned int i = 0; i < p->N; i++) {
    Jij[i] = calloc(p->N, sizeof(double));
  }
  for (unsigned int i = 0; i < p->N; i++) {
    for (unsigned int j = 0; j < p->N; j++) {
      if (i == j) {
        Jij[i][j] = p->kij[i][j] - (1. / (p->gamma[i] * 1000));
      } else {
        Jij[i][j] = p->kij[i][j];
      }
    }
  }
  return Jij;
}

double
trapezoid(double *f, unsigned int n)
{
    double dx = 1./n;
    double complex sum;

    sum = 0.5 * dx * (f[0] + f[n - 1]);
    for (unsigned int i = 1; i < n - 1; i++) {
	sum += dx * f[i];
    }
    return sum;
}
