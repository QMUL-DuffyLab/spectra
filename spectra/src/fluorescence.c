#include "fluorescence.h"

double**
rate_calc(unsigned int N, double **eig, 
          double *eigvals, double** wij, Parameters *p)
{
  /* this function returns the transfer rates between states -
   * it's *not* the transfer matrix! separate function for that */
  unsigned int i, j, k;
  double **kij;
  void *vptr;
  kij = calloc(N, sizeof(double*));
  for (i = 0; i < N; i++) {
    kij[i] = calloc(N, sizeof(double));
  }

  double elem = 0.0;
  /* T = 300K here lol un hardcode at some point! */
  double beta = 1. / 300.0;
  unsigned short print_kij = 0;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      if (i == j) {
        continue;
      }
      for (k = 0; k < N; k++) {
        vptr = &p[k];
        /* 100 cm^-1 = 53 fs^-1 = 53000 ps^-1 : (100. / 53000.)?? */
        double cmperps = 2 * M_PI * CMS * 100 * 1E-12;
        elem = cmperps * (pow(eig[i][k], 2.) * pow(eig[j][k], 2.) *
          p[k].cn((wij[i][j]), vptr));
        kij[i][j] += elem;
      }
      if (print_kij) {
        fprintf(stdout, "%8.6e ", kij[i][j]);
      }
    }
    if (print_kij) {
      fprintf(stdout, "\n");
    }
  }
  return kij;
}

double*
relaxation_rates
(unsigned int N, double *gamma, double **kij)
{
  /* take vector of lifetimes in ns,
   * output relaxation rates in ps^{-1} */
  double *res = calloc(N, sizeof(double));
  for (unsigned int i = 0; i < N; i++) {
    res[i] = (1. / (1000 * gamma[i]));
  }
  return res;
}

double**
transfer_matrix
(unsigned int N, double* relax, double** kij)
{
  /* returns the N x N transfer matrix.
   * this can be optimised a lot later */
  unsigned int i, j, k;
  double **Tij;
  unsigned short print_Tij = 1;
  Tij = calloc(N, sizeof(double*));
  for (i = 0; i < N; i++) {
    Tij[i] = calloc(N, sizeof(double));
  }
  if (print_Tij) {
    fprintf(stdout, "\n---------------\nTRANSFER MATRIX"
                    "\n---------------\n\n");
  }
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      if (i == j) {
        Tij[i][j] = -1. * relax[i];
        for (k = 0; k < N; k++) {
          /* kij will be 0 on the diagonal so we
           * don't need to worry about excluding it */
          Tij[i][j] += (-1. * kij[i][k]);
        }
      } else {
        Tij[i][j] = kij[j][i]; /* incoming rate! */
      }
      if (print_Tij) {
        fprintf(stdout, "%8.6e ", Tij[i][j]);
      }
    }
    if (print_Tij) {
      fprintf(stdout, "\n");
    }
  }
  return Tij;
}

int
jacobian (double t, const double y[], double *dfdy,
          double dfdt[], void *params)
{
  /* in this and odefunc we should be able to use the same
   * transfer matrix as for the steady state case - I
   * deliberately didn't add the incident light source term */
  (void)(t); (void)(y);
  ode_params *p = (ode_params *) params;

  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, p->N, p->N);

  gsl_matrix *m_ptr = &dfdy_mat.matrix;	/* m_ptr points to the matrix */

  unsigned short print_elem = 0;
  for (unsigned int i = 0; i < p->N; i++) {
    for (unsigned int j = 0; j < p->N; j++) {
      gsl_matrix_set (m_ptr, i, j, p->Tij[i][j]);
      if (print_elem) {
        fprintf(stdout, "%d %d %8.6f ", i, j, p->Tij[i][j]);
      }
    }
    if (print_elem) {
      fprintf(stdout, "\n");
    }
    dfdt[i] = 0.0; /* explicit t dependence of f[i] */
  }
  return GSL_SUCCESS;
}

int
odefunc(double x, const double *y, double *f, void *params)
{
  (void)(x); /* suppress unused parameter warning */
  unsigned int i, j;
  ode_params *p = (ode_params *) params;
  for (i = 0; i < p->N; i++) {
    f[i] = 0.; /* lol otherwise it just diverges */
    for (j = 0; j < p->N; j++) {
      f[i] += p->Tij[i][j] * y[j];
    }
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
jacmat (ode_params p)
{
  double **Jij = calloc(p.N, sizeof(double));
  for (unsigned int i = 0; i < p.N; i++) {
    Jij[i] = calloc(p.N, sizeof(double));
  }
  fprintf(stdout, "\n---------------\nJacobian\n"
      "---------------\n\n");
  for (unsigned int i = 0; i < p.N; i++) {
    for (unsigned int j = 0; j < p.N; j++) {
      if (i == j) {
        Jij[i][j] = p.kij[i][j] - p.rates[i];
        fprintf(stdout, "%8.6f ", Jij[i][j]);
      } else {
        Jij[i][j] = p.kij[i][j];
        fprintf(stdout, "%8.6f ", Jij[i][j]);
      }
    }
    fprintf(stdout,"\n");
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
