#include "fluorescence.h"

gsl_matrix*
array_to_gsl_matrix(unsigned int n1, unsigned int n2, double** mat)
{
  /* take n1*n2 double array, return equivalent gsl_matrix * */
  gsl_matrix *res = gsl_matrix_alloc(n1, n2);
  for (unsigned int i = 0; i < n1; i++) {
    for (unsigned int j = 0; j < n2; j++) {
      gsl_matrix_set(res, i, j, mat[i][j]);
    }
  }
  return res;
}

void
check_detailed_balance(unsigned n, double t, double thresh,
                       double **kij, double **wij)
{
  unsigned i, j;
  double beta = 1.439 / t;
  double elem, sum;
  fprintf(stdout, "\n----------------------\n"
                    "DETAILED BALANCE CHECK\n"
                    "----------------------\n");
  sum = 0.;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      elem = fabs(kij[i][j] - (kij[j][i] *
             exp(beta * wij[i][j])));
      sum += elem;
      if (elem > thresh) {
        fprintf(stdout, "detailed balance violated for element"
            " %2d, %2d: difference %8.3e\n", i, j, elem);
        exit(EXIT_FAILURE);
      }
      fprintf(stdout, "%8.3e ", elem);
    }
    fprintf(stdout, "\n");
  }
  fprintf(stdout, "sum = %8.3e\t avg. = %8.3e\n", sum, sum / (n * n));
}

double**
rate_calc(unsigned int N, double **eig, 
          double** wij, Parameters *p)
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
  double cmperps = 2 * M_PI * CMS * 100 * 1E-12;
  unsigned short print_kij = 0;
  unsigned short print_details = 0;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      if (i == j) {
        continue;
      }
      for (k = 0; k < N; k++) {
        vptr = &p[k];
        /* [k][i/j] NOT [i/j][k] because it's the *columns* which
         * hold the normalised eigenvectors, not the rows */
        elem = cmperps * (pow(eig[k][i], 2.) * pow(eig[k][j], 2.) *
          p[k].nu * p[k].cn((wij[i][j]), vptr));
        if (print_details) {
          fprintf(stdout, "\ni j k = %2d %2d %2d:\n"
              "c_k^i^2 = %8.6e\nc_k^j^2 = %8.6e\n"
              "nu_k = %8.6e\nw_ij = %8.6e\nC_n(w_ij) = %8.6e\n"
              "product = %8.6e\n",
              i, j, k, pow(eig[k][i], 2.), pow(eig[k][j], 2.),
              p[k].nu, wij[i][j], p[k].cn(wij[i][j], vptr), elem);
        }
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
(unsigned int N, double *gamma)
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
  unsigned short print_Tij = 0;
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
        /* comment out next line to turn off 
         * relaxation to check we get Boltzmann */
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
bcs (unsigned const int N, const double* eigvals, const double T)
{
  /* calculate t = 0 boltzmannised exciton populations */
  double beta = 1.439 / T;
  double sum = 0.0;
  double *populations = calloc(N, sizeof(double));
  short print_weights = 0;
  for (unsigned int i = 0; i < N; i++) {
    if (print_weights) {
      fprintf(stdout, "%d %16.12e\n", i, exp(-1. * beta * eigvals[i]));
    }
    sum += exp(-1. * beta * eigvals[i]);
  }
  for (unsigned int i = 0; i < N; i++) {
    populations[i] = exp(-1. * beta * eigvals[i]) / sum;
  }
  return populations;
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

double
mean_excitation_lifetime(unsigned n, double **Tij, double *pop)
{
  double **Tij_vr, **Tij_vr_inv, **Tij_wr_mat, 
         *work, *work2, excite, sum;
  unsigned i, j, k;
  int status;
  Tij_vr = calloc(n, sizeof(double*)); /* right eigenvectors of T_ij */
  Tij_vr_inv = calloc(n, sizeof(double*)); /* inverse of above */
  /* inverse eigenvectors of T_ij (matrix form) */
  Tij_wr_mat = calloc(n, sizeof(double*));
  work = calloc(n, sizeof(double));
  work2 = calloc(n, sizeof(double));
  for(i = 0; i < n; i++) {
    Tij_vr[i] = calloc(n, sizeof(double));
    Tij_vr_inv[i] = calloc(n, sizeof(double));
    Tij_wr_mat[i] = calloc(n, sizeof(double));
  }
  double *Tij_wr = calloc(n, sizeof(double));

  status = eig_oop(n, Tij, Tij_vr, Tij_wr);
  /* print_matrix("T_{ij} eigenvectors", n, inv_test); */
  for (i = 0; i < n; i++) {
    Tij_wr_mat[i][i] = 1./(Tij_wr[i]);
  }

  /* now we need to invert the eigenvector matrix and so
   * C * λ^{-1} * C^{-1} * P(0). Multiply from the right! */
  status = invert_matrix_oop(n, Tij_vr, Tij_vr_inv);
  matvec(n, Tij_vr_inv, pop, work);
  matvec(n, Tij_wr_mat, work, work2);

  for (i = 0; i < n; i++) {
    sum = 0.;
    for (k = 0; k < n; k++) {
      sum += Tij_vr[i][k] * work2[k];
    }
    excite -= sum;
  }
  
  free(Tij_vr); free(Tij_vr_inv); free(Tij_wr_mat);
  free(work); free(work2);

  return excite;
}

