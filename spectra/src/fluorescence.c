#include "fluorescence.h"

double*
incident(pulse p, unsigned int tau)
{
  /* returns W(\omega), the (normalised) spectrum of incident light */
  double wn, diff, min;
  unsigned int min_arg;
  double *ww = (double *)calloc(tau, sizeof(double));
  double sum = 0.0; /* lorentzian / gaussian won't be normalised */
  for (unsigned int i = 0; i < tau; i++) {
    wn = (i * 2. * M_PI / tau) * (1. / TOFS);
    if (p.type == FLAT) {
      ww[i] = 1. / tau;
    }
    else if (p.type == LORENTZIAN) {
      ww[i] = (1. / (M_PI * p.width)) * 0.5 * pow(p.width, 2.) / 
              (pow(wn - p.centre, 2.) + pow(0.5 * p.width, 2.));
      sum += ww[i];
    }
    else if (p.type == GAUSSIAN) {
      double sigma = p.width / (2. * sqrt(2. * log(2.)));
      ww[i] = (1. / (sigma * sqrt(2. * M_PI))) *
              exp(-0.5 * pow(fabs(wn - p.centre) / sigma, 2.));
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

double
intensity(double w, double t, pulse p)
{
  double sigma, f;
  double gamma = 1./(2. * p.duration)
               * pow(1./cosh((t - p.t_peak)/p.duration), 2.);
  if (p.type == GAUSSIAN) {
    sigma = p.width / (2. * sqrt(2. * log(2.)));
    f = exp(pow(w - p.centre, 2.) / (2. * pow(sigma, 2.)))
      * (1. / (sigma * sqrt(2. * PI)));
    return p.amplitude * f * gamma;
  } else {
    /* not done this yet lol */
    return NAN;
  }
}

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
  double elem, sum, max;
  fprintf(stdout, "\n----------------------\n"
                    "DETAILED BALANCE CHECK\n"
                    "----------------------\n");
  sum = 0.; max = 0.;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      elem = fabs(kij[i][j] - (kij[j][i] *
             exp(beta * wij[i][j])));
      sum += elem;
      if (elem > max) {
        max = elem;
      }
      if (elem > thresh) {
        fprintf(stdout, "detailed balance violated for element"
            " %2d, %2d: difference %8.3e\n", i, j, elem);
        exit(EXIT_FAILURE);
      }
      fprintf(stdout, "%8.3e ", elem);
    }
    fprintf(stdout, "\n");
  }
  fprintf(stdout, "sum = %8.3e\t avg. = %8.3e\tmax = %8.3e\n",
          sum, sum / (n * n), max);
}

double
redfield_rate(unsigned int N, double **eig, 
          unsigned i, unsigned j, double wij, Parameters *p)
{
  /* single redfield rate - so we can build a matrix by looping
   * over and switching between forster/redfield element-wise */
  double elem = 0.0, rate = 0.;
  double cmperps = 2 * M_PI * CMS * 100 * 1E-12;
  unsigned short print_kij = 1;
  unsigned short print_details = 0;
  void *vptr;
  for (unsigned k = 0; k < N; k++) {
    vptr = &p[k];
    elem = cmperps * (pow(eig[k][i], 2.) * pow(eig[k][j], 2.) *
      p[k].nu * p[k].cn(wij, vptr));
    if (print_details) {
        fprintf(stdout, "\ni j k = %2d %2d %2d:\n"
            "c_k^i = %8.6e\tc_k^j = %8.6e\t"
            "c_k^i^2 = %8.6e\tc_k^j^2 = %8.6e\n"
            "nu_k = %8.6e\nw_ij = %8.6e\tC_n(w_ij) = %8.6e\n"
            "(1 + coth) = %8.6e\tC''(wij) = %8.6e\n"
            "product = %8.6e\n",
            i, j, k, eig[k][i], eig[k][j], 
            pow(eig[k][i], 2.), pow(eig[k][j], 2.),
            p[k].nu, wij, p[k].cn(wij, vptr),
            (1. + (1. / tanh(0.5 * wij * 1.439 / 300.))),
            (p[k].cw(wij, vptr)), elem);
    }
    rate += elem;
  }
  if (print_kij) {
    fprintf(stdout, "%8.3e\n", rate);
  }
  return rate;
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
  kij = (double **)calloc(N, sizeof(double*));
  for (i = 0; i < N; i++) {
    kij[i] = (double *)calloc(N, sizeof(double));
  }

  double elem = 0.0;
  double cmperps = 2 * M_PI * CMS * 100 * 1E-12;
  unsigned short print_kij = 1;
  unsigned short print_details = 0;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      if (i == j) {
        if (print_kij) {
          /* otherwise the diagonal doesn't get printed */
          fprintf(stdout, "%8.3e ", kij[i][j]);
        }
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
                "c_k^i = %8.6e\tc_k^j = %8.6e\t"
                "c_k^i^2 = %8.6e\tc_k^j^2 = %8.6e\n"
                "nu_k = %8.6e\nw_ij = %8.6e\tC_n(w_ij) = %8.6e\n"
                "(1 + coth) = %8.6e\tC''(wij) = %8.6e\n"
                "product = %8.6e\n",
                i, j, k, eig[k][i], eig[k][j], 
                pow(eig[k][i], 2.), pow(eig[k][j], 2.),
                p[k].nu, wij[i][j], p[k].cn(wij[i][j], vptr),
                (1. + (1. / tanh(0.5 * wij[i][j] * 1.439 / 300.))),
                (p[k].cw(wij[i][j], vptr)), elem);
        }
        kij[i][j] += elem;
      }
      if (print_kij) {
        fprintf(stdout, "%8.3e ", kij[i][j]);
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
  double *res = (double *)calloc(N, sizeof(double));
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
  Tij = (double **)calloc(N, sizeof(double*));
  for (i = 0; i < N; i++) {
    Tij[i] = (double *)calloc(N, sizeof(double));
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
  /* note that kB / (h C) in wavenumbers comes out as about 0.69
   * and 1 / (0.69 T) == 1.439 / T */
  double beta = 1.439 / T;
  double sum = 0.0;
  double *populations = (double *)calloc(N, sizeof(double));
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
    double sum;

    sum = 0.5 * dx * (f[0] + f[n - 1]);
    for (unsigned int i = 1; i < n - 1; i++) {
	sum += dx * f[i];
    }
    return sum;
}

/** Do the eigendecomposition of the transfer matrix.
 *
 * Fills the arrays Tij_vr, Tij_vr_inv, Tij_wr
 * with the right eigenvectors of Tij, the inverse of that eigenvector
 * matrix, and the right eigenvalues of Tij in matrix form (diagonal).
 *
 * These are used to calculate the mean excitation lifetime and also
 * project the populations forward as a function of time.
 *
 */
void
decompose_transfer_matrix(unsigned n, double **Tij,
                          double **Tij_vr, double **Tij_vr_inv,
                          double **Tij_wr)
{
  double *Tij_wr_vec;
  unsigned i;
  int status;
  Tij_wr_vec = (double *)calloc(n, sizeof(double));

  status = eig_oop(n, Tij, Tij_vr, Tij_wr_vec);
  if (status != 0) {
    fprintf(stdout, "Eigendecomposition of transfer matrix failed."
        " Error code %d\n", status);
  }

  for (i = 0; i < n; i++) {
    Tij_wr[i][i] = Tij_wr_vec[i];
  }

  status = invert_matrix_oop(n, Tij_vr, Tij_vr_inv);
  if (status != 0) {
    fprintf(stdout, "Inversion of transfer matrix failed."
        " Error code %d\n", status);
  }

  free(Tij_wr_vec);
}

/** Calulcate the mean excitation lifetime.
 *
 * Uses the arrays generated by the decompose_transfer_matrix function.
 * Note that that function prints out the eigenvalues along the diagonal;
 * in this case we need the _inverse_ eigenvalues.
 */
double
mean_excitation_lifetime(unsigned n, double **Tij_vr,
                         double **Tij_vr_inv,
                         double **Tij_wr, double *p0)
{
  double *work, *work2, excite, sum;
  /* need λ^{-1} here but don't want to
   * overwrite Tij_wr - we'll need it later */
  double **Tij_wr_inv;
  unsigned i, j;
  work = (double *)calloc(n, sizeof(double));
  work2 = (double *)calloc(n, sizeof(double));
  Tij_wr_inv = (double **)calloc(n, sizeof(double*));
  for (i = 0; i < n; i++) {
    Tij_wr_inv[i] = (double *)calloc(n, sizeof(double));
  }

  for (i = 0; i < n; i++) {
    Tij_wr_inv[i][i] = 1./(Tij_wr[i][i]);
  }

  /* now we need to invert the eigenvector matrix and so
   * C * λ^{-1} * C^{-1} * P(0). Multiply from the right! */
  matvec(n, Tij_vr_inv, p0, work);
  matvec(n, Tij_wr_inv, work, work2);

  excite = 0.;
  for (i = 0; i < n; i++) {
    sum = 0.;
    if (i == 0) {
      continue;
    } else if (i == 15) {
      continue;
    } else if (i == 63) {
      continue;
    } else if (i == 110) {
      continue;
    } else {
      sum = 0.;
      for (j = 0; j < n; j++) {
        if (j == 0) {
          continue;
        } else if (j == 15) {
          continue;
        } else if (j == 63) {
          continue;
        } else if (j == 110) {
          continue;
        } else {
          sum += Tij_vr[i][j] * work2[j];
          /* fprintf(stdout, "%3u %3u %10.6e %10.6e %10.6e\n", */
              /* i, j, Tij_vr[i][j], work2[j], sum); */
        }
      }
    excite -= sum;
    }
  }

  for (i = 0; i < n; i++) {
    free(Tij_wr_inv[i]);
  }
  free(work); free(work2); free(Tij_wr_inv);

  return excite;
}

/** Calculate the populations P_i(t) in array pt for given time t.
 *
 */
void
population(unsigned n, double t, double *pt, double **Tij_vr,
           double **Tij_vr_inv, double **Tij_wr, double *p0)
{
  double *work, *work2;
  double **exp_wr; /* need e^{λt} for population propagation */
  unsigned i;
  work = (double *)calloc(n, sizeof(double));
  work2 = (double *)calloc(n, sizeof(double));
  exp_wr = (double **)calloc(n, sizeof(double*));
  for (i = 0; i < n; i++) {
    exp_wr[i] = (double *)calloc(n, sizeof(double));
  }

  for (i = 0; i < n; i++) {
    exp_wr[i][i] = exp(Tij_wr[i][i] * t);
  }
  matvec(n, Tij_vr_inv, p0, work);
  matvec(n, exp_wr, work, work2);
  matvec(n, Tij_vr, work2, pt);
  for (i = 0; i < n; i++) {
    free(exp_wr[i]);
  }
  free(exp_wr); free(work); free(work2);
}
