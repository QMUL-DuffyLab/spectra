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
        elem = (1./TOCM1) * (pow(eig[i][k], 2.) * pow(eig[j][k], 2.) *
          p[k].cw(wij[i][j], vptr));
        kij[i][j] += elem;
        if (print_kij) {
          fprintf(stdout, "%d %d %d %10.6e %10.6e %10.6e %10.6e %10.6e ", i, j, k,
              wij[i][j], p[k].cw(wij[i][j], vptr), eig[i][k], eig[j][k], elem);
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

double**
jacobian_calc(unsigned int N, double **kij, double *gamma)
{
  unsigned int i, j;
  double **Jij;
  Jij = calloc(N, sizeof(double));
  for (i = 0; i < N; i++) {
    Jij[i] = calloc(N, sizeof(double));
  }

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      if (i == j) {
        Jij[i][j] = kij[i][j] - gamma[i];
      } else {
        Jij[i][j] = kij[i][j];
      }
    }
  }
  return Jij;
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
        f[i] -= p->gamma[i] * y[i];
      } else {
        f[i] += p->kij[i][j] * y[j];
      }
    }
    /* fprintf(stdout, "i = %d; f[i] = %+f; y[i] = %+f\n", i, f[i], y[i]); */
  }
  return GSL_SUCCESS;
}

