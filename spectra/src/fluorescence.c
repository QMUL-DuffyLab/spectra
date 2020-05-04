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

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      for (k = 0; k < N; k++) {
        vptr = &p[k];
        kij[i][j] += (1./TOCM1) * (pow(eig[i][k], 2.) * pow(eig[j][k], 2.) *
          p[k].cw(wij[i][j], vptr));
        /* fprintf(stdout, "%d %d %d %d %13.10e %13.10e ", i, j, k, p[k].ligand, */
        /*     wij[i][j], p[k].cw(wij[i][j], vptr)); */
      }
      fprintf(stdout, "%d %d %13.10e %13.10e", i, j, wij[i][j],
          kij[i][j]);
    fprintf(stdout, "\n");
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
  (void)(x);
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
    fprintf(stdout, "i = %d; f[i] = %+f; y[i] = %+f\n", i, f[i], y[i]);
  }
  return GSL_SUCCESS;
}

