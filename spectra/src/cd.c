#include "cd.h"

void
cross(double *a, double* b, double *res)
{
  res[0] = a[1] * b[2] - a[2] * b[1];
  res[1] = a[2] * b[0] - a[0] * b[2];
  res[2] = a[0] * b[1] - a[1] * b[0];
}

void
dot(double *a, double* b, double res)
{
  res = (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]);
}

void
cd(unsigned n_chl, unsigned ns, double **chiw, double **mu,
   double **eig, double **com, double *eigvals, double *cd)
{
  double v[3], rnm[3];
  double dd = 0;
  double *chiw_sums = (double *)calloc(n_chl, sizeof(double));
  /* double *g = (double *)calloc(ns, sizeof(double)); */
  for (unsigned n = 0; n < n_chl; n++) {
    double sum = 0.;
    for (unsigned m = 0; m < ns; m++) {
      sum += chiw[n][m];
    }
    chiw_sums[n] = sum;
  }


  for (unsigned n = 0; n < n_chl; n++) {
    for (unsigned m = 0; m < n_chl; m++) {
      cross(mu[n], mu[m], v);
      for (unsigned l = 0; l < 3; l++) {
        rnm[l] = com[m][l] = com[n][l];
      }
      dot(v, rnm, dd);

      for (unsigned k = 0; k < n_chl; k++) {
        for (unsigned l = 0; l < ns; l++) {
          cd[l] += (2 * PI * ((float)l / ns) * dd)
                * ((chiw[k][l] * eig[k][m] * 
                  eig[k][n] * eigvals[k])
                / chiw_sums[k]);
        }
      }

    }
  }
}
