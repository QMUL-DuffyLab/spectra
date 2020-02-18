#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <fftw3.h>

int
main(int *argc, char** argv)
{
  char* input_file;
  unsigned int N, tau, i;
  double *ex, *wi, *gamma, *mu, *lambda, *integral;
  double **gi_array;

  /* need to read these in from files as well */
  wi = calloc(N, sizeof(double));
  gamma = calloc(N, sizeof(double));
  mu = calloc(N, sizeof(double));
  lambda = calloc(N, sizeof(double));

  gi_array = calloc(N, sizeof(double));
  for (i = 0; i < N; i++) {
    gi_array[i] = calloc(tau, sizeof(double));
  }
  gi_array = read_gi(input_file, N, tau);

  /* does it make sense to do it like this? */
  double omega_min = 0.0;
  double omega_max = 30000.0;
  double omega_step = 10.0;
  int num_steps = (int)((omega_max - omega_min)/omega_step);

  integral = calloc(N, sizeof(double));
  for (i = 0; i < N; i++) {
      integral[i] = calloc(num_steps, sizeof(double));
  }

  ex = calloc(tau, sizeof(double));
  for (unsigned int j = 0; j < num_steps; j++) {
    for (i = 0; i < N; i++) {
      w = omega_min + (j * omega_step);
      ex = exponent(w, wi[i], gamma[i], tau, gi_array[i]);
      /* integrate - tau is the number of steps in the integral */
      integral[i][j] = 2.0 * creal(trapezoid(ex, tau)) * pow(mu[i], 2.0);
    }
    /* this is probably a bad way of doing it actually -
     * can put all N integrals in one array, sum by element, then
     * multiply by omega outside the i loop here */
  }
}

double**
read_gi(char *input_file, int N, int tau)
{
  double **gi;
  FILE *fp, *gp;
  unsigned int j;
  char[200] line, file_line;
  gi_array = calloc(N, sizeof(double));
  for (unsigned int i = 0; i < N; i++) {
    gi_array[i] = calloc(tau, sizeof(double));
  }
  j = 0;

  /* read lines from input file. each one is a filename with a
   * line-broadening function saved in it. open each one and sum */
   fp = fopen(input_file, "r");
   if (fp == NULL) {
     fprintf(stdout, "Unable to open input file in read_gi.\n");
     exit(EXIT_FAILURE);
   } else {
      while (fgets(file_line, 199, fp) != NULL) {
        gp = fopen(file_line, "r")
        for (i = 0; i < tau; i++) {
          fgets(line);
          gi_array[j][i] = atof(line);
        }
        j++;
      }
   }
   return gi_array;
}

double complex*
exponent(double w, double w_i, double gamma_i,
	 unsigned int tau, double* gi)
{
  double complex *exponent;
  exponent = calloc(tau, sizeof(double complex));
  for (unsigned int i = 0; i < tau; i++) {
    exponent[i] = cexp(- I * i * (w - w_i) - gi[i] - (0.5 * gamma_i * i));
  }
  return exponent;
}

double complex
trapezoid(double complex *f, unsigned int n)
{
    double dx = 1./n;
    double complex sum;

    sum = 0.5 * dx * (f[0] + f[n - 1]);
    for (unsigned int i = 1; i < n - 1; i++) {
	sum += dx * f[i];
    }
    return sum;
}
