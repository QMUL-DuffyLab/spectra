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
    double integral;
    double *ex, *wi, *gamma;
    double **gi_array;

    wi = malloc(N, sizeof(double));
    gamma = malloc(N, sizeof(double));

    gi_array = malloc(N, sizeof(double));
    for (i = 0; i < N; i++) {
    	gi_array[i] = malloc(tau, sizeof(double));
    }
    gi_array = read_gi(input_file, N, tau);

    ex = malloc(tau, sizeof(double));
    for (i = 0; i < N; i++) {
	ex = exponent(w, wi[i], gamma[i], tau, gi_array[i]);
	/* integrate - tau is the number of steps in the integral */
	integral = trapezoid(ex, tau);
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
          gi_array[j][i] += atof(line);
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
