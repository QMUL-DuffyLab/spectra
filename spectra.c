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
  double musq, w;
  double *ex, *wi, *gamma, *lambda, *integral;
  double **mu, **gi_array;

  /* need to read these in from files as well */
  wi = calloc(N, sizeof(double));
  gamma = calloc(N, sizeof(double));
  lambda = calloc(N, sizeof(double));

  mu = calloc(N, sizeof(double));
  gi_array = calloc(N, sizeof(double));
  for (i = 0; i < N; i++) {
    gi_array[i] = calloc(tau, sizeof(double));
    mu[i] = calloc(3, sizeof(double));
  }
  gi_array = read_gi(input_file, N, tau);

  /* does it make sense to do it like this? */
  double omega_min = 10000.0;
  double omega_max = 30000.0;
  double omega_step = 10.0;
  int num_steps = (int)((omega_max - omega_min)/omega_step);

  integral = calloc(num_steps, sizeof(double));

  ex = calloc(tau, sizeof(double));
  for (i = 0; i < N; i++) {
    musq = pow(mu[i][0], 2.) + pow(mu[i][2], 2.) + pow(mu[i][2], 2.);
    for (unsigned int j = 0; j < num_steps; j++) {
      w = omega_min + (j * omega_step);
      ex = exponent(w, wi[i], gamma[i], tau, gi_array[i]);
      /* integrate - tau is the number of steps in the integral */
      integral[j] += w * musq * 2.0 * creal(trapezoid(ex, tau));
    }
  }

  FILE *fp = fopen(aw_file, "w");
  for (i = 0; i < num_steps; i++) {
    fprintf(fp, "%16.8e %16.8e\n", omega_min + (i * omega_step), integral[i]);
  }

  /* should this maybe be an FFT instead of a trapezoid thing????
   * need to sit down and write it out maybe. have a think */

  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * pr.ns);

  for (unsigned long i = 0; i < tau; i++) {
      in[i][0] = creal(Atv[i]);
      in[i][1] = cimag(Atv[i]);
  }

  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * tau);
  plan = fftw_plan_dft_1d(tau, 
  	 in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(plan); 
  /* unpack the ordering used by FFTW */
  double k = i * 2. * M_PI / (tau);
  double freq = fmod(k, M_PI) - (M_PI * floor(k / M_PI));

  exit(EXIT_SUCCESS);
}

double**
read(char *input_file, int N)
{
  double *arr;
  FILE *fp;
  char[200] line;
  unsigned int i;
  arr = calloc(N, sizeof(double));

  fp = fopen(input_file, "r");
  if (fp == NULL) {
    fprintf(stdout, "Unable to open input file in read function. ",
        "File name was %s\n", input_file);
    exit(EXIT_FAILURE);
  } else {
    for (i = 0; i < N; i++) {
      fgets(line, 199, gp);
      arr[i] = atof(line);
    }
  }
  return arr;
}

double**
read_mu(char *input_file, int N)
{
  double **mu;
  FILE *fp;
  char[200] line;
  unsigned int i;
  mu = calloc(N, sizeof(double));

  for (i = 0; i < N; i++) {
    mu[i] = calloc(3, sizeof(double));
  }

  i = 0;
  fp = fopen(input_file, "r");
  if (fp == NULL) {
    fprintf(stdout, "Unable to open input file in read_gi.\n");
    exit(EXIT_FAILURE);
  } else {
    for (i = 0; i < N; i++) {
      fgets(line, 199, gp);
      /* not right yet */
      mu[i][0] = atof(line);
      mu[i][1] = atof(line);
      mu[i][2] = atof(line);
    }
  }
  return mu;
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

  /* need to check N lines are read in to file_line, i.e. that 
   * the number of g_i's we're reading in is correct.
   * read lines from input file. each one is a filename with a
   * line-broadening function saved in it. open each one and sum */
   fp = fopen(input_file, "r");
   if (fp == NULL) {
     fprintf(stdout, "Unable to open input file in read_gi.\n");
     exit(EXIT_FAILURE);
   } else {
      while (fgets(file_line, 199, fp) != NULL) {
        gp = fopen(file_line, "r")
        for (i = 0; i < tau; i++) {
          fgets(line, 199, gp);
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
