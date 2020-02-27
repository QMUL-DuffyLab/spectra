#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#include <string.h>
/* #include <gsl/gsl_integration.h> */
/* #include <gsl/gsl_errno.h> */
#include <fftw3.h>

#define MAX_PIGMENT_NUMBER 200

typedef struct {
  int N;
  char eigvecs_file[200], eigvals_file[200], mu_file[200],
  lambda_file[200], gamma_file[200];
  char gi_files[MAX_PIGMENT_NUMBER][200];
} Parameters;

double*
read(char *input_file, int N)
{
  double *arr;
  FILE *fp;
  char line[200];
  unsigned int i;
  arr = calloc(N, sizeof(double));

  fp = fopen(input_file, "r");
  if (fp == NULL) {
    fprintf(stdout, "Unable to open input file in read function. "
        "File name was %s\n", input_file);
    exit(EXIT_FAILURE);
  } else {
    for (i = 0; i < N; i++) {
      fgets(line, 199, fp);
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
  char *token;
  unsigned int i, j;
  mu = calloc(N, sizeof(double));
  token = malloc(22 * sizeof(char));

  for (i = 0; i < N; i++) {
    mu[i] = calloc(3, sizeof(double));
  }

  i = 0;
  fp = fopen(input_file, "r");
  if (fp == NULL) {
    fprintf(stdout, "Unable to open input file in read_mu."
        "Input file name was %s\n", input_file);
    exit(EXIT_FAILURE);
  } else {
    for (i = 0; i < N; i++) {
      fgets(token, 19, fp);
      mu[i][0] = atof(token); 
      fgets(token, 20, fp);
      mu[i][1] = atof(token); 
      fgets(token, 22, fp); /* make sure we get to the newline! */
      mu[i][2] = atof(token); 
    }
  }
  return mu;
}

double**
read_eigvecs(char *input_file, int N)
{
  double **eig;
  FILE *fp;
  unsigned int i, j;
  char *token;
  eig = calloc(N, sizeof(double));
  for (unsigned int i = 0; i < N; i++) {
    eig[i] = calloc(N, sizeof(double));
  }
  j = 0;
  token = malloc(22 * sizeof(char));

  fp = fopen(input_file, "r");
  if (fp == NULL) {
    fprintf(stdout, "Unable to open input file in read_gi."
        "Input file was %s\n", input_file);
    exit(EXIT_FAILURE);
  } else {
    for (i = 0; i < N; i++) {
        fgets(token, 19, fp);
        eig[i][0] = atof(token); 
        for (j = 1; j < N - 1; j++) {
          fgets(token, 20, fp);
          eig[i][j] = atof(token); 
        }
        fgets(token, 22, fp); /* make sure we get to the newline! */
        eig[i][N - 1] = atof(token); 
      }
  }
  return eig;
}

double complex**
read_gi(char input_files[MAX_PIGMENT_NUMBER][200], int N, int tau)
{
  double complex **gi;
  FILE *fp;
  unsigned int i, j;
  double real, imag;
  char *line, *token;
  gi = calloc(N, sizeof(double));
  for (unsigned int i = 0; i < N; i++) {
    gi[i] = calloc(tau, sizeof(double));
  }
  j = 0;
  line = malloc(200 * sizeof(char));
  token = malloc(22 * sizeof(char));

  for (i = 0; i < N; i++) {
    fp = fopen(input_files[i], "r");
    if (fp == NULL) {
      fprintf(stdout, "Unable to open input file in read_gi."
          "Input file was %s\n", input_files[i]);
      exit(EXIT_FAILURE);
    } else {
      for (j = 0; j < tau; j++) {
        fgets(token, 19, fp);
        real = atof(token); 
        fgets(token, 22, fp); /* make sure we get to the newline! */
        imag = atof(token); 
        gi[i][j] = real + I * imag;
      }
    }
  }
  return gi;
}

double complex*
exponent(double w, double w_i, double gamma_i,
	 unsigned int tau, double complex* gi)
{
  double complex *exponent;
  exponent = calloc(tau, sizeof(double complex));
  for (unsigned int i = 0; i < tau; i++) {
    /* see Kruger - should this be the line-broadening function
     * or just the lineshape function? the broadening one is 
     * divergent and it's the real part of this integral :S */
    exponent[i] = cexp(- I * i * (w - w_i) - gi[i] - (0.5 * gamma_i * i));
    /* if (i < 100) { */
      /* fprintf(stdout, "%18.10f %18.10f\n", creal(exponent[i]), cimag(exponent[i])); */
    /* } */
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

Parameters
read_input_file(char* filename)
{
  FILE *fp;
  Parameters p;
  unsigned int i;
  char line[200];
  fp = fopen(filename, "r");
  if (fp == NULL) {
    fprintf(stdout, "Unable to open input file.\n");
    exit(EXIT_FAILURE);
  } else {
    fgets(line, 199, fp);
    p.N = atoi(line);
    fgets(line, 199, fp);
    /* from stackoverflow - sets the newline to a null char */
    line[strcspn(line, "\n")] = 0;
    strcpy(p.eigvecs_file, line);
    fgets(line, 199, fp);
    line[strcspn(line, "\n")] = 0;
    strcpy(p.eigvals_file, line);
    fgets(line, 199, fp);
    line[strcspn(line, "\n")] = 0;
    strcpy(p.mu_file, line);
    fgets(line, 199, fp);
    line[strcspn(line, "\n")] = 0;
    strcpy(p.lambda_file, line);
    fgets(line, 199, fp);
    line[strcspn(line, "\n")] = 0;
    strcpy(p.gamma_file, line);
    for (i = 0; i < p.N; i++) {
      fgets(line, 199, fp);
      line[strcspn(line, "\n")] = 0;
      strcpy(p.gi_files[i], line);
    }
  }
  return p;
}

int
main(int argc, char** argv)
{
  unsigned int N, tau, i, j;
  double musq, w;
  double complex *ex, *integral, **gi_array;
  double *wi, *gamma, *lambda, **mu, **eig;

  tau = 2000; /* again probably shouldn't hardcode this but oh well */

  Parameters p = read_input_file(argv[1]);

  /* need to read these in from files as well */
  wi = calloc(p.N, sizeof(double));
  gamma = calloc(p.N, sizeof(double));
  lambda = calloc(p.N, sizeof(double));

  mu = calloc(p.N, sizeof(double));
  gi_array = calloc(p.N, sizeof(double));
  for (i = 0; i < p.N; i++) {
    gi_array[i] = calloc(tau, sizeof(double));
    mu[i] = calloc(3, sizeof(double));
  }
  gi_array = read_gi(p.gi_files, p.N, tau);
  eig = read_eigvecs(p.eigvecs_file, p.N);
  mu = read_mu(p.mu_file, p.N);
  gamma = read(p.gamma_file, p.N);
  lambda = read(p.lambda_file, p.N);
  wi = read(p.eigvals_file, p.N);

  for (i = 0; i < p.N; i++) {
    for (j = 0; j < p.N; j++) {
      fprintf(stdout, "%18.10f ", eig[i][j]); 
    }
    fprintf(stdout, "\n"); 
  }

  /* does it make sense to do it like this? */
  double omega_min = 10000.0;
  double omega_max = 30000.0;
  double omega_step = 10.0;
  int num_steps = (int)((omega_max - omega_min)/omega_step);

  integral = calloc(num_steps, sizeof(double));

  ex = calloc(tau, sizeof(double));
  for (i = 0; i < p.N; i++) {
    musq = pow(mu[i][0], 2.) + pow(mu[i][2], 2.) + pow(mu[i][2], 2.);
    for (unsigned int j = 0; j < num_steps; j++) {
      w = omega_min + (j * omega_step);
      ex = exponent(w, wi[i], gamma[i], tau, gi_array[i]);
      /* integrate - tau is the number of steps in the integral */
      integral[j] += w * musq * 2.0 * creal(trapezoid(ex, tau));
    }
  }

  FILE *fp = fopen("out/aw_test.dat", "w");
  for (i = 0; i < num_steps; i++) {
    fprintf(fp, "%16.8e (%16.8e + %16.8ei)\n", omega_min + (i * omega_step), 
        creal(integral[i]), cimag(integral[i]));
  }

  /* should this maybe be an FFT instead of a trapezoid thing????
   * need to sit down and write it out maybe. have a think */

  /* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * pr.ns); */

  /* for (unsigned long i = 0; i < tau; i++) { */
  /*     in[i][0] = creal(Atv[i]); */
  /*     in[i][1] = cimag(Atv[i]); */
  /* } */

  /* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * tau); */
  /* plan = fftw_plan_dft_1d(tau, */ 
  /* 	 in, out, FFTW_BACKWARD, FFTW_ESTIMATE); */
  /* fftw_execute(plan); */ 
  /* unpack the ordering used by FFTW */
  /* double k = i * 2. * M_PI / (tau); */
  /* double freq = fmod(k, M_PI) - (M_PI * floor(k / M_PI)); */

  exit(EXIT_SUCCESS);
}

