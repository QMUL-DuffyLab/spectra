#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include "../../lineshape/src/parameters.h"
#include "../../lineshape/src/functions.h"
/* #include <gsl/gsl_integration.h> */
/* #include <gsl/gsl_errno.h> */
#include <fftw3.h>

#define MAX_PIGMENT_NUMBER 200

typedef struct {
  unsigned int N;
  char eigvecs_file[200], eigvals_file[200], mu_file[200],
  lambda_file[200], gamma_file[200];
  char *gi_files[];
} Input;

typedef double (*cw_funptr)(double, void*);

double*
read(char *input_file, unsigned int N)
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
read_mu(char *input_file, unsigned int N)
{
  double **mu;
  FILE *fp;
  char *token;
  unsigned int i;
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
read_eigvecs(char *input_file, unsigned int N)
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
    fprintf(stdout, "Unable to open input file in read_eigvecs."
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
read_gi(char *input_files[], 
    unsigned int N, unsigned int tau)
{
  double complex **gi;
  FILE *fp;
  unsigned int i, j;
  double real, imag;
  char *line, *token;
  gi = calloc(N, sizeof(double*));
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
        gi[i][j] = (real + I * imag) * (1. / ((float)CMS * 100. * 1E-15 * 2. * M_PI));
        fprintf(stdout, "%d %d %18.10e %18.10e\n", i, j, gi[i][j]);
      }
    }
  }
  return gi;
}

double complex*
exponent(double w, double w_i, double gamma_i,
	 unsigned int tau, double complex* gi)
{
  double complex c1, c2, c3, *e;
  e = calloc(tau, sizeof(double complex));
  for (unsigned int i = 0; i < tau; i++) {
    /* see Kruger - should this be the line-broadening function
     * or just the lineshape function? the broadening one is 
     * divergent and it's the real part of this integral :S */
    c1 = cexp(- I * i * (w - w_i));
    c2 = - gi[i];
    c3 = - 0.5 * gamma_i * i;
    e[i] = cexp(- I * i * (w - w_i) - gi[i] - (0.5 * gamma_i * i));
    /* if (i < 100) { */
    /*   fprintf(stdout, "%d (%18.10e + %18.10ei) (%18.10e + %18.10ei) " */
    /*       "(%18.10e + %18.10ei)\n", i, creal(c1), cimag(c1), */
    /*        creal(c2), cimag(c2), creal(c3), cimag(c3)); */
    /* } */
  }
  return e;
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

Input*
read_input_file(char* filename)
{
  FILE *fp;
  Input *p;
  unsigned int i, N;
  char line[200];
  fp = fopen(filename, "r");
  if (fp == NULL) {
    fprintf(stdout, "Unable to open input file.\n");
    exit(EXIT_FAILURE);
  } else {
    fgets(line, 199, fp);
    N = atoi(line);
    /* now we know how much space we need for the list of g_is */
    p = malloc(sizeof(Input) + N * 200 * sizeof(char));
    p->N = N;
    fgets(line, 199, fp);
    /* from stackoverflow - sets the newline to a null char */
    line[strcspn(line, "\n")] = 0;
    strcpy(p->eigvecs_file, line);
    fgets(line, 199, fp);
    line[strcspn(line, "\n")] = 0;
    strcpy(p->eigvals_file, line);
    fgets(line, 199, fp);
    line[strcspn(line, "\n")] = 0;
    strcpy(p->mu_file, line);
    fgets(line, 199, fp);
    line[strcspn(line, "\n")] = 0;
    strcpy(p->lambda_file, line);
    fgets(line, 199, fp);
    line[strcspn(line, "\n")] = 0;
    strcpy(p->gamma_file, line);
    for (i = 0; i < p->N; i++) {
      fgets(line, 199, fp);
      line[strcspn(line, "\n")] = 0;
      p->gi_files[i] = strndup(line, 200);
    }
  }
  return p;
}

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
        kij[i][j] = pow(eig[i][k], 2.) * pow(eig[j][k], 2.) *
          p[k].cw(wij[i][j], vptr);
        fprintf(stdout, "%d %d %d %18.10f ", i, j, k,
            p[k].cw(wij[i][j], vptr));
      }
    fprintf(stdout, "\n");
    }
  }
  return kij;
}

int
main(int argc, char** argv)
{
  FILE *fp;
  unsigned int tau, i, j;
  char *line, **lineshape_files;
  double musq, w;
  double complex *ex, **gi_array;
  double *eigvals, *gamma, *lambda, *integral,
         **wij, **kij, **mu, **eig;
  Parameters *line_params;

  if (argc != 3) {
    fprintf(stdout, "Wrong number of arguments - should be 2: "
        "Input file and list of lineshape defs\n");
    exit(EXIT_FAILURE);
  }

  tau = 2000; /* again probably shouldn't hardcode this but oh well */

  Input *p = read_input_file(argv[1]);

  /* malloc 1d stuff, read them in */
  eigvals = calloc(p->N, sizeof(double));
  gamma = calloc(p->N, sizeof(double));
  lambda = calloc(p->N, sizeof(double));
  line_params = malloc(p->N * sizeof(Parameters));
  gamma = read(p->gamma_file, p->N);
  lambda = read(p->lambda_file, p->N);
  eigvals = read(p->eigvals_file, p->N);
  line = malloc(200 * sizeof(char));

  /* test */
  if (0) {
    gamma[0] = 0.1;
    gamma[1] = 4.0;
    gamma[2] = 4.0;
    gamma[3] = 4.0;
  }

  /* malloc 2d stuff */
  lineshape_files = malloc(p->N * sizeof(char*));
  mu = calloc(p->N, sizeof(double*));
  gi_array = calloc(p->N, sizeof(double*));
  eig = calloc(p->N, sizeof(double*));
  wij = calloc(p->N, sizeof(double*));
  kij = calloc(p->N, sizeof(double*));
  for (i = 0; i < p->N; i++) {
    lineshape_files[i] = malloc(200 * sizeof(char));
    gi_array[i] = calloc(tau, sizeof(double));
    mu[i] = calloc(3, sizeof(double));
    eig[i] = calloc(p->N, sizeof(double));
    wij[i] = calloc(p->N, sizeof(double));
    kij[i] = calloc(p->N, sizeof(double));

  }

  gi_array = read_gi(p->gi_files, p->N, tau);
  eig = read_eigvecs(p->eigvecs_file, p->N);
  mu = read_mu(p->mu_file, p->N);

  fp = fopen(argv[2], "r"); /* read in list of lineshape files here */
  for (i = 0; i < p->N; i++) {
    /* now load in the parameters for each ligand */
    fgets(line, 200, fp);
    line[strcspn(line, "\n")] = 0;
    strcpy(lineshape_files[i], line);
    line_params[i] = get_parameters(lineshape_files[i]);
    if (line_params[i].ligand == 0) {
      line_params[i].cw = &cw_chl;
    } else if (line_params[i].ligand == 1) {
      line_params[i].cw = &cw_car;
    } else if (line_params[i].ligand == 2) {
      line_params[i].cw = &cw_odo;
    } else {
      fprintf(stdout, "Ligand number of params %d is %d."
          "No idea what's happened here.\n",
          i, line_params[i].ligand);
    }

    for (j = 0; j < p->N; j++) {
      wij[i][j] = (eigvals[i] - lambda[i]) - (eigvals[j] - lambda[j]);
    }
  }
  int cl = fclose(fp);
  if (cl != 0) {
      fprintf(stdout, "Failed to close list of lineshape files %d.\n", cl);
      exit(EXIT_FAILURE);
  }


  kij = rate_calc(p->N, eig, wij, line_params);

  /* does it make sense to do it like this? */
  double omega_min = -10000.0;
  double omega_max = 10000.0;
  double omega_step = 10.0;
  unsigned int num_steps = (int)((omega_max - omega_min)/omega_step);

  integral = calloc(num_steps, sizeof(double));

  ex = calloc(tau, sizeof(double));
  for (i = 0; i < p->N; i++) {
    musq = pow(mu[i][0], 2.) + pow(mu[i][2], 2.) + pow(mu[i][2], 2.);
    for (unsigned int j = 0; j < num_steps; j++) {
      w = omega_min + (j * omega_step);
      ex = exponent(w, eigvals[i], gamma[i], tau, gi_array[i]);
      /* integrate - tau is the number of steps in the integral */
      integral[j] += musq * 2.0 * creal(trapezoid(ex, tau));
      /* fprintf(stdout, "%d %18.10f %18.10f\n", i, w, integral[j]); */
    }
  }
  for (j = 0; j < tau; j++) {
    integral[j] *= (omega_min + j * omega_step);
  }

  fp = fopen("out/aw_test.dat", "w");
  for (i = 0; i < num_steps; i++) {
    fprintf(fp, "%16.8e %16.8e\n", omega_min + (i * omega_step), 
        integral[i]);
  }
  cl = fclose(fp);
  if (cl != 0) {
      fprintf(stdout, "Failed to close A(w) output file %d.\n", cl);
      exit(EXIT_FAILURE);
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

  free(line); free(lineshape_files); free(ex); free(integral);
  free(gi_array); free(eigvals); free(gamma); free(lambda); free(mu);
  free(eig); free(wij); free(kij); free(p); free(line_params);
  exit(EXIT_SUCCESS);
}
