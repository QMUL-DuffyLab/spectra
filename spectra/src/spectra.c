#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "../../lineshape/src/parameters.h"
#include "../../lineshape/src/functions.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <fftw3.h>

#define MAX_PIGMENT_NUMBER 200
#define EC 1.6022E-19
#define E0 8.854E-12
#define JPERINVCM 1.986E-23
#define KCOUL (1. / (4. * M_PI * E0 * 0.5))
/* #define CONV (2. * M_PI * CMS * 100. * 1E-15) */
#define TOFS (2. * M_PI * CMS * 100. * 1E-15)
#define TOCM1 ((1.295E3 * 8.988E9 * 0.5))

typedef struct {
  unsigned int N;
  char eigvecs_file[200], eigvals_file[200], mu_file[200],
  lambda_file[200], gamma_file[200];
  char *gi_files[];
} Input;

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
  eig = calloc(N, sizeof(double*));
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
        fgets(token, 20, fp);
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
  line = malloc(200 * sizeof(char*));
  token = malloc(22 * sizeof(char*));

  for (i = 0; i < N; i++) {
    fp = fopen(input_files[i], "r");
    if (fp == NULL) {
      fprintf(stdout, "Unable to open input file in read_gi."
          "Input file was %s\n", input_files[i]);
      exit(EXIT_FAILURE);
    } else {
      for (j = 0; j < tau; j++) {
        /* NB: there is definitely a bug in here!!! 
         * fgets will periodically lose its way a little bit
         * and start reading garbage characters from somewhere */
        fgets(token, 20, fp);
        fprintf(stdout, "%d %d %s", i, j, token);
        real = atof(token); 
        fgets(token, 22, fp); /* make sure we get to the newline! */
        fprintf(stdout, " %s", token);
        imag = atof(token); 
        gi[i][j] = (real + I * imag);
        fprintf(stdout, " %12.6f %12.6f\n", real, imag);
        /* gi[i][j] = (real + I * imag) * (1. / ((float)CMS * 100. * 1E-15 * 2. * M_PI)); */
        /* fprintf(stdout, "%d %d %18.10e %18.10e\n", i, j, gi[i][j]); */
      }
    }
  }
  return gi;
}

double complex*
exponent(double w_i, double gamma_i, double l1, double l2,
	 unsigned int tau, double complex* gi)
{
  double t;
  double complex *e;
  e = calloc(tau, sizeof(double complex));
  for (unsigned int i = 0; i < tau; i++) {
    t = (double) i * TOFS;
    /* see Kruger - should this be the line-broadening
     * function or just the lineshape function? :S */
    fprintf(stdout, "%d (%13.10e %+13.10ei)\n",
            i, creal(gi[i]), cimag(gi[i]));
    e[i] = cexp(-I * t * (w_i) - (creal(gi[i]) + I * cimag(gi[i]))
         - I * t * (l1 + l2)
         - (0.5 * (gamma_i) * t));
    /* e[i] = cexp(-I * t * (w_i) - (creal(gi[i]) + I * cimag(gi[i])) */
    /*      - I * t * (l1 + l2)); */
    /* e[i] = cexp(- (1. * gi[i] / (TOFS * 6.4)) */
    /*      - I * i * (l1 + l2) */
         /* - (0.5 * ((1E-6 / (gamma_i))) * i)); */
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
        kij[i][j] = TOCM1 * (pow(eig[i][k], 2.) * pow(eig[j][k], 2.) *
          p[k].cw(TOCM1 * wij[i][j], vptr));
        fprintf(stdout, "%d %d %d %18.10f ", i, j, k,
            p[k].cw(TOCM1 * wij[i][j], vptr));
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
  double musq, kd, freq;
  double complex *ex, **gi_array;
  double *eigvals, *gamma, *lambda, *integral,
         **wij, **kij, **mu, **eig;
  Parameters *line_params;
  fftw_complex *out, *in;
  fftw_plan plan;

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

  /* malloc 2d stuff */
  lineshape_files = malloc(p->N * sizeof(char*));
  mu = calloc(p->N, sizeof(double*));
  gi_array = calloc(p->N, sizeof(double complex*));
  eig = calloc(p->N, sizeof(double*));
  wij = calloc(p->N, sizeof(double*));
  kij = calloc(p->N, sizeof(double*));
  for (i = 0; i < p->N; i++) {
    lineshape_files[i] = malloc(200 * sizeof(char));
    gi_array[i] = calloc(tau, sizeof(double complex));
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
      line_params[i].cw = &cw_car;
    } else if (line_params[i].ligand == 1) {
      line_params[i].cw = &cw_chl;
    } else if (line_params[i].ligand == 2) {
      line_params[i].cw = &cw_odo;
    } else {
      fprintf(stdout, "Ligand number of params %d is %d."
          "No idea what's happened here.\n",
          i, line_params[i].ligand);
    }

    for (j = 0; j < p->N; j++) {
      wij[i][j] = ((eigvals[i] - lambda[i]) - (eigvals[j] - lambda[j]));
    }
  }
  int cl = fclose(fp);
  if (cl != 0) {
      fprintf(stdout, "Failed to close list of lineshape files %d.\n", cl);
      exit(EXIT_FAILURE);
  }

  kij = rate_calc(p->N, eig, wij, line_params);

  integral = calloc(tau, sizeof(double));

  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * tau);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * tau);
  plan = fftw_plan_dft_1d(tau, 
  	 in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

  ex = calloc(tau, sizeof(double complex));
  for (i = 0; i < p->N; i++) {
    musq = pow(mu[i][0], 2.) + pow(mu[i][2], 2.) + pow(mu[i][2], 2.);
    ex = exponent(eigvals[i], gamma[i],
        line_params[i].l1, line_params[i].l2,
        tau, gi_array[i]);

    for (unsigned int j = 0; j < tau; j++) {
      in[j] = At(0.0, creal(gi_array[i][j]),
            cimag(gi_array[i][j]), (double)j * TOFS,
            line_params[i].l1, line_params[i].l2, gamma[i]);
      /* fprintf(stdout, "%d (%18.10e %+18.10ei)" */
      /*     " (%18.10e %+18.10ei)\n", j, creal(ex[j]), cimag(ex[j]), */
      /*       creal(in[j]), cimag(in[j])); */
      /* in[j] = ex[j]; */
    } /* end j loop */

    fftw_execute(plan); 
    for (j = 0; j < tau; j++) {
      /* fprintf(stdout, "%d %18.10e\n", j, integral[j]); */
      /* fprintf(stdout, "%d (%18.10e %+18.10ei) (%18.10e %+18.10ei)\n", */
      /*     j, creal(in[j]), cimag(in[j]), creal(out[j]), cimag(out[j])); */
      integral[j] = musq * 2.0 * creal(out[j]);
    }
  }

  fp = fopen("out/aw_test.dat", "w");
  for (i = 0; i < tau; i++) {
    /* unpack the ordering used by FFTW */
    kd = i * 2. * M_PI / (tau);
    freq = fmod(kd, M_PI) - (M_PI * floor(kd / M_PI));
    fprintf(fp, "%18.10f %18.10f\n", freq / TOFS, 
        integral[i] * TOFS * (1./ sqrt(tau)) * 6.4);
  }
  cl = fclose(fp);
  if (cl != 0) {
      fprintf(stdout, "Failed to close A(w) output file %d.\n", cl);
      exit(EXIT_FAILURE);
  }

  free(line); free(lineshape_files); free(ex); free(integral);
  free(gi_array); free(eigvals); free(gamma); free(lambda); free(mu);
  free(eig); free(wij); free(kij); free(p); free(line_params);
  free(in); free(out);
  exit(EXIT_SUCCESS);
}
