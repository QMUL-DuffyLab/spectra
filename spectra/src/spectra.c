#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "../../lineshape/src/parameters.h"
#include "../../lineshape/src/functions.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
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
  lambda_file[200], gamma_file[200], aw_file[200], fw_file[200];
  char *gi_files[];
} Input;

typedef struct {
  /* this feels like it won't work lol */
  unsigned int N;
  double **kij;
  double *gamma;
} ode_params;

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
    fgets(line, 199, fp);
    line[strcspn(line, "\n")] = 0;
    strcpy(p->aw_file, line);
    fgets(line, 199, fp);
    line[strcspn(line, "\n")] = 0;
    strcpy(p->fw_file, line);
    for (i = 0; i < p->N; i++) {
      fgets(line, 199, fp);
      line[strcspn(line, "\n")] = 0;
      p->gi_files[i] = strndup(line, 200);
    }
  }
  return p;
}

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
  gi = calloc(N, sizeof(double complex*));
  for (unsigned int i = 0; i < N; i++) {
    gi[i] = calloc(tau, sizeof(double complex));
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
        int cl = fscanf(fp, "      %lf      %lf ", &real, &imag);
        if(cl != 2) {
          fprintf(stdout, "fscanf in read_gi failed with error code %d;"
              " line number %d\n", cl, j);
          exit(EXIT_FAILURE);
        }
        gi[i][j] = (real + I * imag);
      }
    }
  }
  return gi;
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

int
main(int argc, char** argv)
{
  FILE *fp;
  unsigned int tau, i, j;
  char *line, **lineshape_files;
  double musq, kd;
  double complex *ex, **gi_array;
  double *eigvals, *gamma, *lambda, *integral,
         **wij, **kij, **Jij, **mu, **eig;
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
  Jij = calloc(p->N, sizeof(double*));
  for (i = 0; i < p->N; i++) {
    lineshape_files[i] = malloc(200 * sizeof(char));
    gi_array[i] = calloc(tau, sizeof(double complex));
    mu[i] = calloc(3, sizeof(double));
    eig[i] = calloc(p->N, sizeof(double));
    wij[i] = calloc(p->N, sizeof(double));
    kij[i] = calloc(p->N, sizeof(double));
    Jij[i] = calloc(p->N, sizeof(double));

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
      /* line_params[i].cw = &cw_chl; */
      line_params[i].cw = &cw_odo;
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
  Jij = jacobian_calc(p->N, kij, gamma);
  fprintf(stdout, "Jacobian matrix:\n\n");
  for (i = 0; i < p->N; i++) {
    for (j = 0; j < p->N; j++) {
      fprintf(stdout, "%12.6e ", Jij[i][j]);
    }
    fprintf(stdout, "\n");
  }

  integral = calloc(tau, sizeof(double));

  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * tau);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * tau);
  plan = fftw_plan_dft_1d(tau, 
  	 in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

  ex = calloc(tau, sizeof(double complex));
  for (i = 0; i < p->N; i++) {
    musq = pow(mu[i][0], 2.) + pow(mu[i][2], 2.) + pow(mu[i][2], 2.);

    for (unsigned int j = 0; j < tau; j++) {
      in[j] = At(eigvals[i], creal(gi_array[i][j]), cimag(gi_array[i][j]),
                 (double)j * TOFS, line_params[i].l1, line_params[i].l2,
                 1. / gamma[i]);
    }

    fftw_execute(plan); 
    for (unsigned int j = 0; j < tau; j++) {
      integral[j] += creal(out[j]) * musq * 2.0;
    }

  }

  ode_params odep;
  odep.N = p->N;
  odep.kij = kij;
  odep.gamma = gamma;
  double *f = calloc(p->N, sizeof(double));
  double *y = calloc(p->N, sizeof(double));
  /* test boundary condition */
  for (i = 0; i < p->N; i++) {
    if (i == 0) {
      y[i] = 1.0;
    } else {
      y[i] = 0.0;
    }
  }
  double xtest = 0.0;
  void *params = &odep;

  int ode_success = odefunc(xtest, y, f, params);
  if (ode_success != GSL_SUCCESS) {
    fprintf(stdout, "ode_success failed!\n");
  }

  fp = fopen(p->aw_file, "w");
  for (i = 0; i < tau; i++) {
    /* unpack the ordering used by FFTW */
    kd = i * 2. * M_PI / (tau);
    fprintf(fp, "%18.10f %18.10f\n", kd / TOFS, 
        creal(integral[i]) * TOFS * (1./ sqrt(tau)) * 6.4);
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
