#include <string.h>
#include "input.h"

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
        fgets(token, 17, fp);
        eig[i][0] = atof(token); 
        /* fprintf(stdout, "%d 0 %10.6e ", i, eig[i][0]); */
        for (j = 1; j < N - 1; j++) {
          fgets(token, 17, fp);
          eig[i][j] = atof(token); 
          /* fprintf(stdout, "%d %d %10.6e ", i, j, eig[i][j]); */
        }
        fgets(token, 18, fp); /* make sure we get to the newline! */
        eig[i][N - 1] = atof(token); 
        /* fprintf(stdout, "%d 13 %10.6e\n", i, eig[i][0]); */
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
  gi = calloc(N, sizeof(double complex*));
  for (unsigned int i = 0; i < N; i++) {
    gi[i] = calloc(tau, sizeof(double complex));
  }
  j = 0;

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
