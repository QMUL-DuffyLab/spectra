#include <string.h>
#include "input.h"

/** Reads in the output from the fortran code, returns Input struct.
 *
 * This is not very portable and to some extent is duplicating effort
 * made in the lineshape code (lineshape/src/parameters.c) - there's
 * probably a way of fixing that / a library I could use instead but
 * whatever.
 *
 * Note that the number of gi_files isn't known till runtime because
 * neither is N, the number of pigments; luckily C99 allows us to leave
 * the last member of a struct with variable size and then malloc it
 * later, which is what I do here.
 *
 * fgets() reads in the newline at the end of the line as well; to
 * stop that newline from being passed to fopen() later I use the
 * line[strcspn()] = 0 trick to set the newline to a null char.
 *
 * strndup is the last remaining gnu99 extension I use in the whole
 * repo; I will eventually get around to fixing that but I'm lazy :)
 *
 * NB: check how to free those malloc'd struct members, I can't
 * remember if there's some fancy stuff you have to do or not.
 */
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
    p->tau = atoi(line);
    fgets(line, 199, fp);
    p->T = atof(line);
    fgets(line, 199, fp);
    /* below is from stackoverflow - sets the newline to a null char */
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
    fgets(line, 199, fp);
    line[strcspn(line, "\n")] = 0;
    strcpy(p->pop_file, line);
    for (i = 0; i < p->N; i++) {
      fgets(line, 199, fp);
      line[strcspn(line, "\n")] = 0;
      p->gi_files[i] = strndup(line, 200);
    }
  }
  return p;
}

/** Reads a file with one double per line and returns it as an array.
 *
 * Note that this and every other function in here returns a pointer
 * to a malloc'd object, so they all need to be freed somewhere!
 */
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

/** Reads transition dipole moments from file and returns array.
 *
 * This one (and the others below) are extremely hacky - the numbers
 * are hardcoded (19, 20, 22) based on the output format of the fortran
 * code. It disgusts me a little bit every time I remember writing it.
 * Not much going on other than that - reads 3 doubles a line, N lines,
 * returns a pointer to the resulting array.
 */
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

/** Reads in N x N eigenvector matrix and returns as an array.
 *
 * It occurs to me as I write these docs that I could do away with
 * separate functions just by taking two ints n and m and reading in
 * an n x m array from those. Would only require a check that both
 * are greater than 1.
 *
 * Other than that same deal, reads N x N floats, outputs N x N array.
 */
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
        for (j = 1; j < N - 1; j++) {
          fgets(token, 17, fp);
          eig[i][j] = atof(token); 
        }
        fgets(token, 18, fp); /* make sure we get to the newline! */
        eig[i][N - 1] = atof(token); 
      }
  }
  return eig;
}

/** Read in the set of line-broadening functions, return as 2d array.
 *
 * The g_i(t) functions are the mixed line-broadening functions output
 * by the fortran code, consisting of tau complex elements each.
 *
 * Honestly it's been a while since I wrote these functions and I can't
 * remember now why this one uses fscanf instead of the fgets etc. from
 * the other functions - maybe I could edit the others to use fscanf
 * instead and save myself some duplication of effort. Either way,
 * it generates gi[i][j], where i ∈ 0 -> N -1, j ∈ 0 -> tau - 1.
 */
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
        /* the space after the second %lf is required - the
         * fortran code adds a space after the imaginary part */
        int cl = fscanf(fp, "      %lf      %lf ", &real, &imag);
        if (cl != 2) {
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
