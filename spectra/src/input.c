#include <string.h>
#include <cmath>
#include "input.h"
#include <complex.h>

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
    p = (Input *)malloc(sizeof(Input) + N * sizeof(char));
    p->N = N;
    fgets(line, 199, fp);
    p->tau = atoi(line);
    fgets(line, 199, fp);
    p->n_chl = atoi(line);
    fgets(line, 199, fp);
    p->n_car = atoi(line);
    /* next block will be useful if/when I can set number of
     * carotenoids before the fortran code - as it is the way it's
     * hacked together in there means this will fail without lut 2 */
    /* if (p->n_chl + p->n_car != p->N) { */
    /*   fprintf(stdout, "Number of chlorophylls + number of" */
    /*       " carotenoids != total. Check the inputs!\n"); */
    /*   exit(EXIT_FAILURE); */
    /* } */
    fgets(line, 199, fp);
    p->T = atof(line);
    fgets(line, 199, fp);
    /* below is from stackoverflow - sets the newline to a null char */
    line[strcspn(line, "\n")] = 0;
    strcpy(p->eigvecs_file, line);
    fgets(line, 199, fp);
    line[strcspn(line, "\n")] = 0;
    strcpy(p->eigvals_file, line);
    /* fgets(line, 199, fp); */
    /* line[strcspn(line, "\n")] = 0; */
    /* strcpy(p->ei_file, line); */
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
    fgets(line, 199, fp);
    line[strcspn(line, "\n")] = 0;
    strcpy(p->jij_file, line);
    fgets(line, 199, fp);
    line[strcspn(line, "\n")] = 0;
    strcpy(p->com_file, line);
    for (i = 0; i < p->N; i++) {
      fgets(line, 199, fp);
      int newline = strcspn(line, "\n");
      line[newline] = 0;
      memcpy(p->gi_files[i], line, newline + 1);
    }
  }
  fclose(fp);
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
  arr = (double *)calloc(N, sizeof(double));

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

  fclose(fp);
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
  unsigned int i, token_size;
  token_size = 20;
  mu    = (double **)calloc(N, sizeof(double));
  token = (char *)malloc(token_size * sizeof(char));

  for (i = 0; i < N; i++) {
    mu[i] = (double *)calloc(3, sizeof(double));
  }

  fp = fopen(input_file, "r");
  if (fp == NULL) {
    fprintf(stdout, "Unable to open input file in read_mu."
        "Input file name was %s\n", input_file);
    exit(EXIT_FAILURE);
  } else {
    for (i = 0; i < N; i++) {
      fgets(token, token_size, fp);
      mu[i][0] = atof(token); 
      fgets(token, token_size, fp);
      mu[i][1] = atof(token); 
      fgets(token, token_size, fp); /* make sure we get to the newline! */
      mu[i][2] = atof(token); 
    }
  }

  fclose(fp);
  free(token);
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
  eig = (double **)calloc(N, sizeof(double*));
  for (unsigned int i = 0; i < N; i++) {
    eig[i] = (double *)calloc(N, sizeof(double));
  }
  j = 0;
  token = (char *)malloc(22 * sizeof(char));

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

  fclose(fp);
  free(token);
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
 * it generates gi[i][j], where \f$ i \in 0 \rightarrow N - 1 \f$, 
 * \f$ j \in 0 \rightarrow \tau - 1 \f$ .
 */
/* double _Complex** */
/* this is ugly - the hardcoded [200]. disgusting. i hate C++ */
fftw_complex**
read_gi(char input_files[N_MAX][200], 
    unsigned int N, unsigned int tau)
{
  fftw_complex** gi;
  FILE *fp;
  unsigned int i, j;
  double real, imag;
  gi = (fftw_complex**) fftw_malloc(
        sizeof(fftw_complex*) * N);
  /* gi = (double _Complex**)calloc(N, sizeof(double _Complex*)); */
  for (unsigned int i = 0; i < N; i++) {
    /* gi[i] = (double _Complex*)calloc(tau, sizeof(double _Complex)); */
    gi[i] = (fftw_complex*) fftw_malloc(
          sizeof(fftw_complex) * tau);
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
        /*
        * compiling with g++ to get the VERA C++ code to work
        * with everything else breaks this line because (I guess)
        * of differences between the C and C++ implementations
        * of complex numbers? even though they're both supposed to
        * just be two doubles? fuck sake. why is it like this.
        * **FORTRAN** is more user friendly than this and it was
        * pretty much invented on stone tablets!!!
        */
        /* std::complex<double> elem(real, imag); */
        /* gi[i][j] = reinterpret_cast<double(&)[2]>(elem); */
        gi[i][j][0] = real;
        gi[i][j][1] = imag;
        /* fprintf(stdout, "%10.6f\t%10.6f\n", */
        /*         creal(gi[i][j]), cimag(gi[i][j])); */
      }
    }
    fclose(fp);
  }
  return gi;
}

int
generate_filename(unsigned size, char *src,
                  const char *find, const char *replace)
{
  char *pch;
  pch = strstr(src, find);
  if(pch == NULL) {
    fprintf(stdout, "Cannot find string %s in string %s."
        " Source string will not be modified!\n", find, src);
    return 1;
  }
  if (strlen(src) - strlen(find) + strlen(replace) > size) {
    fprintf(stdout, "Length of modified string is greater than"
        " its size. String %s will not be modified!\n", src);
    return 2;
  }
  
  /* pch + strlen(find) points to the rest of the string after find.
   * we want to move this to pch + strlen(replace).
   * strlen(pch + strlen(find)) is the length of the string 
   * left after find, so this moves the end of the string
   * (strlen(replace) - strlen(find)) bytes along. The +1 is for
   * the terminating null character, which strlen ignores. */
  memmove(pch + strlen(replace), pch + strlen(find),
          strlen(pch + strlen(find)) + 1);
  /* then this slots replace in the gap */
  memcpy(pch, replace, strlen(replace));
  return 0;
}

