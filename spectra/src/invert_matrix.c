#include "invert_matrix.h"

void
print_matrix(FILE *stream, char* name, unsigned n, double **matrix)
{
  unsigned i, j;
  if (name) {
    fprintf(stream, "Printing matrix %s:\n", name);
  }
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      fprintf(stream, "%8.4e ", matrix[i][j]);
    }
    fprintf(stream, "\n");
  }
}

void
print_vector(FILE *stream, char* name, unsigned n, double *vector)
{
  unsigned i;
  if (name) {
  fprintf(stream, "Printing vector %s:\n", name);
  }
  for (i = 0; i < n; i++) {
    fprintf(stream, "%8.4e\n", vector[i]);
  }
}

/** interchange between an aritrary 2d array and a contiguous one.
 *
 *  the arrays passed to LAPACK must be contiguous in memory.
 *  this means that a declaration of the form
 *  double **arr;
 *  arr = calloc(N, sizeof(double*)):
 *  for (i = 0; i < N; i++) {
 *    arr[i] = calloc(N, sizeof(double));
 *  }
 *  will NOT work.
 *  This function will either (with direction = 'F') assign
 *  a contiguous 1d array from the 2d one, or (with direction = 'B')
 *  do the reverse.
 */
void
interchange_2d_contiguous(char direction, unsigned n,
              double** array_2d, double* contiguous)
{
  unsigned i, j;
  if (direction == 'F') {
    for(i = 0; i < n; i++) {
      for(j = 0; j < n; j++) {
        contiguous[(i * n) + j] = array_2d[i][j];
      }
    }
  } else if (direction == 'B') {
    for(i = 0; i < n; i++) {
      for(j = 0; j < n; j++) {
        array_2d[i][j] = contiguous[(i * n) + j];
      }
    }
  } else {
    fprintf(stdout, "interchange_2d_contiguous function"
        " is broken. direction variable should be B or F but"
        " it's %c. ???\n", direction);
    exit(EXIT_FAILURE);
  }
}

/** very simple wrapper for LAPACKE double matrix inversion.
 *
 *  Does LU decomposition and then inversion
 *  of an n x n double matrix, out-of-place.
 *  NB: not necessary atm but I had problems with non-square
 *  matrices. Something to do with LDA and ipiv??? Maybe worth
 *  testing at some point later, can't be that hard
 */
int
invert_matrix_oop(unsigned n, double** in, double** out)
{
  lapack_int info, lda;
  lapack_int* ipiv;
  double *contiguous = (double *)calloc(n * n, sizeof(double));
  interchange_2d_contiguous('F', n, in, contiguous); 

  ipiv = (lapack_int*)calloc(n, sizeof(lapack_int));
  lda = n;

  info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, contiguous, lda, ipiv);
  if (info != 0) { /* if info != 0 dgetri won't work anyway */
    if (info < 0) {
      fprintf(stdout, "invert_matrix: dgetrf: argument %2d has an "
          "illegal value. Try again.\n", info);
      return info;
    } else if (info > 0) {
      fprintf(stdout, "invert_matrix: dgetrf: element (%2d, %2d) = 0."
          "Matrix won't be invertible.\n", info, info);
      return info;
    }
  }
  info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, contiguous, lda, ipiv);
  interchange_2d_contiguous('B', n, out, contiguous);
  free(contiguous);
  free(ipiv);
  return info;
}


/** very simple wrapper for inplace matrix inversion.
 *
 *  Does LU decomposition and then inversion
 *  of an n x n double matrix, inplace.
 */
int
invert_matrix_ip(unsigned n, double** inout)
{
  lapack_int info, lda;
  lapack_int* ipiv;
  double *contiguous = (double *)calloc(n * n, sizeof(double));

  ipiv = (lapack_int *)calloc(n, sizeof(lapack_int));
  lda = n;
  interchange_2d_contiguous('F', n, inout, contiguous);

  info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, contiguous, lda, ipiv);
  if (info != 0) {
    if (info < 0) {
      fprintf(stdout, "invert_matrix_inplace:"
          " dgetrf: argument %2d has an illegal value."
          " Try again.\n", info);
      return info;
    } else if (info > 0) {
      fprintf(stdout, "invert_matrix_inplace:"
          " dgetrf: element (%2d, %2d) = 0."
          " Matrix won't be invertible.\n", info, info);
      return info;
    }
  }
  info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, contiguous, lda, ipiv);
  interchange_2d_contiguous('B', n, inout, contiguous);
  free(ipiv);
  free(contiguous);
  return info;
}

/** wrapper for dgeev - right eigenvectors and eigenvalues.
 *
 * Note that vl shouldn't be necessary but it didn't seem to
 * work unless I allocated it and set it as a parameter, so I
 * just made it identical to vr for my own convenience. This
 * is pretty inefficient though and can probably be fixed.
 */
int
eig_oop(unsigned n, double** in, double** out, double* eigvals)
{
  lapack_int info, lda, ldvl, ldvr;
  unsigned i;
  double *wr, *wi, *vl;
  double *contiguous = (double *)calloc(n * n, sizeof(double));
  double *vr = (double *)calloc(n * n, sizeof(double));
  wr = (double *)calloc(n, sizeof(double));
  wi = (double *)calloc(n, sizeof(double));
  vl = (double *)calloc(n * n, sizeof(double));
  lda = n;
  ldvl = n; ldvr = n;
  interchange_2d_contiguous('F', n, in, contiguous);
  info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', n,
                       contiguous, lda, wr, wi, vl,
                       ldvl, vr, ldvr);

  for (i = 0; i < n; i++) {
    if(wi[i] > 1e-6) {
      fprintf(stdout, "eigenvalue %d is complex: = %12.8e + %+12.8e i\n",
          i, wr[i], wi[i]);
    }
    eigvals[i] = wr[i];
  }

  interchange_2d_contiguous('B', n, out, vr);
  free(wr); free(wi); free(vl); free(vr); free(contiguous);
  return info;
}

int
eig_ip(unsigned n, double** inout, double* eigvals)
{
  lapack_int info, lda, ldvl, ldvr;
  unsigned i;
  double *wr, *wi, *vl;
  double *contiguous = (double *)calloc(n * n, sizeof(double));
  double *vr = (double *)calloc(n * n, sizeof(double));
  wr = (double *)calloc(n, sizeof(double));
  wi = (double *)calloc(n, sizeof(double));
  vl = (double *)calloc(n * n, sizeof(double));
  lda = n;
  ldvl = n; ldvr = n;
  interchange_2d_contiguous('F', n, inout, contiguous);
  info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', n,
                       contiguous, lda, wr, wi, vl,
                       ldvl, vr, ldvr);

  for (i = 0; i < n; i++) {
    if(wi[i] > 1e-6) {
      fprintf(stdout, "eigenvalue %d is complex: = %12.8e + %+12.8e i\n",
          i, wr[i], wi[i]);
    }
    eigvals[i] = wr[i];
  }

  interchange_2d_contiguous('B', n, inout, vr);
  free(wr); free(wi); free(vl); free(vr); free(contiguous);
  return info;
}

void
matvec(unsigned n, double **mat, double *vec, double *res)
{
  unsigned i, j;
  double sum;
  for (i = 0; i < n; i++) {

    sum = 0.;

    for (j = 0; j < n; j++) {
      sum += mat[i][j] * vec[j];
      /* fprintf(stdout, "%2u %2u %10.6e %10.6e %10.6e\n", i, j, mat[i][j], vec[i], sum); */
    }
    res[i] = sum;

  }

}

void
matmul(unsigned n, double **a, double **b, double **ab)
{
  unsigned i, j, k;
  double sum;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {

      sum = 0.;

      for (k = 0; k < n; k++) {
        sum += a[i][k] * b[k][j];
      }

      ab[i][j] = sum;

    }
  }
}

