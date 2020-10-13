#include "invert_matrix.h"

void
print_matrix(char* name, unsigned n, double **matrix)
{
  unsigned i, j;
  fprintf(stdout, "Printing matrix %s:\n", name);
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      fprintf(stdout, "%10.6f ", matrix[i][j]);
    }
    fprintf(stdout, "\n");
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
        " it's %s. ???\n", direction);
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
  double *contiguous = calloc(n * n, sizeof(double));
  interchange_2d_contiguous('F', n, in, contiguous); 

  ipiv = calloc(n, sizeof(lapack_int));
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
  double *contiguous = calloc(n * n, sizeof(double));

  ipiv = calloc(n, sizeof(lapack_int));
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
  return info;
}
