#include "invert_matrix.h"

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
  unsigned i, j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      out[i][j] = in[i][j];
    }
  }

  ipiv = calloc(n, sizeof(lapack_int));
  lda = n;

  info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, *out, lda, ipiv);
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
  info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, *out, lda, ipiv);
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

  ipiv = calloc(n, sizeof(lapack_int));
  lda = n;

  info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, *inout, lda, ipiv);
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
  info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, *inout, lda, ipiv);
  return info;
}
