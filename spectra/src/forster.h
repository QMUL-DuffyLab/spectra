#include <fftw3.h>
#include <math.h>
#include "input.h"

#ifndef __FORSTER_H__
#define __FORSTER_H__

typedef struct forster_struct {
  unsigned ns;
  double w00;
  double lambda;
  double rate;
  double complex gi[];
} forster_struct;

double forster_rate(forster_struct A, forster_struct F, double J);

forster_struct* create_forster_struct(unsigned ns);
void assign_forster_struct(forster_struct* f,
                                      Parameters line_parameters);

static double
forster_overlap(forster_struct* A, forster_struct* F)
{
  fftw_complex *Atime, *Ftime, *Aw, *Fw;
  fftw_plan plan;

  if (A->ns != F->ns) {
    /* error */
  }

  Atime = (fftw_complex*) fftw_malloc(
           sizeof(fftw_complex) * A->ns);
  Ftime = (fftw_complex*) fftw_malloc(
           sizeof(fftw_complex) * F->ns);
  Aw    = (fftw_complex*) fftw_malloc(
           sizeof(fftw_complex) * A->ns);
  Fw    = (fftw_complex*) fftw_malloc(
           sizeof(fftw_complex) * F->ns);

  for (unsigned i = 0; i < A->ns; i++) {
    /* NB: rate or lifetime? also TOFS will be fine once
     * the includes are added. some problem with complex also */
    Atime[i] = At(A->w00,
               creal(A->gi[i]), cimag(A->gi[i]),
               (double)i * TOFS,
               1. / (1000000. * A->rate));
    Ftime[i] = Ft(F->w00,
               creal(F->gi[i]), cimag(F->gi[i]),
               F->lambda, (double)i * TOFS,
               1. / (1000000. * F->rate));
  }

  plan = fftw_plan_dft_1d(A->ns, 
  	 Atime, Aw, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(plan);
  plan = fftw_plan_dft_1d(A->ns, 
  	 Ftime, Fw, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(plan);

  double integral = 0.;
  for (unsigned i = 0; i < A->ns; i++) {
    /* NB: need to check this - fine to take real part? */
    integral += Aw[i][0] * Fw[i][0];
  }

  fftw_free(plan);
  fftw_free(Atime); fftw_free(Ftime);
  fftw_free(Aw); fftw_free(Fw);

  return integral;
}

#endif
