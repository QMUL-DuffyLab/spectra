#include "forster.h"

chromophore*
create_chromophore(unsigned ns)
{
  chromophore* c = malloc(sizeof(chromophore)
                 +   ns * sizeof(double complex));
  c->ns     = ns;
  c->w00    = NAN;
  c->lambda = NAN;
  c->rate   = NAN;
  for (unsigned i = 0; i < ns; i++) {
    c->gi[i] = NAN + I * NAN;
  }

  return c;
}

void
assign_chromophore(chromophore* c, Parameters line_parameters)
{
  /* c->w00 = line_parameters.energy; */
  /* c->lambda = line_parameters.lambda; */
  /* c->rate = 1./(1000. * line_parameters.gamma); */
  /* ns is assigned when the struct is created */
  for (unsigned i = 0; i < c->ns; i++) {
    /* c->gi[i] = line_parameters.gi[i]; */
  }
}

void
free_chromophore(chromophore* c)
{
  free(c);
  /* NB: I thought it was required to free the variable-length member
   * separately when you have one in a struct, but this free call
   * doesn't work. Is it wrong? Check valgrind */
  /* free(c->gi); */
}

static double
forster_overlap(chromophore* A, chromophore* F)
{
  fftw_complex *Atime, *Ftime, *Aw, *Fw;
  fftw_plan plan;

  if (A->ns != F->ns) {
    /* this should never happen but will break the FT so check */
    fprintf(stderr, "Error in forster_overlap: two g(t) lengths"
        "should be the same but are %u and %u for acceptor and"
        "donor respectively.\n", A->ns, F->ns);
    return NAN;
  }

  Atime = (fftw_complex*) fftw_malloc(
           sizeof(fftw_complex) * A->ns);
  Ftime = (fftw_complex*) fftw_malloc(
           sizeof(fftw_complex) * F->ns);
  Aw    = (fftw_complex*) fftw_malloc(
           sizeof(fftw_complex) * A->ns);
  Fw    = (fftw_complex*) fftw_malloc(
           sizeof(fftw_complex) * F->ns);

  double complex elem;
  for (unsigned i = 0; i < A->ns; i++) {
    /* NB: rate or lifetime? also TOFS will be fine once
     * the includes are added. some problem with complex also */
    elem = At(A->w00,
               creal(A->gi[i]), cimag(A->gi[i]),
               (double)i * TOFS,
               A->rate / (1000.));
    Atime[i][0] = creal(elem);
    Atime[i][1] = cimag(elem);

    elem = Ft(F->w00,
               creal(F->gi[i]), cimag(F->gi[i]),
               F->lambda, (double)i * TOFS,
               F->rate / (1000.));
    Ftime[i][0] = creal(elem);
    Ftime[i][1] = cimag(elem);
  }

  plan = fftw_plan_dft_1d(A->ns, 
  	 Atime, Aw, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(plan);
  /* calling fftw_plan_dft_1d again without destroying the plan
   * first seems to lead to a memory leak (so valgrind tells me) */
  fftw_destroy_plan(plan);

  plan = fftw_plan_dft_1d(F->ns, 
  	 Ftime, Fw, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(plan);

  /* find minima - need to subtract minimum values from the arrays */
  double A_min = 0., F_min = 0.;
  unsigned A_loc = 0, F_loc = 0;
  for (unsigned i = 0; i < A->ns; i++) {
    if (Aw[i][0] < Aw[A_loc][0]) {
      A_loc = i;
    }
    if (Fw[i][0] < Fw[F_loc][0]) {
      F_loc = i;
    }
  }
  A_min = Aw[A_loc][0]; F_min = Fw[F_loc][0];
  double A_sum = 0., F_sum = 0.;
  for (unsigned i = 0; i < A->ns; i++) {
    Aw[i][0] -= A_min;
    Fw[i][0] -= F_min;
    A_sum += Aw[i][0];
    F_sum += Fw[i][0];
  }
  fprintf(stdout, "A_min = %18.10f at w = %18.10f; Aw[A_loc][0] = %18.10f\n",
          A_min, index_to_frequency(A_loc, A->ns), Aw[A_loc][0]);
  fprintf(stdout, "F_min = %18.10f at w = %18.10f; Fw[F_loc][0] = %18.10f\n",
          F_min, index_to_frequency(F_loc, F->ns), Fw[F_loc][0]);

  /* now get norms so we can normalise the spectra */
  /* FFTW doesn't normalise but the Ns gives spurious results */
  /* double norm = 1. / (sqrt(A->ns) * fabs(sum)); */
  double A_norm = 1. / (fabs(A_sum));
  double F_norm = 1. / (fabs(F_sum));

  double integral = 0.;
  double dx = index_to_frequency(1, A->ns)
            - index_to_frequency(0, A->ns);

  /* NB: need to check this - fine to take real part? */
  integral = 0.5 * dx * (A_norm * F_norm) * ((Aw[0][0] * Fw[0][0])
           + (Aw[A->ns - 1][0] * Fw[F->ns - 1][0]));
  FILE *fp = fopen("out/aw_forster_test.dat", "w");
  FILE *gp = fopen("out/fw_forster_test.dat", "w");
  FILE *hp = fopen("out/forster_integral.dat", "w");
  for (unsigned i = 0; i < A->ns; i++) {
    fprintf(fp, "%18.10f\t%18.10f\n",
            index_to_frequency(i, A->ns), A_norm * Aw[i][0]);
    fprintf(gp, "%18.10f\t%18.10f\n",
            index_to_frequency(i, F->ns), F_norm * Fw[i][0]);
    fprintf(hp, "%18.10f\t%18.10f\t%18.10f\t%18.10f\n",
            index_to_frequency(i, F->ns), dx, A_norm * Aw[i][0], F_norm * Fw[i][0]);
    integral += dx * (A_norm * F_norm) * Aw[i][0] * Fw[i][0];
  }
  fclose(fp);
  fclose(gp);
  fclose(hp);
  
  fftw_free(Atime); fftw_free(Ftime);
  fftw_free(Aw); fftw_free(Fw);
  fftw_destroy_plan(plan);
  fftw_cleanup();

  fprintf(stdout, "raw integral = %12.8f\n", integral);

  return integral;
}

double
index_to_frequency(unsigned i, unsigned max_index)
{
  /* for some reason if i just return the expression the compiler
   * complains about control reaching end of non-void function? */
  double freq = ((double)i * 2. * PI) / ((double)max_index * TOFS);
  /* i think these give the exact same answer as the line above */
  /* double k = i * 2. * PI / (max_index); */
  /* double freq = (fmod(k, PI) - (PI * floor(k / PI))) / TOFS; */
  return freq;
}

double
forster_rate(chromophore* A, chromophore* F, double J)
{
  return 2. * pow(J, 2.) * forster_overlap(A, F);
}

double
forster_test()
{
  /* manually construct two chromophores corresponding
   * to CLA610 and LUT620, then check what rate we get */
  double step, real, imag;
  chromophore *CLA610 = create_chromophore((unsigned)2048);
  chromophore *LUT620 = create_chromophore((unsigned)2048);

  CLA610->w00    = 15100.0;
  CLA610->lambda = 37.0;
  CLA610->rate   = 1. / 4000.;

  LUT620->w00    = 14900.0;
  LUT620->lambda = 450.0;
  LUT620->rate   = 1. / 10.;

  double J610_620 = -5.5;

  FILE *fp = fopen("lineshape/out/CLA_gt.dat", "r");
  FILE *gp = fopen("lineshape/out/LUT_gt.dat", "r");
  if (fp == NULL) {
    fprintf(stdout, "Unable to open CLA g(t) file in forster test\n");
    exit(EXIT_FAILURE);
  } else {
    for (unsigned j = 0; j < CLA610->ns; j++) {
      /* NB: this fscanf call does not need the exact whitespace */
      int cl = fscanf(fp, "   %lf     %lf     %lf", &step, &real, &imag);
      if (cl != 3) {
        fprintf(stdout, "fscanf reading CLA_gt failed with error code %d;"
            " line number %d\n", cl, j);
        exit(EXIT_FAILURE);
      }
      CLA610->gi[j] = (real + I * imag);
      cl     = fscanf(gp, "   %lf     %lf     %lf", &step, &real, &imag);
      if (cl != 3) {
        fprintf(stdout, "fscanf reading LUT_gt failed with error code %d;"
            " line number %d\n", cl, j);
        exit(EXIT_FAILURE);
      }
      LUT620->gi[j] = (real + I * imag);
    }
  }
  fclose(fp);
  fclose(gp);

  double rate = forster_rate(CLA610, LUT620, J610_620);
  free_chromophore(CLA610);
  free_chromophore(LUT620);

  /* units??? */
  return (200. * PI * CVAC * 1E-12) * rate;

}
