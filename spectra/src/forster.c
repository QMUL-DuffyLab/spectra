#include "forster.h"

double
forster_rate(forster_struct* A, forster_struct* F, double J)
{
  return 2 * M_PI * pow(J, 2.) * forster_overlap(A, F);
}

forster_struct*
create_forster_struct(unsigned ns)
{
  forster_struct* f = malloc(sizeof(forster_struct)
                      + ns * sizeof(double complex));
  f->ns     = ns;
  f->w00    = NAN;
  f->lambda = NAN;
  f->rate   = NAN;
  for (unsigned i = 0; i < ns; i++) {
    f->gi[i] = NAN + I * NAN;
  }

  return f;
}

void
assign_forster_struct(forster_struct* f, Parameters line_parameters)
{
  f->w00 = line_parameters.energy;
  f->lambda = line_parameters.lambda;
  f->rate = 1./(1000. * line_parameters.gamma);
  /* ns is assigned when the struct is created */
  for (unsigned i = 0; i < f->ns; i++) {
    f->gi[i] = line_parameters.gi[i];
  }
}

void
free_forster_struct(forster_struct* f)
{
  free(f->gi);
  free(f);
}
