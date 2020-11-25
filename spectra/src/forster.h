#include <fftw3.h>
#include <math.h>
#include "input.h"

#ifndef __FORSTER_H__
#define __FORSTER_H__

/** Struct to hold all the info we need to calculate Redfield/Forster rates
 *
 * it's called chromophore because it can equally describe a
 * pigment or a superposition of them (exciton)
 */
typedef struct chromophore {
  unsigned ns;
  double w00;
  double lambda;
  double rate; /* rate in inverse picoseconds!! */
  double complex gi[];
} chromophore;

chromophore* create_chromophore(unsigned ns);
void assign_chromophore(chromophore* c, Parameters line_parameters);
void free_chromophore(chromophore* c);

double forster_rate(chromophore* A, chromophore* F, double J);
double index_to_frequency(unsigned i, unsigned max_index);

double forster_test();

#endif
