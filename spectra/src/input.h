#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "../../lineshape/src/parameters.h"
#include "../../lineshape/src/functions.h"

#ifndef __INPUT_H__
#define __INPUT_H__

#define MAX_PIGMENT_NUMBER 200
#define EC 1.6022E-19
#define E0 8.854E-12
#define JPERINVCM 1.986E-23
#define KCOUL (1. / (4. * M_PI * E0 * 0.5))
/* #define CONV (2. * M_PI * CMS * 100. * 1E-15) */
#define TOFS (2. * M_PI * CMS * 100. * 1E-15)
#define TOCM1 ((1.295E3 * 8.988E9 * 0.5))

/** Input struct for important parameters and filenames where
 * the exciton-basis data is to be read in from.
 * 
 * This struct mainly holds a list of filenames which have been
 * output by the fortran code, listing where all the data is for
 * the eigenstates of our Hamiltonian (the exciton basis data).
 * Also holds a few parameters that aren't specific
 *  to any particular pigment or exciton.
 */
typedef struct {
  unsigned int N; /**< Number of pigments/excitons */
  double T; /**< Temperature */
  unsigned int tau; /**< number of steps in g(t) arrays (\equiv fs) */
  char eigvecs_file[200], eigvals_file[200], mu_file[200],
  lambda_file[200], gamma_file[200], aw_file[200], fw_file[200],
  pop_file[200];
  char *gi_files[];
} Input;

Input* read_input_file(char* filename);
double* read(char *input_file, unsigned int N);
double** read_mu(char *input_file, unsigned int N);
double** read_eigvecs(char *input_file, unsigned int N);

double complex** read_gi(char *input_files[], 
		 unsigned int N, unsigned int tau);

#endif
