#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#ifndef PARAMETERS_H
#define PARAMETERS_H

#define CMS 299792458 /* speed of light _C_ in _M_etres per _S_econd */

/** Defines which spectral density ansatz to use.
 *
 * OBO    -> simple overdamped Brownian oscillator,
 * Renger -> ansatz from M\"{u}h and Renger derived from FLN data
 * BIG    -> from Novoderzhkin via Mancal - OBO with 48 
 *           high-frequency underdamped modes added on
 * CAR    -> carotenoid ansatz from Chris
 *
 * Note that currently I don't use the option CAR - at the moment,
 * the function `get_parameters` looks for the pigment code, and
 * assigns a variable `Parameters.ligand` based on what it finds.
 * If it's chlorophyll c we use the OBO, if it's a carotenoid we use
 * the carotenoid; only if it's a chlorophyll a or b do we use this
 * enum and the function `choose_ansatz`. I've left the CAR logic in
 * here just in case we ever want to add more spectral densities.
 */
typedef enum ansatz {
  OBO = 0,
  RENGER = 1,
  BIG = 2,
  CAR = 3,
} ansatz;

typedef struct {
    long unsigned int ns;
    double T;
    ansatz chl_ansatz;
} Protocol;

typedef struct {
    double s0, s1, s2, g0, g1, g2, l0, l1, l2,
           offset, w1, w2, ti, T, nu;
    double (* cw)(double, void *);
    double (* cn)(double, void *);
    double gsw[3][48];
    ansatz ans; /* 1 for chlorophyll, 0 for carotenoid */
    char aw_file[200], gt_file[200];
    char fw_file[200], lambda_file[200];
    char offset_file[200];
} Parameters;

Parameters get_parameters(char* filename, ansatz chl_ansatz);
Protocol get_protocol(char* filename);
Parameters* fortran_wrapper(int ligand);

#endif
