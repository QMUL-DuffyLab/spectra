#ifndef __HYBRID_VERA_H__
#define __HYBRID_VERA_H__
#include "vera.h"


std::vector<double>
k_i_xa_hybrid(VERA x, unsigned n_chl, unsigned n_car, unsigned tau,
               double **eig, double *eigvals, double **car_rates,
               double **Jij, double **normed_ai, double **normed_fi,
               pulse v_abs, double beta, double *car_decays);

double*
hybrid_boltz(unsigned n_chl, unsigned n_car, double beta,
             double *eigvals, VERA *x);
double**
hybrid_transfer(unsigned n_chl, unsigned n_car, VERA *x,
    double *gamma, double **Jij, std::vector<double> k_i_delta,
    double **redfield_rates, double *car_decays);

#endif
