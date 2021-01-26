#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <stdlib.h>
#include "helper.h"
#include "LSODA.h"
#include "forster.h"
#include "steady_state.h"

#ifndef __VERA_H__
#define __VERA_H__

class VERA {
  public:
    VERA(size_t elec, size_t norm, size_t vib,
                double beta,
                double* wi, size_t n_wi,
                double* w_mode, size_t n_mode,
                double** l_ic, size_t side_l,
                double* l_ivr, size_t n_l,
                double** g_ic, size_t side_g,
                double* g_ivr, size_t n_g,
                double*** di, size_t n_mode_di, size_t e_l, size_t e_u
                );
    VERA(size_t elec, size_t norm, size_t vib,
                double beta,
                std::vector<double> w_elec,
                std::vector<double> w_mode,
                std::vector<double> l_ic,
                std::vector<double> l_ivr,
                std::vector<double> g_ic,
                std::vector<double> g_ivr,
                std::vector<double> di
                );
    /* this needs changing - should the set functions be public? */
    void set_extents();
    void set_w_elec(double* wi, size_t n_wi);
    void set_w_elec(std::vector<double> wi);
    void set_w_normal(double* w_mode, size_t n_mode);
    void set_w_normal(std::vector<double> w_mode);
    void set_l_ic_ij(double** l_ic, size_t size);
    void set_l_ic_ij(std::vector<double> l_ic);
    void set_l_ivr_i(double* l_ivr, size_t n);
    void set_l_ivr_i(std::vector<double> l_ivr);
    void set_g_ic_ij(double** g_ic, size_t size);
    void set_g_ic_ij(std::vector<double> g_ic);
    void set_g_ivr_i(double* g_ivr, size_t n);
    void set_g_ivr_i(std::vector<double> g_ivr);
    void set_disp(double*** di, size_t n_mode, size_t e_l, size_t e_u);
    void set_disp(std::vector<double> di);
    void fc_calc();
    void set_k_ivr(double beta);
    void set_k_ic(double beta);
    std::vector<double> dndt(double *population,
                             double t, pulse pump);
    double get_w_elec(size_t i);
    double get_w_normal(size_t i);
    std::vector<size_t> get_extents();
    std::vector<size_t> get_pop_extents();
    std::vector<double> get_w_normal();
    double get_w_vib(size_t i);
    double get_l_ic_ij(size_t i, size_t j);
    double get_l_ivr_i(size_t i);
    double get_g_ic_ij(size_t i, size_t j);
    double get_g_ivr_i(size_t i);
    /* The indexing is disp[mode][elec_lower][elec_upper] */
    double get_disp(std::vector<size_t> subscripts);
    double get_fc(std::vector<size_t> subscripts);
    void get_k_ivr(double beta, size_t elec, size_t norm);
    size_t n_elec;
    size_t n_normal;
    size_t n_vib;
    double beta;

  private:
    std::vector<size_t> extents;
    std::vector<size_t> ivr_extents;
    std::vector<size_t> ic_extents;
    std::vector<size_t> disp_extents;
    std::vector<size_t> fc_extents;
    std::vector<size_t> population_extents;
    std::vector<double> w_elec;
    std::vector<double> w_normal;
    /* std::vector<double> w_vib; */
    /** \lambda_ic_ij
     *
     * Reorganisation energy for interconversion from electronic state j -> i.
     * Doesn't need to include transitions from S_n.
     */
    std::vector<double> l_ic_ij;
    /** \lambda_ivr_i
     *
     * Reorganisation energy for Intramolecular Vibrational Redistribution
     * on a single electronic state. This does not need to include Sn
     */
    std::vector<double> l_ivr_i;
    /** G^{IC}_{ij}, G^{IVR}_{i}
     *
     * Damping times for IC and IVR (in fs) for IC from electronic state j -> i. 
     * Convert to wavenumbers later.
     * They are assumed to be identical for each normal
     * mode, and again we don't need to include S_n.
    */
    std::vector<double> g_ic_ij;
    std::vector<double> g_ivr_i;
    std::vector<double> disp;
    std::vector<double> fc;
    std::vector<double> k_ivr; /* rates for IVR */
    std::vector<double> k_ic;  /* rates for IC */
};

/*
 * LSODA struct and function
 */

typedef struct vera_lsoda_data {
  VERA *chromo;
  pulse *pump;
} vera_lsoda_data;

double C_OBO(double w, double l, double g);
size_t sub2ind(std::vector<size_t> subscripts,
       std::vector<size_t> extents);
std::vector<size_t> ind2sub(size_t index, std::vector<size_t> extents);
VERA create_VERA_from_file(char *filename);
void func(double t, double *y, double *ydot, void *data);
double thermal_osc(std::vector<size_t> a,
                   std::vector<double> w_alpha, double beta);
double 
k_inter(std::vector<size_t> e_ij,
    std::vector<size_t> a, std::vector<size_t> b,
    std::vector<double> w_e, std::vector<double> w_alpha,
    double beta, double lambda, double gamma);

double k_calc(double w, double beta, double lambda, double gamma);

#endif
