#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <stdlib.h>
#include "helper.h"
#include "LSODA.h"
/* #include "forster.h" */
/* #include "steady_state.h" */

#ifndef __VER_H__
#define __VER_H__

/* this stuff is from elsewhere - i want this to compile on its own lol
 * worry about this stuff later */

class Chromophore;

#define CVAC 299792458.0
#define CM_PER_PS (200. * M_PI * CVAC * 1E-12)
typedef enum pulse_type {
  FLAT = 0,
  LORENTZIAN = 1,
  GAUSSIAN = 2,
  DELTA = 3,
} pulse_type;

typedef struct pulse {
  pulse_type type;
  size_t target_state;
  double amplitude;
  double centre;
  double width;
  double t_peak;
  double duration;
} pulse;

double
intensity(double w, double t, pulse p)
{
  double sigma, f;
  double gamma = 1./(2. * p.duration)
               * pow(1./cosh((t - p.t_peak)/p.duration), 2.);
  if (p.type == GAUSSIAN) {
    sigma = p.width / (2. * sqrt(2. * log(2.)));
    f = exp(-1. * pow(w - p.centre, 2.) / (2. * pow(sigma, 2.)))
      * (1. / (sigma * sqrt(2. * M_PI)));
    return p.amplitude * f * gamma;
  } else {
    /* not done this yet lol */
    return NAN;
  }
}


/** convert a set of tensor subscripts to a flattened index.
 *
 */
size_t
sub2ind(std::vector<size_t> subscripts, std::vector<size_t> extents)
{
  /* std::cout << "subscripts: " << subscripts[subscripts.size() - 1] */
  /*   << ", "; */
  size_t index = subscripts[subscripts.size() - 1];
  if (subscripts[subscripts.size() - 1] 
      >= extents[extents.size() - 1]) {
    std::cout << "sub2ind: subscripts[-1] > extents[-1]" << std::endl;
    index = -1;
  }
  if (subscripts.size() != extents.size()) {
    std::cout << "sub2ind error: subscript and extents vectors are"
      " different sizes." << std::endl;
    index = -1;
  } else {
    for (size_t i = 0; i < subscripts.size() - 1; i++) {
      /* std::cout << " i = " << i; */
      if (subscripts[i] >= extents[i]) {
        std::cerr << "sub2ind error: subscripts[" << i
          << "] > extents[" << i << "]" << std::endl;
        index = -1;
      } else {
        size_t prod = 1;
        for (size_t j = i + 1; j < subscripts.size(); j++) {
          /* std::cout << " , j = " << j; */
          prod *= extents[j];
        }
        index += prod * subscripts[i];
      }
      /* std::cout << subscripts[i] << ", "; */ 
    }
  }
  /* std::cout << "index = " << index << std::endl; */
  return index;
}

/** reverse of sub2ind.
 *
 */
std::vector<size_t>
ind2sub(size_t index, std::vector<size_t> extents)
{
  std::vector<size_t> subscripts(extents.size());
  size_t total = 1;
  /* index can't be bigger than product of extents */
  for (size_t i = 0; i < extents.size(); i++) {
    total *= extents[i];
  }
  if (index >= total) {
    std::cerr << "ind2sub error: index is larger than product of"
      " dimensions given" << std::endl;
    for (size_t i = 0; i < subscripts.size(); i++) {
      subscripts[i] = -1;
    }
  } else if (index < 0) {
    std::cerr << "ind2sub error: index is < 0" << std::endl;
    for (size_t i = 0; i < subscripts.size(); i++) {
      subscripts[i] = -1;
    }
  } else {
    subscripts[extents.size() - 1] = index % extents[extents.size() - 1];
    /* declare outside loop so we can use it for the last subscript */
    /* declare as 1 here! otherwise it's uninitialised for a
     * two-index vector, causes all kinds of hell later */
    size_t prod = 1;
    for (size_t i = extents.size() - 2; i > 0; i--) {
      prod = 1;
      for (size_t j = extents.size() - 1; j > i; j--) {
        prod *= extents[j];
      }
      /* int division discards remainder */
      subscripts[i] = (index / prod) % extents[i];
    }
    /* prod is currently equal to the product of all the extents
     * from the end to the third-to-last - now multiply by e[1] */
    prod *= extents[1];
    subscripts[0] = index / prod;
  }

  return subscripts;
}

/* this won't work - was thinking about a C only version */
/*
#define EL_STRUCT(N_NORMAL, N_VIB)\
  struct e_state_#N_NORMAL_#N_VIB {\
  static double n_normal = N_NORMAL;\
  static double n_vib = N_VIB;\
  double w;\
  double w_normal[N_NORMAL];\
  double** l_ic_ij;\
  double l_ivr_i[N_VIB];\
  double** g_ic_ij;\
  double g_ivr_i[N_VIB];\
  double disp[N_VIB + 1][N_VIB + 1];\
  double fc[N_NORMAL][N_VIB + 1][N_VIB + 1];\
  } e_state_#N_NORMAL_#N_VIB;

EL_STRUCT(2,3);
*/

class Chromophore {
  public:
    Chromophore(size_t elec, size_t norm, size_t vib,
                double beta,
                double* wi, size_t n_wi,
                double* w_mode, size_t n_mode,
                double** l_ic, size_t side_l,
                double* l_ivr, size_t n_l,
                double** g_ic, size_t side_g,
                double* g_ivr, size_t n_g,
                double*** di, size_t n_mode_di, size_t e_l, size_t e_u
                );
    /* this needs changing - should the set functions be public? */
    void set_extents();
    void set_w_elec(double* wi, size_t n_wi);
    void set_w_normal(double* w_mode, size_t n_mode);
    void set_l_ic_ij(double** l_ic, size_t size);
    void set_l_ivr_i(double* l_ivr, size_t n);
    void set_g_ic_ij(double** g_ic, size_t size);
    void set_g_ivr_i(double* g_ivr, size_t n);
    void set_disp(double*** di, size_t n_mode, size_t e_l, size_t e_u);
    void fc_calc();
    void set_k_ivr(double beta);
    void set_k_ic(double beta);
    std::vector<double> dndt(double *population,
                             double t, pulse pump);
    double get_w_elec(size_t i);
    double get_w_normal(size_t i);
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

  private:
    size_t n_elec;
    size_t n_normal;
    size_t n_vib;
    std::vector<size_t> extents;
    std::vector<size_t> ivr_extents;
    std::vector<size_t> ic_extents;
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
    double beta;
    /* std::vector<double> dn;  /1* rates for IC *1/ */
};

Chromophore::Chromophore(size_t elec, size_t norm, size_t vib,
    double bet,
    double* wi, size_t n_wi,
    double* w_mode, size_t n_mode,
    double** l_ic, size_t side_l,
    double* l_ivr, size_t n_l,
    double** g_ic, size_t side_g,
    double* g_ivr, size_t n_g,
    double*** di, size_t n_mode_di, size_t e_l, size_t e_u
    )
  : n_elec(elec),
    n_normal(norm),
    n_vib(vib),
    beta(bet)
{
  set_extents();
  set_w_elec(wi, n_wi);
  set_w_normal(w_mode, n_mode);
  set_l_ic_ij(l_ic, side_l);
  set_l_ivr_i(l_ivr, n_l);
  set_g_ic_ij(g_ic, side_g);
  set_g_ivr_i(g_ivr, n_g);
  set_disp(di, n_mode_di, e_l, e_u);
  set_k_ivr(beta);
  set_k_ic(beta);
  fc_calc();
}

void
Chromophore::set_w_elec(double* wi, size_t n_wi)
{
  /* note - could specify the wi as a vector too. HMMM */
  if (n_wi != n_elec + 1) {
    std::cerr << "length of wi array does not match n_elec + 1"
              << std::endl;
  }
  for (size_t i = 0; i < n_wi; i++) {
    w_elec.push_back(wi[i]);
  }
}

void
Chromophore::set_w_normal(double* w_mode, size_t n_mode)
{
  /* note - could specify the wi as a vector too. HMMM */
  if (n_mode != n_normal) {
    std::cerr << "length of w_mode array does not match n_normal"
              << std::endl;
  }
  for (size_t i = 0; i < n_mode; i++) {
    w_normal.push_back(w_mode[i]);
  }
}

void
Chromophore::set_l_ic_ij(double** l_ic, size_t side)
{
  if (side != n_elec) {
    std::cerr << "side of l_ic_ij array does not match n_elec"
              << std::endl;
  }
  for (size_t i = 0; i < side; i++) {
    for (size_t j = 0; j < side; j++) {
      l_ic_ij.push_back(l_ic[i][j]);
    }
  }
}

void
Chromophore::set_l_ivr_i(double* l_ivr, size_t n)
{
  /* note - could specify the wi as a vector too. HMMM */
  if (n != n_elec) {
    std::cerr << "length of l_ivr array does not match n_elec"
              << std::endl;
  }
  for (size_t i = 0; i < n; i++) {
    l_ivr_i.push_back(l_ivr[i]);
  }
}

void
Chromophore::set_g_ic_ij(double** g_ic, size_t side)
{
  if (side != n_elec) {
    std::cerr << "side of g_ic_ij array does not match n_elec"
              << std::endl;
  }
  for (size_t i = 0; i < side; i++) {
    for (size_t j = 0; j < side; j++) {
      g_ic_ij.push_back(g_ic[i][j]);
    }
  }
}

void
Chromophore::set_g_ivr_i(double* g_ivr, size_t n)
{
  /* note - could specify the wi as a vector too. HMMM */
  if (n != n_elec) {
    std::cerr << "length of l_ivr array does not match n_elec"
              << std::endl;
  }
  for (size_t i = 0; i < n; i++) {
    g_ivr_i.push_back(g_ivr[i]);
  }
}

void 
Chromophore::set_disp(double*** di,
                           size_t n_mode,
                           size_t e_l,
                           size_t e_u
                          )
{
  if (n_mode != n_normal) {
    std::cerr << "mode index of disp array does not match n_normal"
              << std::endl;
  }
  if (e_l != n_elec + 1) {
    std::cerr << "e_lower for disp array does not match n_elec + 1"
              << std::endl;
  }
  if (e_u != n_elec + 1) {
    std::cerr << "e_upper for disp array does not match n_elec + 1"
              << std::endl;
  }
  for (size_t i = 0; i < n_mode; i++) {
    for (size_t j = 0; j < e_l; j++) {
      for (size_t k = 0; k < e_u; k++) {
        disp.push_back(di[i][j][k]);
      }
    }
  }
}

void
Chromophore::fc_calc()
{
  size_t index;
  /* NB: this assumes all modes have the same number of
   * vibrational levels, but we assume that everywhere else too */
  /* also for some reason this didn't work with c-style casts */
  /* double **fc_array = static_cast<double**>(malloc(n_vib + 1 * sizeof(double*))); */
  /* if (fc_array == NULL) { */
  /*   std::cout << "fc_array couldn't be allocated? top level" */
  /*     << std::endl; */
  /* } else { */
  /*   for (size_t row = 0; row < n_vib + 1; row++) { */
  /*     fc_array[row]   = static_cast<double *>(calloc(n_vib + 1, sizeof(double))); */
  /*     if (fc_array[row] == NULL) { */
  /*       std::cout << "fc_array couldn't be allocated? bottom level." */
  /*         << " row = " << row << std::endl; */
  /*     } */
  /*   } */
  /* } */
  /* double **fc_array = new double*[n_vib + 1]; */
  /* for (size_t row = 0; row < n_vib + 1; row++) { */
  /*   fc_array[row] = new double[n_vib + 1]; */
  /* } */
  std::vector<std::vector<double>> fc_array(n_vib + 1,
              std::vector<double>(n_vib + 1, 0.));

  for (size_t i = 0; i < n_elec + 1; i++) {
    /* std::cout << "i = " << i; */
    for (size_t j = 0; j < n_elec + 1; j++) {
      if (i != j) {

        /* std::cout << ", j = " << j; */
        for (size_t alpha = 0; alpha < n_normal; alpha++) {
          double displacement = get_disp({alpha, i, j});

          fc_array[0][0] = exp(-0.5 * pow(displacement, 2.));

          for (size_t a = 1; a < n_vib + 1; a++) {
            fc_array[0][a] = ((- 1. * displacement) /
                           sqrt((double)a)) * fc_array[0][a - 1];
            fc_array[a][0] = (double)(pow(-1., a)) * fc_array[0][a];
          }

          for (size_t a = 1; a < n_vib + 1; a++) {
            for (size_t b = a; b < n_vib + 1; b++) {
              fc_array[a][b] = ((-1. * displacement) /  
                (sqrt((double)b)) * fc_array[a][b - 1])
              + (sqrt((double)a)/sqrt((double)b)) * fc_array[a - 1][b - 1];
              fc_array[b][a] = (double)(pow(-1., (double)(b - a))) * fc_array[a][b];
            }
          }

          /* this is how chris's code does it - weird. ask */
          for (size_t a = 0; a < n_vib + 1; a++) {
            for (size_t b = 0; b < n_vib + 1; b++) {
              std::vector<size_t> subscripts = {i, j, alpha, a, b};
              index = sub2ind(subscripts, fc_extents);
              fc[index] = fc_array[a][b];
            }
          }

        }

      }
    }
  }
  /* free(fc_array); */
  /* for (size_t row = 0; row < n_vib + 1; row++) { */
  /*   delete [] fc_array[row]; */
  /* } */
  /* delete [] fc_array; */

}

void
Chromophore::set_extents()
{
  size_t size = 1;
  extents.push_back(n_elec);
  population_extents.push_back(n_elec);
  extents.push_back(n_normal);
  ic_extents = {n_elec, n_elec};
  size *= n_elec * n_elec;
  for (size_t i = 0; i < n_normal; i++) {
    extents.push_back(n_vib);
    population_extents.push_back(n_vib + 1);
    ic_extents.push_back(n_vib + 1);
    ic_extents.push_back(n_vib + 1);
    size *= (n_vib + 1);
  }
  k_ic.resize(size, 0.);
  size = 1;

  ivr_extents.push_back(n_elec);
  ivr_extents.push_back(n_normal);
  ivr_extents.push_back(2);
  size = n_elec * n_normal * 2;
  /* k_ivr.resize(size, 0.); */

  fc_extents.push_back(n_elec + 1);
  fc_extents.push_back(n_elec + 1);
  fc_extents.push_back(n_normal);
  fc_extents.push_back(n_vib + 1);
  fc_extents.push_back(n_vib + 1);
  size = pow(n_elec + 1, 2) * pow(n_vib + 1, 2) * n_normal;
  fc.resize(size, 0.);

}

/* getting things */

double
Chromophore::get_w_elec(size_t elec_i)
{
  if (elec_i > n_elec) {
    std::cout << "invalid index passed to get_w_elec: "
    << elec_i << " - max is " << n_elec + 1 << std::endl;
    return NAN;
  }
  return w_elec[elec_i];
}

std::vector<double>
Chromophore::get_w_normal()
{
  return w_normal;
}

double
Chromophore::get_w_normal(size_t norm)
{
  if (norm >= n_normal) {
    std::cout << "invalid index passed to get_w_elec: "
    << norm << " - max is " << n_normal << std::endl;
    return NAN;
  }
  return w_normal[norm];
}

double
Chromophore::get_l_ic_ij(size_t elec_i, size_t elec_j)
{
  return l_ic_ij[elec_i + (n_elec * elec_j)];
}

double
Chromophore::get_l_ivr_i(size_t elec_i)
{
  return l_ivr_i[elec_i];
}

double
Chromophore::get_g_ic_ij(size_t elec_i, size_t elec_j)
{
  return g_ic_ij[elec_i + (n_elec * elec_j)];
}

double
Chromophore::get_g_ivr_i(size_t elec_i)
{
  return g_ivr_i[elec_i];
}

double
Chromophore::get_disp(std::vector<size_t> subscripts)
{
  /* the ordering's a bit off here - extents lists the electronic
   * extents first, then the normal modes, whereas the disp array
   * indexes them with the mode first. that could be changed, then
   * this would just be return disp[sub2ind()] etc */
  size_t index = sub2ind(subscripts, {n_normal, n_elec + 1, n_elec + 1});
  return disp[index];
}

double
Chromophore::get_fc(std::vector<size_t> subscripts)
{
  return fc[sub2ind(subscripts, fc_extents)];
}

/* this could be replaced with the existing cw_obo but
 * would require something like (void *) params etc. */
double
C_OBO(double w, double l, double g)
{
  if (w != INFINITY) {
    return 2.0 * l * g * w / (pow(w, 2.) + pow(g, 2.));
  } else {
    return 0.0;
  }
}


/** Vibrational relaxation rates on the same mode/electronic state. 
 *
 * wrapper function because we'll use this multiple times.
 * Computed by (secular, Markovian) Redfield theory.
 * w is the frequency (e.g of normal mode or energy gap between states).
 * lambda and gamma are parameters for the spectral density.
 *
 */
double
k_calc(double w, double beta, double lambda, double gamma)
{
  return C_OBO(w, lambda, gamma) * ((1.0 / tanh(0.5 * beta * w)) + 1);
}

/** Interconversion rate between vibronic states on different manifolds.
 *
 * e_ij and w_e are 2-vectors; e_ij[0] is state i, e_ij[1] is state j.
 * w_e are their respective energies.
 * a and b are tuples of vibrational quantum numbers on each normal mode
 * and w_alpha are the frequencies of the normal modes - these should all
 * have .size() == n_normal.
 */
double 
k_inter(std::vector<size_t> e_ij,
    std::vector<size_t> a, std::vector<size_t> b,
    std::vector<double> w_e, std::vector<double> w_alpha,
    double beta, double lambda, double gamma)
{
  double k;
  double w_ab = 0.;
  if ((a.size() != b.size()) || (a.size() != w_alpha.size())) {
    std::cerr << "size mismatch in k_inter: a.size() = "
      << a.size() << ", b.size() = " << b.size()
      << ", w_alpha.size() = " << w_alpha.size() << std::endl;
  }
  for (size_t i = 0; i < a.size(); i++) {
    /* note - casting int here might be a problem if
     * a[i] or b[i] get very large. but i think we're fine */
    w_ab += (w_alpha[i] * (int)(a[i] - b[i]));
  }
  double w_ij = fabs(w_e[0] - w_e[1]);
  double delta_ij_ab = w_ij + w_ab;
  if (e_ij[0] > e_ij[1]) {
    k = k_calc(-delta_ij_ab, beta, lambda, gamma);
  } else {
    k = k_calc(delta_ij_ab, beta, lambda, gamma);
  }
  return k;
}

double thermal_osc(std::vector<size_t> a,
                   std::vector<double> w_alpha, double beta)
{
  double eps_a = 0.;
  if (a.size() != w_alpha.size()) {
    std::cerr << "thermal_osc error: number of quantum numbers a = "
      << a.size() << ", number of normal modes = " << w_alpha.size()
      << ", should be identical." << std::endl;
  }
  for (size_t i = 0; i < a.size(); i++) {
    eps_a += w_alpha[i] * (a[i] + 0.5);
  }
  double n_a = exp(- beta * eps_a);
  double Z = 1.;
  for (size_t i = 0; i < w_alpha.size(); i++) {
    Z *= exp(-0.5 * beta * w_alpha[i]) / (1 - exp(- beta * w_alpha[i]));
  }
  return n_a / Z;
}

void
Chromophore::set_k_ivr(double beta)
{
  for (size_t i = 0; i < n_elec; i++) {
    for (size_t alpha = 0; alpha < n_normal; alpha++) {
      /* NB: python code puts these in backwards - the downward rate
       * first then the forward one, or so i thought.
       * this is the same ordering but check it!!
       */
      k_ivr.push_back(CM_PER_PS *
          k_calc(w_normal[alpha], beta, l_ivr_i[i], g_ivr_i[i]));
      k_ivr.push_back(CM_PER_PS *
          k_calc(-1. * w_normal[alpha], beta, l_ivr_i[i], g_ivr_i[i]));
    }
  }
}

void
Chromophore::set_k_ic(double beta)
{
  std::vector<size_t> vib_extents;
  size_t total_vib_rates = 1;
  for (size_t i = 0; i < n_normal; i++) {
    /* one for acceptor, one for donor */
    vib_extents.push_back(n_vib + 1);
    vib_extents.push_back(n_vib + 1);
    total_vib_rates *= pow((n_vib + 1), 2);
  }
  /* make it the correct length. can we zero it out here? */
  k_ic.resize(pow(n_elec, 2) * total_vib_rates, 0.);

  /* std::cout << "k_ic:" << std::endl; */
  for (size_t i = 0; i < n_elec; i++) {
    for (size_t j = 0; j < n_elec; j++) {

      if (i != j) {
        std::vector<size_t> e_ij = {i, j};
        size_t ij_index = sub2ind(e_ij, {n_elec, n_elec});
        std::vector<double> w_e = {w_elec[i], w_elec[j]};
        /* chris did this by creating an array with dimensions
         * k_ic_extents and then making an iterator over the last
         * (2 * n_normal) dimensions of that array. but i feel like
         * there's an easier way of doing it */
        for (size_t vib = 0; vib < total_vib_rates; vib++) {
          std::vector<size_t> subs = ind2sub(vib, vib_extents);
          if (subs.size() != 2 * n_normal) {
            std::cerr << "set_k_ic: subscript length != 2 * n_normal"
              << std::endl;
          } else {
            std::vector<size_t> a, b;
            for (size_t k = 0; k < subs.size(); k++) {
              if (k < subs.size()/2 ) {
                a.push_back(subs[k]);
              } else {
                b.push_back(subs[k]);
              }
            }
            /* now push the electronic indices onto the front of subs */
            subs.insert(subs.begin(), j);
            subs.insert(subs.begin(), i);
            size_t index = sub2ind(subs, ic_extents);
            k_ic[index] = k_inter(e_ij, a, b, w_e, w_normal, beta,
                l_ic_ij[ij_index],
                g_ic_ij[ij_index]) * CM_PER_PS;
          }
        }

      }
    }
  }
}

std::vector<double>
Chromophore::dndt(double *population, double t, pulse pump)
{
  std::vector<size_t> n_extents;
  std::vector<size_t> vib_extents;
  n_extents.push_back(n_elec);
  size_t n_total = n_elec;
  for (size_t i = 0; i < n_normal; i++) {
    n_extents.push_back(n_vib + 1);
    vib_extents.push_back(n_vib + 1);
    n_total *= (n_vib + 1);
  }
  
  std::vector<double> dndt_ivr(n_total, 0.);
  std::vector<double> dndt_ic(n_total, 0.);
  std::vector<double> dndt_pump(n_total, 0.);
  std::vector<double> dndt(n_total, 0.);

  std::vector<size_t> subscripts(n_extents.size(), 0);
  std::vector<size_t> sub_lower(n_extents.size(), 0),
                      sub_upper(n_extents.size(), 0);
  size_t i_lower = 0, i_upper = 0;
  for (size_t i = 0; i < n_total; i++) {
    subscripts = ind2sub(i, n_extents);  
     
    /* IVR. can probably be optimised */
    for (size_t alpha = 0; alpha < n_normal; alpha++) {
      if (subscripts[alpha + 1] > 0) {
        /* loss of population to lower vibronic state */
        dndt_ivr[i] -= subscripts[alpha + 1]
                    * k_ivr[sub2ind({subscripts[0], alpha, 0},
                      {n_elec, n_normal, 2})]
                    * population[i];
        sub_lower = subscripts;
        sub_lower[alpha + 1]--;
        i_lower = sub2ind(sub_lower, n_extents);

        /* gain of population from level below */
        dndt_ivr[i] += subscripts[alpha + 1]
                    * k_ivr[sub2ind({subscripts[0], alpha, 1},
                      {n_elec, n_normal, 2})]
                    * population[i_lower];
      }

      if (subscripts[alpha + 1] < n_vib) {
        /* loss of population to upper state */
        dndt_ivr[i] -= (subscripts[alpha + 1] + 1.)
                    * k_ivr[sub2ind({subscripts[0], alpha, 1},
                      {n_elec, n_normal, 2})]
                    * population[i];
        sub_upper = subscripts;
        sub_upper[alpha + 1]++;
        i_upper = sub2ind(sub_upper, n_extents);
        /* gain from the upper state */
        dndt_ivr[i] += (subscripts[alpha + 1] + 1.)
                    * k_ivr[sub2ind({subscripts[0], alpha, 0},
                      {n_elec, n_normal, 2})]
                    * population[i_upper];
      }
    }

  }

  /* interconversion */
  /* can probably be folded into above loop but
   * i haven't figured out how yet exactly */

  for (size_t i = 0; i < n_elec; i++) {
    for (size_t j = 0; j < (pow((n_vib + 1), n_normal)); j++) {
      std::vector<size_t> a = ind2sub(j, vib_extents);
      for (size_t k = 0; k < (pow((n_vib + 1), n_normal)); k++) {
        std::vector<size_t> b = ind2sub(k, vib_extents);

        /* indices */
        std::vector<size_t> n_i, n_i_plus, n_i_minus;
        std::vector<size_t> k_ba_ji, k_ba_ij, k_ab_ji, k_ab_ij;
        k_ba_ji = {i + 1, i};
        k_ba_ij = {i, i + 1};
        k_ab_ji = {i - 1, i};
        k_ab_ij = {i, i - 1};
        n_i.push_back(i); 
        n_i_plus.push_back(i + 1); 
        n_i_minus.push_back(i - 1); 
        for (size_t ai = 0; ai < a.size(); ai++) {
          n_i.push_back(a[ai]);
          k_ab_ji.push_back(a[ai]);
          k_ab_ij.push_back(a[ai]);
        }
        for (size_t bi = 0; bi < b.size(); bi++) {
          n_i_plus.push_back(b[bi]);
          n_i_minus.push_back(b[bi]);
          k_ba_ji.push_back(b[bi]);
          k_ba_ij.push_back(b[bi]);
        }

        /* two sets of vibrational indices to worry about */
        for (size_t ai = 0; ai < a.size(); ai++) {
          k_ba_ji.push_back(a[ai]);
          k_ba_ij.push_back(a[ai]);
        }
        for (size_t bi = 0; bi < b.size(); bi++) {
          k_ab_ji.push_back(b[bi]);
          k_ab_ij.push_back(b[bi]);
        }

        if (i == 0) { /* electronic ground state */
          double fc_i_plus = 1.;
          for (size_t alpha = 0; alpha < n_normal; alpha++) {
            fc_i_plus *= pow(fc[sub2ind({i, i + 1, alpha, 
                      a[alpha], b[alpha]}, fc_extents)], 2.);
          }

          dndt_ic[sub2ind(n_i, n_extents)] -= fc_i_plus
          * k_ic[sub2ind(k_ba_ji, ic_extents)]
          * population[sub2ind(n_i, n_extents)];
          dndt_ic[sub2ind(n_i, n_extents)] += fc_i_plus
          * k_ic[sub2ind(k_ba_ij, ic_extents)]
          * population[sub2ind(n_i_plus, n_extents)];

        } else if (i == n_elec - 1) { /* highest electronic state */
          double fc_i_minus = 1.;
          for (size_t alpha = 0; alpha < n_normal; alpha++) {
            fc_i_minus *= pow(fc[sub2ind({i - 1, i, alpha, 
                      b[alpha], a[alpha]}, fc_extents)], 2.);
          }

          dndt_ic[sub2ind(n_i, n_extents)] -= fc_i_minus
          * k_ic[sub2ind(k_ab_ji, ic_extents)]
          * population[sub2ind(n_i, n_extents)];
          dndt_ic[sub2ind(n_i, n_extents)] += fc_i_minus
          * k_ic[sub2ind(k_ab_ij, ic_extents)]
          * population[sub2ind(n_i_minus, n_extents)];
        } else { /* intermediate */
          double fc_i_plus = 1.; double fc_i_minus = 1.;
          for (size_t alpha = 0; alpha < n_normal; alpha++) {
            fc_i_plus *= pow(fc[sub2ind({i, i + 1, alpha, 
                      a[alpha], b[alpha]}, fc_extents)], 2.);
            fc_i_minus *= pow(fc[sub2ind({i - 1, i, alpha, 
                      b[alpha], a[alpha]}, fc_extents)], 2.);
          }

          dndt_ic[sub2ind(n_i, n_extents)] -= fc_i_plus
          * k_ic[sub2ind(k_ba_ji, ic_extents)]
          * population[sub2ind(n_i, n_extents)];
          dndt_ic[sub2ind(n_i, n_extents)] += fc_i_plus
          * k_ic[sub2ind(k_ba_ij, ic_extents)]
          * population[sub2ind(n_i_plus, n_extents)];

          dndt_ic[sub2ind(n_i, n_extents)] -= fc_i_minus
          * k_ic[sub2ind(k_ab_ji, ic_extents)]
          * population[sub2ind(n_i, n_extents)];
          dndt_ic[sub2ind(n_i, n_extents)] += fc_i_minus
          * k_ic[sub2ind(k_ab_ij, ic_extents)]
          * population[sub2ind(n_i_minus, n_extents)];

        } /* end of electronic state if/else */

      }
    }
  }

  /* pumping. again this seems to be duplication of above loops */
  for (size_t j = 0; j < (n_normal * (n_vib + 1)); j++) {
    std::vector<size_t> a = ind2sub(j, vib_extents);
    for (size_t k = 0; k < (n_normal * (n_vib + 1)); k++) {
      std::vector<size_t> b = ind2sub(k, vib_extents);

      double delta_x0_ba = w_elec[pump.target_state] - w_elec[0];
      double fc_sq = 1.;
      for (size_t alpha = 0; alpha < n_normal; alpha++) {
        delta_x0_ba += w_normal[alpha] * (int)(b[alpha] - a[alpha]);
        fc_sq *= pow(fc[sub2ind({0,
              pump.target_state, alpha, a[alpha], b[alpha]},
              fc_extents)], 2.);
      }

      std::vector<size_t> gs_sub = {0};
      std::vector<size_t> exc_sub = {pump.target_state};
      for (size_t ai = 0; ai < a.size(); ai++) {
        gs_sub.push_back(a[ai]);
        exc_sub.push_back(b[ai]);
      }
      dndt_pump[sub2ind(gs_sub, n_extents)]
            -= intensity(delta_x0_ba, t, pump) * fc_sq
            * population[sub2ind(gs_sub, n_extents)];
      dndt_pump[sub2ind(exc_sub, n_extents)]
            += intensity(delta_x0_ba, t, pump) * fc_sq
            * population[sub2ind(gs_sub, n_extents)];
    }
  }

  double ivr_sum = 0., ic_sum = 0., pump_sum = 0.;
  for (size_t i = 0; i < n_total; i++) {
    ivr_sum  += dndt_ivr[i];
    ic_sum   += dndt_ic[i];
    pump_sum += dndt_pump[i];
    dndt[i] = dndt_ivr[i] + dndt_ic[i] + dndt_pump[i];
    /* fprintf(stdout, "%5lu  %18.10f  %18.10f  %18.10f  %18.10f\n", */
    /*     i, dndt_ivr[i], dndt_ic[i], dndt_pump[i], dndt[i]); */ 
  }
  /* fprintf(stdout, "%18.10f  %18.10f  %18.10f\n", ivr_sum, ic_sum, pump_sum); */ 
  return dndt;
}

/* OTHER THOUGHTS:
 *
 * for populations we need a set of [N_chromophores][extents]
 *
 */

typedef struct vera_lsoda_data {
  Chromophore *chromo;
  /* std::vector<double> *population; */
  pulse *pump;
} vera_lsoda_data;


void
func(double t, double *y, double *ydot, void *data)
{
  vera_lsoda_data *vls = (vera_lsoda_data *)data;
  Chromophore *chromo = vls->chromo;
  std::vector<double> ydot_vec = chromo->dndt(y, t, (*vls->pump));
  size_t i;
  for (i = 0; i < ydot_vec.size(); i++) {
    ydot[i] = ydot_vec[i];
  }
}
 

int
main(int argc, char** argv)
{

  chrono::steady_clock::time_point begin = chrono::steady_clock::now();

  size_t n_elec = 3, n_normal = 2, n_vib = 3;
  double beta = 1./200.; /* temp in wavenumbers i think */

  /* these will need to be calloc'd and filled from a file eventually */
  double *wi = (double *)calloc(n_elec + 1, sizeof(double));
  wi[0] = 0.; wi[1] = 10000.0; wi[2] = 20000.0; wi[3] = 30000.0;

  double *wnorm = (double *)calloc(n_normal, sizeof(double));
  wnorm[0] = 1100.; wnorm[1] = 1500.0;

  double **licij = (double **)calloc(n_elec, sizeof(double**));
  if (licij == NULL) {
    std::cout << "licij array calloc didn't work" << std::endl;
  } else {
    for (size_t i = 0; i < n_elec; i++) {
      licij[i] = (double *)calloc(n_elec, sizeof(double*));
      if (licij[i] == NULL) {
        std::cout << "licij[" << i <<"] calloc didn't work" << std::endl;
      }
    }
  }
  licij[0][1] = 200.;
  licij[1][2] = 1000.;

  double *livri = (double *)calloc(n_elec, sizeof(double));
  if (livri == NULL) {
    std::cout << "livri array calloc didn't work" << std::endl;
  } else {
    livri[0] = 50.0; livri[1] = 100.0; livri[2] = 300.0;
  }

  /* these need to be converted from fs to cm^-1 */
  double **gicij = (double **)calloc(n_elec, sizeof(double**));
  if (gicij == NULL) {
    std::cout << "gicij array calloc didn't work" << std::endl;
  } else {
    for (size_t i = 0; i < n_elec; i++) {
      gicij[i] = (double *)calloc(n_elec, sizeof(double*));
      if (licij[i] == NULL) {
        std::cout << "gicij[" << i <<"] calloc didn't work" << std::endl;
      }
    }
  }
  gicij[0][1] = 163.6;
  gicij[1][2] = 163.6;

  double *givri = (double *)calloc(n_elec, sizeof(double));
  if (givri == NULL) {
    std::cout << "livri array calloc didn't work" << std::endl;
  } else {
    givri[0] = 163.6; givri[1] = 163.6; givri[2] = 163.6;
  }
  for (size_t i = 0; i < n_elec; i++) {
    /* the damping times were given in femtoseconds */
    givri[i] = 1. / (givri[i] * 1E-3 * CM_PER_PS);
    for (size_t j = 0; j < n_elec; j++) {
      if (gicij[i][j] != 0) {
        gicij[i][j] = 1. / (gicij[i][j] * 1E-3 * CM_PER_PS);
      }
    }
  }
  /* for (size_t i = 0; i < n_elec; i++) { */
  /*   std::cout << givri[i] << std::endl; */
  /* } */
  /* for (size_t i = 0; i < n_elec; i++) { */
  /*   for (size_t j = 0; j < n_elec; j++) { */
  /*     std::cout << gicij[i][j] << std::endl; */
  /*   } */
  /* } */

  double ***disp = (double ***)calloc(n_normal, sizeof(*disp));
  if (disp == NULL) {
    std::cout << "disp array calloc didn't work" << std::endl;
  } else {
    for (size_t i = 0; i < n_normal; i++) {
      disp[i] = (double **)calloc(n_elec + 1, sizeof(*(disp[i])));
      if (disp[i] == NULL) {
        std::cout << "disp[" << i <<"] calloc didn't work" << std::endl;
      } else {
        for (size_t j = 0; j < n_elec + 1; j++) {
          disp[i][j] = (double *)calloc(n_elec + 1, sizeof(*(disp[i][j])));
          if (disp[i][j] == NULL) {
            std::cout << "disp[" << i <<"][" << j 
              << "] calloc didn't work" << std::endl;
          }
        }
      }
    }
  }
  disp[0][0][1] = 1.;
  disp[1][0][1] = 1.;
  disp[0][1][2] = 1.;
  disp[1][1][2] = 1.;
  disp[0][0][2] = 0.7;
  disp[1][0][2] = 0.7;
  disp[0][1][3] = 0.4;
  disp[1][1][3] = 0.4;

  Chromophore lutein = Chromophore(n_elec, n_normal, n_vib,
      beta, wi, n_elec + 1,
      wnorm, n_normal,
      licij, n_elec,
      livri, n_elec,
      gicij, n_elec,
      givri, n_elec,
      disp, n_normal, n_elec + 1, n_elec + 1
      );

  for (size_t i = 0; i < n_normal; i++) {
    for (size_t j = 0; j < n_elec + 1; j++) {
      free(disp[i][j]);
    }
  }
  for (size_t i = 0; i < n_elec; i++) {
    free(gicij[i]);
    free(licij[i]);
    if (i < n_normal) {
      free(disp[i]);
    }
  }
  free(disp);
  free(givri);
  free(gicij);
  free(livri);
  free(licij);
  free(wnorm);
  free(wi);

  pulse pump = {
    .type = GAUSSIAN,
    .target_state = 2,
    .amplitude = 100.,
    .centre = 20000.,
    .width = 10.,
    .t_peak = 0.1,
    .duration = 70.0E-3,
  };

  /* this should be available via c.get_n_extents() or something */
  /* NB: this whole thing should be a method - calculate initial
   * populations; only parameter is beta */
  std::vector<size_t> pop_extents = {n_elec};
  std::vector<size_t> vib_extents;

  size_t pop_total = n_elec;
  size_t vib_total = 1;
  for (size_t i = 0; i < n_normal; i++) {
    pop_extents.push_back(n_vib + 1);
    vib_extents.push_back(n_vib + 1);
    pop_total *= n_vib + 1;
    vib_total *= n_vib + 1;
  }
  double *pop = (double *)calloc(pop_total, sizeof(double));
  std::vector<double> n0(pop_total, 0.);
  /* vec is because LSODA requires vector argument
   * - this could maybe be fixed somehow with vector.data()? */
  std::vector<double> y(pop_total, 0.);

  for (size_t i = 0; i < vib_total; i++) {
    std::vector<size_t> tot_subs(pop_extents.size(), 0);
    std::vector<size_t> vib_subs(vib_extents.size(), 0);
    vib_subs = ind2sub(i, vib_extents);
    for (size_t j = 0; j < vib_subs.size(); j++) {
      tot_subs[j + 1] = vib_subs[j];
    }
    size_t index = sub2ind(tot_subs, pop_extents);
    n0[index] = thermal_osc(vib_subs, lutein.get_w_normal(), beta);
    pop[index] = n0[index];
    y[index] = pop[index];
  }

  /* for (size_t i = 0; i < pop_total; i++) { */
  /*   fprintf(stdout, "%lu\t%18.10f\n", i, pop[i]); */
  /* } */

  vera_lsoda_data *vls = (vera_lsoda_data *)malloc(sizeof(vera_lsoda_data));
  vls->chromo = &lutein;
  vls->pump = &pump;
  void *data = (void *)vls;

  LSODA lsoda;
  std::vector<double> yout;
  int istate = 1;

  size_t n_steps = 10000;
  double ti = 0., tf = 10., dt = (tf - ti) / n_steps;

  double t = ti, tout = dt;
  double ysum = 0.;


  for (size_t i = 0; i < n_steps; i++) {

    std::cout << t << ' ' << setprecision(8) << y[0] << ' ' 
              << setprecision(8) 
              << y[sub2ind({1, 0, 0}, pop_extents)] << ' ' 
              << setprecision(8) 
              << y[sub2ind({2, 0, 0}, pop_extents)] << ' '
              << setprecision(8) 
              << y[sub2ind({0, 0, 1}, pop_extents)] << ' ' 
              << setprecision(8) 
              << y[sub2ind({0, 1, 0}, pop_extents)] << ' ' 
              << setprecision(8) << ysum
              << std::endl;

    ysum = 0.;
    /* 1.49e-8 is the default value for rtol and atol
     * in the scipy odeint, which also uses lsoda afaik.
     * i was trying to match the output from that but
     * the tolerances don't seem to make much difference? */
    lsoda.lsoda_update(func, pop_total, y, yout,
        &t, tout, &istate, data, 1.49e-8, 1.49e-8);
    tout += dt;

    for (size_t j = 0; j < pop_total; j++) {
      y[j] = yout[j + 1];
      ysum += y[j];
    }
    /* std::cout << "Step " << i << ": " << std::endl; */
    /* if (fabs(ysum - 1.) > 1E-4) { */
    /*   fprintf(stdout, "Sum of populations .ne. 1: ysum = %12.8f," */
    /*       " step number = %5lu\n", ysum, i); */
    /* } */

    if (istate <= 0)
    {
       std::cerr << "error istate = " <<  istate << std::endl;
       throw runtime_error( "Failed to compute the solution." );
    }


  }

  free(pop);
  free(vls);

  chrono::steady_clock::time_point end = chrono::steady_clock::now();
  std::cout << "|| Time taken (ms)= " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << std::endl;

  return 0;

}

#endif
