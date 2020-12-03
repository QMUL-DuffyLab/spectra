#include <iostream>
#include <vector>
#include <cmath>
/* #include "forster.h" */

#ifndef __VER_H__
#define __VER_H__

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
                double* wi, size_t n_wi,
                double* w_mode, size_t n_mode,
                double** l_ic, size_t side_l,
                double* l_ivr, size_t n_l,
                double** g_ic, size_t side_g,
                double* g_ivr, size_t n_g,
                double*** di, size_t n_mode_di, size_t e_l, size_t e_u
                );
    void set_w_elec(double* wi, size_t n_wi);
    void set_w_normal(double* w_mode, size_t n_mode);
    void set_l_ic_ij(double** l_ic, size_t size);
    void set_l_ivr_i(double* l_ivr, size_t n);
    void set_g_ic_ij(double** g_ic, size_t size);
    void set_g_ivr_i(double* g_ivr, size_t n);
    void set_disp(double*** di, size_t n_mode, size_t e_l, size_t e_u);
    void fc_calc();
    double get_w_elec(size_t i);
    double get_w_normal(size_t i);
    double get_w_vib(size_t i);
    double get_l_ic_ij(size_t i, size_t j);
    double get_l_ivr_i(size_t i);
    double get_g_ic_ij(size_t i, size_t j);
    double get_g_ivr_i(size_t i);
    /* The indexing is disp[mode][elec_lower][elec_upper] */
    double get_disp(size_t mode, size_t e_lower, size_t e_upper);
    double get_fc(
             size_t elec_i,
             size_t elec_j,
             size_t mode,
             size_t vib_a,
             size_t vib_b
             );

  private:
    size_t n_elec;
    size_t n_normal;
    size_t n_vib;
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
};

Chromophore::Chromophore(size_t elec, size_t norm, size_t vib,
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
    n_vib(vib)
{
  set_w_elec(wi, n_wi);
  set_w_normal(w_mode, n_mode);
  set_l_ic_ij(l_ic, side_l);
  set_l_ivr_i(l_ivr, n_l);
  set_g_ic_ij(g_ic, side_g);
  set_g_ivr_i(g_ivr, n_g);
  set_disp(di, n_mode_di, e_l, e_u);
  fc_calc();
}

void Chromophore::set_w_elec(double* wi, size_t n_wi)
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

void Chromophore::set_w_normal(double* w_mode, size_t n_mode)
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

void Chromophore::set_l_ic_ij(double** l_ic, size_t side)
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

void Chromophore::set_l_ivr_i(double* l_ivr, size_t n)
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

void Chromophore::set_g_ic_ij(double** g_ic, size_t side)
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

void Chromophore::set_g_ivr_i(double* g_ivr, size_t n)
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

void Chromophore::set_disp(double*** di,
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

void Chromophore::fc_calc()
{
  double **fc_array = (double **)calloc(n_vib, sizeof(double*));
  for (size_t row = 0; row < n_vib + 1; row++) {
    fc_array[row]   = (double *) calloc(n_vib, sizeof(double));
  }

  for (size_t i = 0; i < n_elec + 1; i++) {
    for (size_t j = 0; j < n_elec + 1; j++) {
      if (i != j) {

        for (size_t alpha = 0; alpha < n_normal; alpha++) {
          double displacement = get_disp(alpha, i, j);

          fc_array[0][0] = exp(-0.5 * pow(displacement, 2.));

          for (size_t a = 1; a < n_vib + 1; a++) {

            fc_array[0][a] = (- 1. * displacement) /
                           sqrt((double)a * fc_array[0][a - 1]);
            fc_array[a][0] = (pow(-1.0, a)) * fc_array[0][a];

            for (size_t b = 1; b < n_vib + 1; b++) {

              fc_array[a][b] = (- 1. * displacement) /  
                (sqrt((double)b) * fc_array[a][b - 1])
              + (sqrt((double)a)/sqrt((double)b)) * fc_array[a - 1][b - 1];
              fc_array[b][a] = (pow(-1.0, b - a)) * fc_array[a][b];

              fc.push_back(fc_array[a][b]);
            }
          }
        }

      } else {
        /* this is to ensure fc is of the right length so the
         * get_fc function works properly re: indexing */
        for (size_t dummy = 0;
             dummy < n_normal * n_vib + 1 * n_vib + 1;
             dummy++) {
          fc.push_back(0.0);
        }
      }
    }
  }
  free(fc_array);
}

/* getting things */

double Chromophore::get_w_elec(size_t elec_i)
{
  if (elec_i > n_elec) {
    std::cout << "invalid index passed to get_w_elec: "
    << elec_i << " - max is " << n_elec + 1 << std::endl;
    return NAN;
  }
  return w_elec[elec_i];
}

double Chromophore::get_w_normal(size_t norm)
{
  if (norm >= n_normal) {
    std::cout << "invalid index passed to get_w_elec: "
    << norm << " - max is " << n_normal << std::endl;
    return NAN;
  }
  return w_normal[norm];
}

/*
double Chromophore::get_w_vib(size_t vib_a)
{
  return w_vib[vib_a];
}
*/

double Chromophore::get_l_ic_ij(size_t elec_i, size_t elec_j)
{
  return l_ic_ij[elec_i + (n_elec * elec_j)];
}

double Chromophore::get_l_ivr_i(size_t elec_i)
{
  return l_ivr_i[elec_i];
}

double Chromophore::get_g_ic_ij(size_t elec_i, size_t elec_j)
{
  return g_ic_ij[elec_i + (n_elec * elec_j)];
}

double Chromophore::get_g_ivr_i(size_t elec_i)
{
  return g_ivr_i[elec_i];
}

double Chromophore::get_disp(
             size_t mode,
             size_t e_lower,
             size_t e_upper
)
{
  size_t index = (mode   + n_normal * 
                 (e_lower  + (n_elec + 1)  * e_upper));
  return disp[index];
}

double Chromophore::get_fc(
             size_t elec_i,
             size_t elec_j,
             size_t mode,
             size_t vib_a,
             size_t vib_b
)
{
  size_t index = elec_i + (n_elec + 1) *
                (elec_j + (n_elec + 1) * 
                (mode   + n_normal * 
                (vib_a  + (n_vib + 1)  * vib_b)));
  return fc[index];
}

/* this could be replaced with the existing cw_obo but
 * would require something like (void *) params etc. */
double C_OBO(double w, double l, double g)
{
  if (w != INFINITY) {
    return 2.0 * l * g * w / (pow(w, 2.) + pow(g, 2.));
  } else {
    return 0.0;
  }
}

/* wrapper function because we'll use this multiple times */
double k_calc(double w, double beta, double lambda, double gamma)
{
  return C_OBO(w, lambda, gamma) * (1.0 / (tanh(0.5 * beta * w)) + 1);
}

/** Vibrational relaxation rates on the same mode/electronic state. 
 *
 * Computed by (secular, Markovian) Redfield theory.
 * w_alpha is the frequency of normal mode alpha.
 * lambda and gamma are parameters for the spectral density.
 * returns (k_up, k_down)
 */
std::vector<double> k_alpha(double w_alpha, double beta,
                            double lambda, double gamma)
{
  std::vector<double> k;
  k.push_back(k_calc(-1. * w_alpha, beta, lambda, gamma));
  k.push_back(k_calc(w_alpha, beta, lambda, gamma));
  /* k.push_back(C_OBO((-1.) * w_alpha, lambda, gamma) * */ 
  /*            (1.0 / (tanh(-0.5 * beta * w_alpha)) + 1)); */
  /* k.push_back(C_OBO(w_alpha, lambda, gamma) * */ 
  /*            (1.0 / (tanh(0.5 * beta * w_alpha)) + 1)); */
  return k;
}

/** Interconversion rate between vibronic states on different manifolds.
 *
 * e_ij and w_e are 2-vectors; e_ij[0] is state i, e_ij[1] is state j.
 * w_e are their respective energies.
 * a and b are tuples of vibrational quantum numbers on each normal mode
 * and w_alpha are the frequencies of the normal modes - these should all
 * have .size() == n_normal.
 */
double k_inter(std::vector<size_t> e_ij,
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
    w_ab += (w_alpha[i] * (a[i] - b[i]));
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

/** convert a set of tensor subscripts to a flattened index.
 *
 * haven't tested it, just wrote it down on pen and paper
 * so hopefully it's correct lmao
 */
size_t sub2ind(std::vector<size_t> subscripts, std::vector<size_t> extents)
{
  size_t index = subscripts[0];
  if (subscripts.size() != extents.size()) {
    std::cerr << "sub2ind error: subscript and extents vectors are"
      " different sizes." << std::endl;
    index = -1;
  } else {
    for (size_t i = 1; i < subscripts.size(); i++) {
      size_t prod = 1;
      for (size_t j = 0; j < i; j++) {
        prod *= extents[j];
      }
      index += prod * subscripts[i];
    }
  }
  return index;
}

/** reverse of sub2ind.
 *
 * again i haven't checked this yet, just spent ten minutes
 * writing it out by hand and i think this is how it should look
 */
std::vector<size_t> ind2sub(size_t index, std::vector<size_t> extents)
{
  std::vector<size_t> subscripts(extents.size());
  /* should check product of extents */
  size_t remainder = index;
  size_t total = 1;
  for (size_t i = 0; i < extents.size(); i++) {
    total *= extents[i];
  }
  if ((index >= total) || (index < 0)) {
    std::cerr << "ind2sub error: index is larger than product of"
      " dimensions given" << std::endl;
    return subscripts; // shouldn't really do this lol
  } else {
    for (size_t i = 1; i < subscripts.size(); i++) {
      size_t prod = 1;
      for (size_t j = 0; j < i; j++) {
        prod *= extents[j];
      }
      subscripts[i] = remainder / prod;
      remainder -= subscripts[i] * prod;
    }
    subscripts[0] = remainder;
  }
  return subscripts;
}

#endif
