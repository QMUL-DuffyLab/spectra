#include "vera.h"

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
        std::cerr << "subscripts = ";
        for (unsigned j = 0; j < subscripts.size(); j++) {
          std::cerr << subscripts[i];
        }
        std::cerr << std::endl;

        std::cerr << "extents = ";
        for (unsigned j = 0; j < subscripts.size(); j++) {
          std::cerr << extents[i];
        }
        std::cerr << std::endl;
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
  } else if ((int)index < 0) {
    std::cerr << "ind2sub error: index underflowed" << std::endl;
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

void
func(double t, double *y, double *ydot, void *data)
{
  vera_lsoda_data *vls = (vera_lsoda_data *)data;
  VERA *chromo = vls->chromo;
  std::vector<double> ydot_vec = chromo->dndt(y, t, (*vls->pump));
  size_t i;
  for (i = 0; i < ydot_vec.size(); i++) {
    ydot[i] = ydot_vec[i];
  }
}
 
/*
 * input file stuff
 */

VERA
create_VERA_from_file(char *filename)
{
  size_t i, n_elec, n_normal, n_vib = 0;
  double beta = 0.; double mu_ratio = 0.; double s2_stokes = 0.;
  std::string line, rest;
  std::string::size_type sz;
  std::vector<double> w_elec, w_normal, widths,
    l_ic_ij, l_ivr_i, g_ic_ij, g_ivr_i, disp;
  std::ifstream param_file;
  param_file.open(filename);
  if (param_file.is_open()) {
    while(param_file.good()) {
      std::getline(param_file, line);
      if (line.find("n_elec")!=std::string::npos) {
        sz = line.find("=") + 1; /* start after the equals! */
        n_elec = std::stoul(line.substr(sz));
      }
      if (line.find("n_normal")!=std::string::npos) {
        sz = line.find("=") + 1;
        n_normal = std::stoul(line.substr(sz));
      }
      if (line.find("n_vib")!=std::string::npos) {
        sz = line.find("=") + 1;
        n_vib = std::stoul(line.substr(sz));
      }
      if (line.find("T")!=std::string::npos) {
        sz = line.find("=") + 1;
        double t = std::stod(line.substr(sz));
        /* i want to be able to specify temperature in Kelvin
         * for convenience but the VERA functions use beta in
         * wavenumbers because that's how I wrote them. just
         * do the conversion here and be done with it :) */
        t *= KB / (2 * PI * HBAR * 100 * CVAC);
        beta = 1. / t;
      }
      if (line.find("mu_ratio")!=std::string::npos) {
        sz = line.find("=") + 1;
        mu_ratio = std::stod(line.substr(sz));
      }
      if (line.find("s2_stokes")!=std::string::npos) {
        sz = line.find("=") + 1;
        s2_stokes = std::stod(line.substr(sz));
      }
      /* arrays - given as name = space-separated list e.g. 
       * w_elec = 0.0 10000.0 20000.0 30000.0 */
      if (line.find("w_elec")!=std::string::npos) {
        sz = line.find("=") + 1;
        rest = line.substr(sz);
        for (i = 0; i < n_elec + 1; i++) {
          w_elec.push_back(std::stod(rest, &sz));
          rest = rest.substr(sz);
        }
      }
      if (line.find("w_normal")!=std::string::npos) {
        sz = line.find("=") + 1;
        rest = line.substr(sz);
        for (i = 0; i < n_normal; i++) {
          w_normal.push_back(std::stod(rest, &sz));
          rest = rest.substr(sz);
        }
      }
      if (line.find("widths")!=std::string::npos) {
        sz = line.find("=") + 1;
        rest = line.substr(sz);
        while (!rest.empty()) {
          widths.push_back(std::stod(rest, &sz));
          rest = rest.substr(sz);
        }
      }
      if (line.find("l_ic_ij")!=std::string::npos) {
        sz = line.find("=") + 1;
        rest = line.substr(sz);
        for (i = 0; i < (pow(n_elec, 2)); i++) {
          l_ic_ij.push_back(std::stod(rest, &sz));
          rest = rest.substr(sz);
        }
      }
      if (line.find("l_ivr_i")!=std::string::npos) {
        sz = line.find("=") + 1;
        rest = line.substr(sz);
        for (i = 0; i < n_elec; i++) {
          l_ivr_i.push_back(std::stod(rest, &sz));
          rest = rest.substr(sz);
        }
      }
      if (line.find("g_ic_ij")!=std::string::npos) {
        sz = line.find("=") + 1;
        rest = line.substr(sz);
        for (i = 0; i < (pow(n_elec, 2)); i++) {
          g_ic_ij.push_back(std::stod(rest, &sz));
          rest = rest.substr(sz);
        }
      }
      if (line.find("g_ivr_i")!=std::string::npos) {
        sz = line.find("=") + 1;
        rest = line.substr(sz);
        for (i = 0; i < n_elec; i++) {
          g_ivr_i.push_back(std::stod(rest, &sz));
          rest = rest.substr(sz);
        }
      }
      if (line.find("disp")!=std::string::npos) {
        sz = line.find("=") + 1;
        rest = line.substr(sz);
        for (i = 0; i < (n_normal * pow(n_elec + 1, 2)); i++) {
          disp.push_back(std::stod(rest, &sz));
          rest = rest.substr(sz);
        }
      }

    }
  }
  for (i = 0; i < g_ic_ij.size(); i++) {
    if (g_ic_ij[i] != 0.) {
      g_ic_ij[i] = 1. / (g_ic_ij[i] * 1E-3 * CM_PER_PS);
    }
  }
  for (i = 0; i < g_ivr_i.size(); i++) {
    if (g_ivr_i[i] != 0.) {
      g_ivr_i[i] = 1. / (g_ivr_i[i] * 1E-3 * CM_PER_PS);
    }
  }

  /* could add pumping parameters here if necessary */
  std::cout << "n_elec = "    << n_elec     << std::endl;
  std::cout << "n_normal = "  << n_normal   << std::endl;
  std::cout << "n_vib = "     << n_vib      << std::endl;
  std::cout << "beta = "      << beta       << std::endl;
  std::cout << "mu_ratio = "  << mu_ratio   << std::endl;
  std::cout << "s2_stokes = " << s2_stokes  << std::endl;
  std::cout << "w_elec = ";
  for (i = 0; i < w_elec.size(); i++ ) {
    std::cout << w_elec[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "w_normal = ";
  for (i = 0; i < w_normal.size(); i++ ) {
    std::cout << w_normal[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "widths = ";
  for (i = 0; i < widths.size(); i++ ) {
    std::cout << widths[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "l_ic_ij = ";
  for (i = 0; i < l_ic_ij.size(); i++ ) {
    std::cout << l_ic_ij[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "l_ivr_i = ";
  for (i = 0; i < l_ivr_i.size(); i++ ) {
    std::cout << l_ivr_i[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "g_ic_ij = ";
  for (i = 0; i < g_ic_ij.size(); i++ ) {
    std::cout << g_ic_ij[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "g_ivr_i = ";
  for (i = 0; i < g_ivr_i.size(); i++ ) {
    std::cout << g_ivr_i[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "disp =            ";
  for (i = 0; i < disp.size(); i++ ) {
    std::cout << disp[i] << " ";
  }
  std::cout << std::endl;
  VERA result = VERA(n_elec, n_normal, n_vib, beta,
      mu_ratio, s2_stokes, w_elec, w_normal, widths,
      l_ic_ij, l_ivr_i, g_ic_ij, g_ivr_i, disp);
  return result;
}

VERA::VERA(size_t elec, size_t norm, size_t vib,
    double bet, double mu_r, double s2_s,
    double* wi, size_t n_wi,
    double* w_mode, size_t n_mode,
    double* wid, size_t n_wid,
    double** l_ic, size_t side_l,
    double* l_ivr, size_t n_l,
    double** g_ic, size_t side_g,
    double* g_ivr, size_t n_g,
    double*** di, size_t n_mode_di, size_t e_l, size_t e_u
    )
  : n_elec(elec),
    n_normal(norm),
    n_vib(vib),
    beta(bet),
    mu_ratio(mu_r),
    s2_stokes(s2_s)
{
  set_extents();
  set_w_elec(wi, n_wi);
  set_w_normal(w_mode, n_mode);
  set_widths(wid, n_wid);
  set_l_ic_ij(l_ic, side_l);
  set_l_ivr_i(l_ivr, n_l);
  set_g_ic_ij(g_ic, side_g);
  set_g_ivr_i(g_ivr, n_g);
  set_disp(di, n_mode_di, e_l, e_u);
  set_k_ivr(beta);
  set_k_ic(beta);
  fc_calc();
}

VERA::VERA(size_t elec, size_t norm, size_t vib,
    double bet, double mu_r, double s2_s,
    std::vector<double> wi,
    std::vector<double> w_mode,
    std::vector<double> wid,
    std::vector<double> l_ic,
    std::vector<double> l_ivr,
    std::vector<double> g_ic,
    std::vector<double> g_ivr,
    std::vector<double> di
    )
  : n_elec(elec),
    n_normal(norm),
    n_vib(vib),
    beta(bet),
    mu_ratio(mu_r),
    s2_stokes(s2_s)
{
  set_extents();
  set_w_elec(wi);
  set_w_normal(w_mode);
  set_widths(wid);
  set_l_ic_ij(l_ic);
  set_l_ivr_i(l_ivr);
  set_g_ic_ij(g_ic);
  set_g_ivr_i(g_ivr);
  set_disp(di);
  set_k_ivr(beta);
  set_k_ic(beta);
  fc_calc();
}

void
VERA::set_w_elec(double* wi, size_t n_wi)
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
VERA::set_w_elec(std::vector<double> wi)
{
  /* note - could specify the wi as a vector too. HMMM */
  if (wi.size() != n_elec + 1) {
    std::cerr << "length of wi array does not match n_elec + 1"
              << std::endl;
  }
  for (size_t i = 0; i < wi.size(); i++) {
    w_elec.push_back(wi[i]);
  }
}

void
VERA::set_w_normal(double* w_mode, size_t n_mode)
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
VERA::set_w_normal(std::vector<double> w_mode)
{
  /* note - could specify the wi as a vector too. HMMM */
  if (w_mode.size() != n_normal) {
    std::cerr << "length of w_mode array does not match n_normal"
              << std::endl;
  }
  for (size_t i = 0; i < w_mode.size(); i++) {
    w_normal.push_back(w_mode[i]);
  }
}

void
VERA::set_widths(double* wid, size_t n_widths)
{
  /* note - could specify the wi as a vector too. HMMM */
  for (size_t i = 0; i < n_widths; i++) {
    widths.push_back(wid[i]);
  }
}

void
VERA::set_widths(std::vector<double> wid)
{
  /* note - could specify the wi as a vector too. HMMM */
  for (size_t i = 0; i < wid.size(); i++) {
    widths.push_back(wid[i]);
  }
}

void
VERA::set_l_ic_ij(double** l_ic, size_t side)
{
  if (side != n_elec) {
    std::cerr << "side of l_ic_ij array does not match n_elec"
              << std::endl;
  }
  for (size_t i = 0; i < side; i++) {
    for (size_t j = 0; j < side; j++) {
      l_ic_ij[sub2ind({i, j}, {side, side})] = l_ic[i][j];
    }
  }
}

void
VERA::set_l_ic_ij(std::vector<double> l_ic)
{
  if (l_ic.size() != pow(n_elec, 2)) {
    std::cerr << "side of l_ic_ij array does not match n_elec"
              << std::endl;
  }
  for (size_t i = 0; i < l_ic.size(); i++) {
    l_ic_ij[i] = l_ic[i];
  }
}

void
VERA::set_l_ivr_i(double* l_ivr, size_t n)
{
  /* note - could specify the wi as a vector too. HMMM */
  if (n != n_elec) {
    std::cerr << "length of l_ivr array does not match n_elec"
              << std::endl;
  }
  for (size_t i = 0; i < n; i++) {
    l_ivr_i[i] = l_ivr[i];
  }
}

void
VERA::set_l_ivr_i(std::vector<double> l_ivr)
{
  /* note - could specify the wi as a vector too. HMMM */
  if (l_ivr.size() != n_elec) {
    std::cerr << "length of l_ivr array does not match n_elec"
              << std::endl;
  }
  for (size_t i = 0; i < l_ivr.size(); i++) {
    l_ivr_i[i] = l_ivr[i];
  }
}

void
VERA::set_g_ic_ij(double** g_ic, size_t side)
{
  if (side != n_elec) {
    std::cerr << "side of g_ic_ij array does not match n_elec"
              << std::endl;
  }
  for (size_t i = 0; i < side; i++) {
    for (size_t j = 0; j < side; j++) {
      g_ic_ij[sub2ind({i, j}, {side, side})] = g_ic[i][j];
    }
  }
}

void
VERA::set_g_ic_ij(std::vector<double> g_ic)
{
  if (g_ic.size() != pow(n_elec, 2)) {
    std::cerr << "side of g_ic_ij array does not match n_elec"
              << std::endl;
  }
  for (size_t i = 0; i < g_ic.size(); i++) {
    g_ic_ij[i] = g_ic[i];
  }
}

void
VERA::set_g_ivr_i(double* g_ivr, size_t n)
{
  /* note - could specify the wi as a vector too. HMMM */
  if (n != n_elec) {
    std::cerr << "length of l_ivr array does not match n_elec"
              << std::endl;
  }
  for (size_t i = 0; i < n; i++) {
    g_ivr_i[i] = g_ivr[i];
  }
}

void
VERA::set_g_ivr_i(std::vector<double> g_ivr)
{
  /* note - could specify the wi as a vector too. HMMM */
  if (g_ivr.size() != n_elec) {
    std::cerr << "length of l_ivr array does not match n_elec"
              << std::endl;
  }
  for (size_t i = 0; i < g_ivr.size(); i++) {
    g_ivr_i[i] = (g_ivr[i]);
  }
}

void 
VERA::set_disp(double*** di,
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
        disp[sub2ind({i, j, k}, disp_extents)] = di[i][j][k];
      }
    }
  }
}

void 
VERA::set_disp(std::vector<double> di)
{
  if (di.size() != (n_normal * pow(n_elec + 1, 2))) {
    std::cerr << "size of disp vector doesn't match extents"
              << std::endl;
  }
  for (size_t i = 0; i < di.size(); i++) {
    disp[i] = di[i];
  }
}

void
VERA::fc_calc()
{
  size_t index;
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
              + (sqrt((double)a)/sqrt((double)b))
              * fc_array[a - 1][b - 1];
              fc_array[b][a] = (double)(pow(-1., (double)(b - a))) 
                             * fc_array[a][b];
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
}

void
VERA::set_extents()
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
  l_ic_ij.resize(size, 0.);
  g_ic_ij.resize(size, 0.);
  size = 1;

  ivr_extents.push_back(n_elec);
  ivr_extents.push_back(n_normal);
  ivr_extents.push_back(2);
  size = n_elec * n_normal * 2;
  k_ivr.resize(size, 0.);
  l_ivr_i.resize(size, 0.);
  g_ivr_i.resize(size, 0.);

  /* pretty sure this is always true, regardless
   * of the values of n_normal or n_elec */
  disp_extents = {n_normal, n_elec + 1, n_elec + 1};
  disp.resize((n_normal * pow(n_elec + 1, 2)), 0.);

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
VERA::get_w_elec(size_t elec_i)
{
  if (elec_i > n_elec + 1) {
    std::cout << "invalid index passed to get_w_elec: "
    << elec_i << " - max is " << n_elec + 1 << std::endl;
    return NAN;
  }
  return w_elec[elec_i];
}

std::vector<double>
VERA::get_w_normal()
{
  return w_normal;
}

double
VERA::get_w_normal(size_t norm)
{
  if (norm >= n_normal) {
    std::cout << "invalid index passed to get_w_normal: "
    << norm << " - max is " << n_normal << std::endl;
    return NAN;
  }
  return w_normal[norm];
}

double
VERA::get_l_ic_ij(size_t elec_i, size_t elec_j)
{
  return l_ic_ij[elec_i + (n_elec * elec_j)];
}

double
VERA::get_l_ivr_i(size_t elec_i)
{
  return l_ivr_i[elec_i];
}

double
VERA::get_g_ic_ij(size_t elec_i, size_t elec_j)
{
  return g_ic_ij[elec_i + (n_elec * elec_j)];
}

double
VERA::get_g_ivr_i(size_t elec_i)
{
  return g_ivr_i[elec_i];
}

double
VERA::get_disp(std::vector<size_t> subscripts)
{
  size_t index = sub2ind(subscripts, {n_normal, n_elec + 1, n_elec + 1});
  return disp[index];
}

double
VERA::get_fc(std::vector<size_t> subscripts)
{
  return fc[sub2ind(subscripts, fc_extents)];
}

std::vector<size_t>
VERA::get_extents()
{
  return extents;
}

std::vector<size_t>
VERA::get_pop_extents()
{
  return population_extents;
}

void
VERA::set_k_ivr(double beta)
{
  for (size_t i = 0; i < n_elec; i++) {
    for (size_t alpha = 0; alpha < n_normal; alpha++) {
      /* NB: python code puts these in backwards - the downward rate
       * first then the forward one, or so i thought.
       * this is the same ordering but check it!!
       */
      size_t index = sub2ind({i, alpha, 0}, ivr_extents);
      k_ivr[index] = CM_PER_PS *
          k_calc(w_normal[alpha], beta, l_ivr_i[i], g_ivr_i[i]);
      index = sub2ind({i, alpha, 1}, ivr_extents);
      k_ivr[index] = CM_PER_PS *
          k_calc(-1. * w_normal[alpha], beta, l_ivr_i[i], g_ivr_i[i]);
      /* k_ivr.push_back(CM_PER_PS * */
      /*     k_calc(w_normal[alpha], beta, l_ivr_i[i], g_ivr_i[i])); */
      /* k_ivr.push_back(CM_PER_PS * */
      /*     k_calc(-1. * w_normal[alpha], beta, l_ivr_i[i], g_ivr_i[i])); */
    }
  }
}

void
VERA::set_k_ic(double beta)
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
VERA::dndt(double *population, double t, pulse pump)
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
  for (size_t j = 0; j < (pow(n_vib + 1, n_normal)); j++) {
    std::vector<size_t> a = ind2sub(j, vib_extents);
    for (size_t k = 0; k < (pow(n_vib + 1, n_normal)); k++) {
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
    dndt[i] = dndt_ivr[i] + dndt_ic[i];
    /* dndt[i] = dndt_ivr[i] + dndt_ic[i] + dndt_pump[i]; */
    /* fprintf(stdout, "%5lu  %18.10f  %18.10f  %18.10f  %18.10f\n", */
    /*     i, dndt_ivr[i], dndt_ic[i], dndt_pump[i], dndt[i]); */ 
  }
  /* fprintf(stdout, "%18.10f  %18.10f  %18.10f\n", ivr_sum, ic_sum, pump_sum); */ 
  return dndt;
}

std::vector<double>
k_i_xa(VERA x, unsigned n_chl, unsigned n_car,
               unsigned tau, double **eig, double *eigvals,
               double **Jij, double **normed_ai, double **normed_fi,
               pulse v_abs, double beta)
{
  size_t n_total = x.get_pop_extents()[0];
  double ji = 0., ji_work = 0.;
  double *abs = (double *)calloc(tau, sizeof(double));
  double *fi_ad  = (double *)calloc(tau, sizeof(double));
  double *ai_fd  = (double *)calloc(tau, sizeof(double));
  std::vector<double> k_i_xa;
  unsigned short print_ji = 0;
  unsigned short print_delta_fc = 0;

  for (unsigned ii = 1; ii < x.get_pop_extents().size(); ii++) {
    n_total *= x.get_pop_extents()[ii];
  }

  /* note - the 48 * 48 rate things are the same for every chlorophyll
   * and in our case the same for both carotenoids, since we use the
   * same parameters for both. so we could define a 48 * 48 array
   * of those and look them up; would probably be faster overall. */

  for (unsigned chl_index = 0; chl_index < n_chl; chl_index++) {
    for (unsigned carotenoid = n_chl;
        carotenoid < n_chl + n_car; carotenoid++) {
      /* assumes the carotenoids are at the end of the arrays;
       * fortran code puts them there so shouldn't be an issue */

      ji = 0.;
      for (unsigned n = 0; n < n_chl; n++) {
        for (unsigned m = 0; m < n_chl; m++) {
          ji += eig[m][chl_index] * eig[n][chl_index]
             * Jij[m][carotenoid] * Jij[n][carotenoid];
        }
      }

      if (print_ji) {
        fprintf(stdout, "%5u %12.8e\n", chl_index, ji);
      }

      for (unsigned ii = 0; ii < n_total; ii++) {
        std::vector<size_t> a = ind2sub(ii, x.get_pop_extents());
        double e_xa = 0.;
        double chl_car = 0.;
        double car_chl = 0.;

        if (a[0] == 0) {
          k_i_xa.push_back(0.);
          k_i_xa.push_back(0.);
          continue; /* excitons cannot couple to g/s */
        } else {
          for (unsigned kk = 0; kk < n_total; kk++) {
            if (ii == kk) continue;

            std::vector<size_t> b = ind2sub(kk, x.get_pop_extents());
            if (a[0] == b[0]) {
              continue;
            } else {
              double delta_xy_ba = x.get_w_elec(a[0]) - x.get_w_elec(b[0]);
              e_xa = x.get_w_elec(a[0]);
              /* NB: check with chris that this is 1200 and not like 800 */
              /* FWHM should be 1150 or 1200 */
              v_abs.width = 1200.0;
              double fc_sq = 1.;

              for (unsigned alpha = 0; alpha < x.n_normal; alpha++) {
                delta_xy_ba += x.get_w_normal(alpha)
                             * (int)(a[alpha + 1] - b[alpha + 1]);
                e_xa += x.get_w_normal(alpha) * a[alpha + 1];
                fc_sq *= pow(x.get_fc({a[0], b[0], alpha,
                      b[alpha + 1], a[alpha + 1]}), 2.);
              }

              /* now we have to calculate the rate between every delta_xy_ba */
              v_abs.centre = delta_xy_ba;
              if (v_abs.width != 0 && delta_xy_ba != 0.) {
                abs = incident(v_abs, tau);

                /* this is probably not right - placeholder!!
                 * need to know what the reorganisation is for
                 * each transition and then do the reverse rate */
                /* v_abs.centre += 2.* x.get_l_ic_ij(a[0], b[0]); */
                /* flu = incident(v_abs, tau); */
              }

              ji_work = ji * fc_sq;
              /* the coupling is actually to S1!! */
              /* i've done this this way because it's possible that
               * we might want couplings to other electronic states
               * as well; in that case, we'd keep the sums like this.
               * if we never want that we could limit the ii sum to
               * a[0] == 1 only */
              if (a[0] != 1) {
                ji_work = 0.;
              }

              for (unsigned step = 0; step < tau; step++) {
                fi_ad[step] = normed_fi[chl_index][step] * abs[step];
                /* abs in the next line was flu!! */
                ai_fd[step] = normed_ai[chl_index][step] * abs[step];
              }
              
              if (print_delta_fc) {
                fprintf(stdout, "%5u %5u %12.8e %12.8e %12.8e\n", ii, kk, 
                    delta_xy_ba, fc_sq, ji_work);
              }

              chl_car += pow(ji_work, 2.) * trapezoid(fi_ad, tau);
              car_chl += pow(ji_work, 2.) * trapezoid(ai_fd, tau);

              /* ki_delta_xy_ba.push_back(CM_PER_PS * (2 * PI) * */
              /*     pow(ji_work, 2.) * trapezoid(fi_ad, tau)); */
              /* ki_delta_xy_ba.push_back(CM_PER_PS * (2 * PI) * */
              /*     pow(ji_work, 2.) * trapezoid(ai_fd, tau)); */
            /* } */
            }
          } // kk
        }

        /* here we can do the detailed balance check - the push back
         * will be here once the sum is sorted out */
        chl_car *= CM_PER_PS * 2 * PI;
        car_chl *= CM_PER_PS * 2 * PI;
        /* penalise the upward rate - detailed balance */
        if (eigvals[chl_index] > e_xa) {
          car_chl *= exp(-beta * (eigvals[chl_index] - e_xa));
        } else {
          chl_car *= exp(-beta * (e_xa - eigvals[chl_index]));
        }
        k_i_xa.push_back(chl_car);
        k_i_xa.push_back(car_chl);

        bool print_details = false;
        if (print_details) {
          fprintf(stdout, "chl = %2u, car = %2u,"
          " state = (%2lu, %2lu, %2lu): e_chl = %10.6e,"
          " e_car = %10.6e, forward rate = %10.6e,"
          " backward_rate = %10.6e\n", chl_index, carotenoid,
          a[0], a[1], a[2], eigvals[chl_index], e_xa, chl_car, car_chl);
        }

      } // ii
    }
  } // chl_index
  free(abs);
  free(fi_ad);
  free(ai_fd);

  return k_i_xa;
}

std::vector<double>
VERA::intra_rates()
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
  
  std::vector<double> rates(pow(n_total, 2), 0.);

  std::vector<size_t> subscripts(n_extents.size(), 0);
  std::vector<size_t> sub_lower(n_extents.size(), 0),
                      sub_upper(n_extents.size(), 0);
  size_t i_lower = 0, i_upper = 0;

  bool print_ivr = false;
  if (print_ivr) {
    fprintf(stdout, "IVR rates:\n");
    for (size_t i = 0; i < k_ivr.size(); i++) {
      fprintf(stdout, "%2lu %10.6e\n", i, k_ivr[i]);
    }
  }


  for (size_t i = 0; i < n_total; i++) {
    subscripts = ind2sub(i, n_extents);  
     
    /* IVR. can probably be optimised */
    for (size_t alpha = 0; alpha < n_normal; alpha++) {
      if (subscripts[alpha + 1] > 0) {

        sub_lower = subscripts;
        sub_lower[alpha + 1]--;
        i_lower = sub2ind(sub_lower, n_extents);

        /* loss of population to lower vibronic state */
        rates[sub2ind({i, i}, {n_total, n_total})] -=
          subscripts[alpha + 1]
        * k_ivr[sub2ind({subscripts[0], alpha, 0},
          {n_elec, n_normal, 2})];

        /* gain of population from level below */
        rates[sub2ind({i, i_lower}, {n_total, n_total})] +=
          subscripts[alpha + 1]
        * k_ivr[sub2ind({subscripts[0], alpha, 1},
          {n_elec, n_normal, 2})];
      }

      if (subscripts[alpha + 1] < n_vib) {
        sub_upper = subscripts;
        sub_upper[alpha + 1]++;
        i_upper = sub2ind(sub_upper, n_extents);

        /* loss of population to upper state */
        rates[sub2ind({i, i}, {n_total, n_total})] -=
          (subscripts[alpha + 1] + 1.)
        * k_ivr[sub2ind({subscripts[0], alpha, 1},
          {n_elec, n_normal, 2})];
        /* gain from the upper state */
        rates[sub2ind({i, i_upper}, {n_total, n_total})] +=
          (subscripts[alpha + 1] + 1.)
        * k_ivr[sub2ind({subscripts[0], alpha, 0},
          {n_elec, n_normal, 2})];
      }
    }

  }

  if (print_ivr) {
    fprintf(stdout, "IVR total rates:\n");
    for (size_t i = 0; i < n_total; i++) {
      for (size_t j = 0; j < n_total; j++) {
        fprintf(stdout, "%2lu %2lu %10.6e\n", i, j,
            rates[sub2ind({i, j}, {n_total, n_total})]);
      }
    }
  }

  /* interconversion */
  /* i make the IC vector all zeroes and then just fill in
   * the non-zero ones, so we can just list them */

  bool print_ic = false;
  if (print_ic) {
    fprintf(stdout, "IC rates:\n");
  }
  for (size_t i = 0; i < n_elec; i++) {
    for (size_t j = 0; j < n_elec; j++) {
      if (i != j) {
        for (size_t k = 0; k < (pow((n_vib + 1), n_normal)); k++) {
          std::vector<size_t> a = ind2sub(k, vib_extents);
          for (size_t l = 0; l < (pow((n_vib + 1), n_normal)); l++) {
            std::vector<size_t> b = ind2sub(l, vib_extents);

            /* indices */
            std::vector<size_t> n_i, n_j;
            std::vector<size_t> k_ba_ji, k_ba_ij, k_ab_ji, k_ab_ij;

            n_i.push_back(i);
            n_j.push_back(j);
            k_ba_ji = {j, i};
            k_ba_ij = {i, j};
            k_ab_ji = {j, i};
            k_ab_ij = {i, j};
            for (size_t ai = 0; ai < a.size(); ai++) {
              n_i.push_back(a[ai]);
              k_ab_ij.push_back(a[ai]);
              k_ab_ji.push_back(a[ai]);
            }
            for (size_t bi = 0; bi < b.size(); bi++) {
              n_j.push_back(b[bi]);
              k_ba_ji.push_back(b[bi]);
              k_ba_ij.push_back(b[bi]);
            }

            /* two sets of vibrational indices to worry about */
            for (size_t ai = 0; ai < a.size(); ai++) {
              k_ba_ji.push_back(a[ai]);
              k_ba_ij.push_back(a[ai]);
            }
            for (size_t bi = 0; bi < b.size(); bi++) {
              k_ab_ij.push_back(b[bi]);
              k_ab_ji.push_back(b[bi]);
            }

            size_t n_i_ind = sub2ind(n_i, n_extents);
            size_t n_j_ind = sub2ind(n_j, n_extents);

            double fc_ij = 1., fc_ji = 1.;
            for (size_t alpha = 0; alpha < n_normal; alpha++) {
              if (j > i) {
                fc_ij *= pow(fc[sub2ind({i, j, alpha, 
                          a[alpha], b[alpha]}, fc_extents)], 2.);
              } else {
                fc_ij *= pow(fc[sub2ind({j, i, alpha, 
                          b[alpha], a[alpha]}, fc_extents)], 2.);
              }
            }

            if (print_ic) {
              fprintf(stdout, "%1lu %1lu %1lu <-> %1lu %1lu %1lu "
                  " %10.6e %10.6e\n",
                  n_i[0], n_i[1], n_i[2], 
                  n_j[0], n_j[1], n_j[2], 
                  fc_ij * k_ic[sub2ind(k_ba_ji, ic_extents)],
                  fc_ij * k_ic[sub2ind(k_ab_ij, ic_extents)]
                  );
            }

            if (j > i) {
              rates[sub2ind({n_i_ind, n_i_ind}, {n_total, n_total})] -=
                fc_ij * k_ic[sub2ind(k_ba_ji, ic_extents)];
              rates[sub2ind({n_i_ind, n_j_ind}, {n_total, n_total})] +=
                fc_ij * k_ic[sub2ind(k_ba_ij, ic_extents)];
            } else {
              rates[sub2ind({n_i_ind, n_i_ind}, {n_total, n_total})] -=
                fc_ij * k_ic[sub2ind(k_ab_ji, ic_extents)];
              rates[sub2ind({n_i_ind, n_j_ind}, {n_total, n_total})] +=
                fc_ij * k_ic[sub2ind(k_ab_ij, ic_extents)];
            }

          } /* end of electronic state if/else */

        }
      } /* i != j */
    }
  }

  return rates;
}

void
intra_rate_test(double *population, std::vector<double> rates,
size_t n_total, double *res)
{
  double **rate_matrix = (double **)calloc(n_total, sizeof(double*));
  for (unsigned i = 0; i < n_total; i++) {
    rate_matrix[i] = (double *)calloc(n_total, sizeof(double));
  }
  for (unsigned i = 0; i < n_total; i++) {
    for (unsigned j = 0; j < n_total; j++) {
      rate_matrix[i][j] = rates[sub2ind({i, j}, {n_total, n_total})];
    }
  }

  matvec(n_total, rate_matrix, population, res);

  for (unsigned i = 0; i < n_total; i++) {
    free(rate_matrix[i]);
  }
  free(rate_matrix);
}

double **total_rates(unsigned n_chl, VERA car, unsigned n_car,
unsigned n_s_car, double *gamma, double **Jij, std::vector<double> k_i_delta,
double **redfield_rates)
{
  unsigned n_tot = n_chl + 1 + (n_car * n_s_car);
  double **k_tot = (double **)calloc(n_tot, sizeof(double*));
  for (unsigned i = 0; i < n_tot; i++) {
    k_tot[i] = (double *)calloc(n_tot, sizeof(double));
  }

  std::vector<size_t> extents = {n_chl, n_car, n_s_car, 2};
  std::vector<double> intra_car = car.intra_rates();

  /* note - might not need 0 < j < n_tot; might just be
   * i < j < n_tot and then fill in the rest with detailed balance etc.
   * but i need to think about that - redfield rates have it built in,
   * i think we have to do it manually here when summing over the
   * chl-car rates? */

  /* also: the intra-carotenoid rates and the redfield rates matrix
   * are already transfer matrices, so we just fill them in as
   * k_total[i][j] = blah blah blah. The chlorophyll-carotenoid rates
   * aren't, so we switch the indices */
  for (unsigned i = 0; i < n_tot; i++) {
    for (unsigned j = 0; j < n_tot; j++) {
      bool i_rgs = (i == 0);
      bool j_rgs = (j == 0);
      bool i_chls = (i > 0 && i <= n_chl);
      bool j_chls = (j > 0 && j <= n_chl);
      bool i_620 = (i > n_chl && i < n_chl + 1 + n_s_car);
      bool j_620 = (j > n_chl && j < n_chl + 1 + n_s_car);
      bool i_621 = (i >= n_chl + 1 + n_s_car);
      bool j_621 = (j >= n_chl + 1 + n_s_car);

      if (i_rgs) {
        if (j_rgs) {
          continue;
        }
        if (j_chls) {
          k_tot[i][j] += (1. / (1000 * gamma[j - 1]));
        }
        if (j_620) {
          continue;
        }
        if (j_621) {
          continue;
        }
      }

      if (j_rgs) {
        /* nothing comes out of Redfield ground state */
        continue;
      }

      if (i_chls) {
        if (j_rgs) {
          continue;
        }
        if (j_chls) {
          /* redfield rate: i/j - 1 because of ground state */
          /* [i][j] not [j][i] because this is already a transfer matrix */
          k_tot[i][j] = redfield_rates[i - 1][j - 1];
          if (i == j) {
            /* need to subtract the outward rates to the carotenoids */
            for (unsigned k = 0; k < n_s_car; k++) {
              k_tot[i][j] -= k_i_delta[sub2ind({i - 1, 0, k, 0}, extents)];
              k_tot[i][j] -= k_i_delta[sub2ind({i - 1, 1, k, 0}, extents)];
            }
          }
        }
        if (j_620) {
          k_tot[j][i] = k_i_delta[sub2ind({i - 1, 0, j - (n_chl + 1), 0},
                        extents)];
        }
        if (j_621) {
          k_tot[j][i] = k_i_delta[sub2ind({i - 1, 1,
              j - (n_s_car + n_chl + 1), 0}, extents)];
        }
      }

      if (i_620) {
        if (j_rgs) {
          continue;
        }
        if (j_chls) {
          k_tot[j][i] = k_i_delta[sub2ind({j - 1, 0, i - (n_chl + 1), 1},
                        extents)];
        }
        if (j_620) {
          k_tot[i][j] = intra_car[sub2ind({i - (n_chl + 1),
              j - (n_chl + 1)}, {n_s_car, n_s_car})];
          if (i == j) {
            /* need to subtract the outward rates to the chls */
            for (unsigned k = 0; k < n_chl; k++) {
              k_tot[i][j] -= k_i_delta[sub2ind({k, 0, j - (n_chl + 1), 1},
                            extents)];
            }
          }
        }
        if (j_621) {
          continue; /* inter-car rates assumed to be zero */
        }
      }

      if (i_621) {
        if (j_rgs) {
          continue;
        }
        if (j_chls) {
          k_tot[j][i] = k_i_delta[sub2ind({j - 1, 1,
              i - (n_s_car + n_chl + 1), 1}, extents)];
        }
        if (j_620) {
          continue;
        }
        if (j_621) {
          /* sum over IVR/IC rates */
          k_tot[i][j] = intra_car[sub2ind({i - (n_s_car + n_chl + 1),
              j - (n_s_car + n_chl + 1)}, {n_s_car, n_s_car})];
          if (i == j) {
            for (unsigned k = 0; k < n_chl; k++) {
              k_tot[i][j] -= k_i_delta[sub2ind({k, 1,
                  j - (n_s_car + n_chl + 1), 1}, extents)];
            }
          }
        }
      }

    }
  }
  return k_tot;
}

