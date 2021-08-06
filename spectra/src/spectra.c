#include <fftw3.h>
#include <stdio.h>
#include <cmath>
#include <complex>
#include <gsl/gsl_eigen.h>
#include "input.h"
#include "cd.h"
#include "steady_state.h"
#include "forster.h"
#include "vera.h"
#include "hybrid_vera.h"

int
main(int argc, char** argv)
{
  FILE *fp;
  unsigned int i, j;
  int status;
  char *line, **lineshape_files;
  double kd;
  fftw_complex **gi_array;
  double *eigvals, *gamma, *rates, *musq, *chi_p,
         *lambda, *integral, *chiw_ints, *ww,
         **wij, **kij, **Jij, **mu, **eig, **chiw, **pump;
  Parameters *line_params;
  fftw_complex *in, *ft_in, *out, *ft_out;
  fftw_plan plan, ft_plan;

  fprintf(stdout, "== Spectra calculation code ==\n");

  if (argc != 4) {
    fprintf(stdout, "Wrong number of arguments - should be 3: "
        "Input file, protocol file and list of lineshape defs\n");
    exit(EXIT_FAILURE);
  }

  Input *p = read_input_file(argv[1]);
  Protocol protocol = get_protocol(argv[2]);
  fprintf(stdout, "τ = %5d, T = %9.4f\n", p->tau, p->T);
  fprintf(stdout, "n_chl = %2d, n_car = %2d\n", p->n_chl, p->n_car);

  /* specifies form of incident light for source term in P_i eqns */
  pulse pump_properties = {
    .type=GAUSSIAN,
    .target_state = 2,
    .amplitude = 1.,
    .centre=15873.,
    .width=50.,
    .t_peak = 0.1,
    .duration = 70.0E-3
  };

  pulse VERA_absorption = {
    .type = GAUSSIAN,
    .target_state = 0,
    .amplitude = 0.,
    .centre = 0.,
    .width = 0.,
    .t_peak = 0.1,
    .duration = 70.0E-3,
  };

  bool calculate_CD = true;

  /* form of P_i(0) - check steady_state.h for details */
  ss_init population_guess = MUSQ;

  /* this is a kind of hack to get rid of some compiler warnings:
   * really, I should read in the filenames for the VERA parameters
   * as const char** and then iterate over that list to create
   * the vector of carotenoids here. but that's unnecessary rn */
  const char* lut1 = "in/LUT620.dat";
  const char* lut2 = "in/LUT621.dat";

  std::vector<VERA> cars;
  if (p->n_car == 1) {
    cars.push_back(create_VERA_from_file(lut1));
  } else {
    cars.push_back(create_VERA_from_file(lut1));
    cars.push_back(create_VERA_from_file(lut2));
  }

  /* malloc 1d stuff, read them in */
  eigvals     = (double *)calloc(p->n_chl, sizeof(double));
  gamma       = (double *)calloc(p->n_chl, sizeof(double));
  rates       = (double *)calloc(p->n_chl, sizeof(double));
  musq        = (double *)calloc(p->n_chl, sizeof(double));
  lambda      = (double *)calloc(p->n_chl, sizeof(double));
  line_params = (Parameters *)malloc(p->n_chl * sizeof(Parameters));
  chiw_ints   = (double *)calloc(p->n_chl, sizeof(double));
  ww          = (double *)calloc(p->tau, sizeof(double));
  integral    = (double *)calloc(p->tau, sizeof(double));
  in          = (fftw_complex *)fftw_malloc(p->tau * sizeof(fftw_complex));
  out         = (fftw_complex *)fftw_malloc(p->tau * sizeof(fftw_complex));
  ft_in       = (fftw_complex *)fftw_malloc(p->tau * sizeof(fftw_complex));
  ft_out      = (fftw_complex *)fftw_malloc(p->tau * sizeof(fftw_complex));

  /* malloc 2d stuff */
  lineshape_files = (char **)malloc(p->n_chl * sizeof(char*));
  mu              = (double **)calloc(p->n_chl, sizeof(double*));
  gi_array        = (fftw_complex**)
                    fftw_malloc(p->n_chl * sizeof(fftw_complex*));
  wij             = (double **)calloc(p->n_chl, sizeof(double*));
  kij             = (double **)calloc(p->n_chl, sizeof(double*));
  chiw            = (double **)calloc(p->n_chl, sizeof(double*));
  pump            = (double **)calloc(p->n_chl, sizeof(double*));

  eig             = (double **)calloc(p->N, sizeof(double*));
  Jij             = (double **)calloc(p->N, sizeof(double*));
  for (i = 0; i < p->N; i++) {
    Jij[i]             = (double *)calloc(p->N, sizeof(double));
    eig[i]             = (double *)calloc(p->N, sizeof(double));
  }
                      
  for (i = 0; i < p->n_chl; i++) {
    lineshape_files[i] = (char *)malloc(200 * sizeof(char));
    gi_array[i]        = (fftw_complex*)
                          fftw_malloc(p->tau * sizeof(fftw_complex));
    mu[i]              = (double *)calloc(3, sizeof(double));
    wij[i]             = (double *)calloc(p->n_chl, sizeof(double));
    kij[i]             = (double *)calloc(p->n_chl, sizeof(double));
    chiw[i]            = (double *)calloc(p->tau, sizeof(double));
    pump[i]            = (double *)calloc(p->tau, sizeof(double));
  }

  /* read function reads a set of doubles: here we read in
   * the rates (gamma), reorganisation energies (lambda) 
   * and exciton energies (eigvals - they're the eigenvalues
   * of the diagonalised Hamiltonian */
  gamma   = read(p->gamma_file, p->N);
  lambda  = read(p->lambda_file, p->N);
  eigvals = read(p->eigvals_file, p->N);

  /* read 2d stuff */
  gi_array = read_gi(p->gi_files, p->n_chl, p->tau);
  eig      = read_eigvecs(p->eigvecs_file, p->N);
  mu       = read_mu(p->mu_file, p->n_chl);
  ww       = incident(pump_properties, p->tau);

  Jij      = read_eigvecs(p->jij_file, p->N);

  /* fprintf(stdout, "EIG CHECK!\n"); */
  /* for (unsigned i = 0; i < p->n_chl; i++) { */
  /*   for (unsigned j = 0; j < p->n_chl; j++) { */
  /*     fprintf(stdout, "eig[%2u][%2u] = %8.5f\n", i, j, eig[i][j]); */
  /*   } */
  /* } */

  fp = fopen(argv[3], "r"); /* read in list of lineshape files here */
  line = (char *)malloc(200 * sizeof(char));
  for (i = 0; i < p->n_chl; i++) {
    /* now load in the parameters for each ligand */
    fgets(line, 200, fp);
    line[strcspn(line, "\n")] = 0;
    strcpy(lineshape_files[i], line);
    line_params[i] = get_parameters(lineshape_files[i],
                     protocol.chl_ansatz);
    if (line_params[i].ans == 3) {
      /* this'll have to be changed somehow */
      /* vera = create_VERA_from_file("in/vera_params.dat"); */
    }
    line_params[i].cw = choose_ansatz(line_params[i].ans);
    line_params[i].cn = &c_n;
    /* this is also kinda ugly - the temperature's a global parameter
     * really, not a lineshape one - but it avoids setting the
     * temperature in every lineshape file, at least */
    line_params[i].T = protocol.T;

    for (j = 0; j < p->n_chl; j++) {
      /* wij is \omega_{ij} - the gap between the 0-0 lines of two
       * excitons. It's used in calculating the Redfield rates later */
      wij[i][j] = ((eigvals[i] - lambda[i]) - (eigvals[j] - lambda[j]));
    }
  }
  int cl = fclose(fp);
  if (cl != 0) {
    fprintf(stdout, "Failed to close list of lineshape files %d.\n", cl);
    exit(EXIT_FAILURE);
  }

  kij = rate_calc(p->n_chl, eig, wij, line_params);

  /* fprintf(stdout, "\n\nNEXT EIG CHECK:\n\n"); */
  /* for (unsigned i = 0; i < p->n_chl; i++) { */
  /*   for (unsigned j = 0; j < p->n_chl; j++) { */
  /*     fprintf(stdout, "eig[%2u][%2u] = %8.5f\n", i, j, eig[i][j]); */
  /*   } */
  /* } */

  check_detailed_balance(p->n_chl, protocol.T, 1e-5, kij, wij);
  rates = relaxation_rates(p->n_chl, gamma);
  char fn[200];
  strcpy(fn, p->fw_file);
  status = generate_filename(sizeof(fn), fn, "fw", "rates");
  fp = fopen(fn, "w");
  print_matrix(fp, NULL, p->n_chl, kij);
  cl = fclose(fp);
  if (cl != 0) {
    fprintf(stdout, "Failed to close list of lineshape files %d.\n", cl);
    exit(EXIT_FAILURE);
  }


  plan = fftw_plan_dft_1d(p->tau, 
  	 in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  ft_plan = fftw_plan_dft_1d(p->tau, 
  	 ft_in, ft_out, FFTW_BACKWARD, FFTW_ESTIMATE);

  double **normed_ai = (double **)calloc(p->n_chl, sizeof(double*));
  double **normed_fi = (double **)calloc(p->n_chl, sizeof(double*));
  double ai_sum = 0., fi_sum = 0.;
  for (unsigned i = 0; i < p->n_chl; i++) {
    normed_ai[i] = (double *)calloc(p->tau, sizeof(double));
    normed_fi[i] = (double *)calloc(p->tau, sizeof(double));
  }

  COMPLEX *Atv = (COMPLEX *)calloc(p->tau, sizeof(COMPLEX));
  COMPLEX *Ftv = (COMPLEX *)calloc(p->tau, sizeof(COMPLEX));
  for (i = 0; i < p->n_chl; i++) {
    ai_sum = 0.; fi_sum = 0;
    musq[i] = pow(mu[i][0], 2.) + pow(mu[i][1], 2.) + pow(mu[i][2], 2.);

    for (unsigned int j = 0; j < p->tau; j++) {
      Atv[j] = At(eigvals[i], gi_array[i][j][0], gi_array[i][j][1],
                 (double)j * TOFS,
                 1. / (1000000. * gamma[i]));
      Ftv[j] = Ft(eigvals[i],
                  gi_array[i][j][0], gi_array[i][j][1],
                  lambda[i], (double)j * TOFS,
                  1. / (1000000. * gamma[i]));

      in[j][0] = REAL(Atv[j]);
      in[j][1] = IMAG(Atv[j]);
      ft_in[j][0] = REAL(Ftv[j]);
      ft_in[j][1] = IMAG(Ftv[j]);
    }

    fftw_execute(plan); 
    fftw_execute(ft_plan); 
    for (unsigned int j = 0; j < p->tau; j++) {
      chiw[i][j]   = out[j][0]  * musq[i];
      pump[i][j]   = chiw[i][j] * ww[j];
      integral[j] += out[j][0]  * musq[i];
      /* double omega = j * 2. * M_PI / (TOFS * p->tau); */
      /* normed_ai[i][j] = (out[j][0] * omega); */
      /* normed_fi[i][j] = (ft_out[j][0] * pow(omega, 3.)); */
      normed_ai[i][j] = out[j][0];
      normed_fi[i][j] = ft_out[j][0];
    }
    /* 2. * pi / TOFS * p->tau is the (constant) dx for the integral:
     * the 2. pi is due to the definition of TOFS!!! */
    /* ai_sum = trapezoid(normed_ai[i], 1. / (TOFS * p->tau), p->tau); */
    /* fi_sum = trapezoid(normed_fi[i], 1. / (TOFS * p->tau), p->tau); */
    ai_sum = trapezoid(normed_ai[i], 2. * M_PI / (TOFS * p->tau), p->tau);
    fi_sum = trapezoid(normed_fi[i], 2. * M_PI / (TOFS * p->tau), p->tau);
    for (unsigned int j = 0; j < p->tau; j++) {
      /* NB: the 2\pi is because FFTW outputs unnormalised transforms */
      /* normed_ai[i][j] = normed_ai[i][j] / (2. * M_PI * ai_sum); */
      /* normed_fi[i][j] = normed_fi[i][j] / (2. * M_PI * fi_sum); */
      normed_ai[i][j] = normed_ai[i][j] / (ai_sum);
      normed_fi[i][j] = normed_fi[i][j] / (fi_sum);
    }

    chi_p = pump[i];
    chiw_ints[i] = trapezoid(chi_p, 2. * M_PI / (TOFS * p->tau), p->tau);

  }
  free(Atv);
  free(ww);

  fprintf(stdout, "\nWriting A(w) file\n");
  fp = fopen(p->aw_file, "w");
  for (i = 0; i < p->tau; i++) {
    /* unpack the ordering used by FFTW */
    kd = i * 2. * M_PI / (p->tau);
    fprintf(fp, "%18.10f %18.10f\n", kd / TOFS, 
        kd * integral[i] * (1./ sqrt(p->tau)));
  }
  cl = fclose(fp);
  if (cl != 0) {
    fprintf(stdout, "Failed to close A(w) output file %d.\n", cl);
    exit(EXIT_FAILURE);
  }

  fprintf(stdout, "\nWriting χ_i(w) files\n");
  /* replace g_ with chi_ in gi filenames so we can print out the 
   * individual parts of the spectra. strings in C are so annoying.
   * NB: assumes 200 is enough bytes for gi filename + a few chars! */
  for (i = 0; i < p->n_chl; i++) {
    strcpy(fn, p->gi_files[i]);
    status = generate_filename(sizeof(fn), fn, "g_", "chi_");
    if(status != 0) {
      fprintf(stdout, "Cannot find string 'g_' in "
          "g_i filename no. %d, needed to print out chi_i. "
          "Something got renamed accidentally?\n", i);
      break;
    }
    fp = fopen(fn, "w");
    for (j = 0; j < p->tau; j++) {
      /* unpack the ordering used by FFTW */
      kd = j * 2. * M_PI / (p->tau);
      fprintf(fp, "%18.10f %18.10f\n", kd / TOFS, 
          kd * chiw[i][j] * (1./ sqrt(p->tau)));
    }
    cl = fclose(fp);
    if (cl != 0) {
      fprintf(stdout, "Failed to close χ_i(w) "
          "output file no. %d, error no. %d.\n", i, cl);
      exit(EXIT_FAILURE);
    }
  }

  /* one with the highest oscillator strength gets excited? */
  unsigned int max = 0;
  double musq_max = musq[0];
  double musq_sum = 0.;
  double chiw_sum = 0.;
  fprintf(stdout, "\n----------------------------------\n"
      "OSC. STRENGTHS AND CHI(W) INTEGRAL\n"
      "----------------------------------\n\n");
  fprintf(stdout, "Pigment        |μ^2|      ∫χ_i(w)\n");
  for (i = 0; i < p->n_chl; i++) {
    fprintf(stdout, "%7d %10.6e %10.6e\n", i + 601, musq[i], chiw_ints[i]);
    if (musq[i] > musq_max) {
      max = i;
      musq_max = musq[i];
      musq_sum += musq[i];
      chiw_sum += chiw_ints[i];
    }
  }

  ode_params odep;
  odep.N = p->n_chl;
  odep.kij = kij;
  odep.rates = rates;
  odep.chiw = chiw_ints;
  odep.Tij = transfer_matrix(p->n_chl, rates, kij);
  double *f = (double *)calloc(p->n_chl, sizeof(double));
  /* double *y = (double *)calloc(p->n_chl, sizeof(double)); */
  /* check convergence */
  double *boltz = (double *)calloc(p->n_chl, sizeof(double));

  boltz = bcs(p->n_chl, eigvals, protocol.T);
  fprintf(stdout, "\n-----------------\nBOLTZMANN WEIGHTS\n"
                  "-----------------\n\n");
  fprintf(stdout, "Pigment    p_i   |μ^2|*p_i\n");
  for (i = 0; i < p->n_chl; i++) {
    fprintf(stdout, "%7d %8.6f %8.6f\n", i + 601, boltz[i],
        boltz[i] * musq[i]);
    /* possible initial values for transient absorption? */
    /* if (i == max) { */
    /*   y[i] = 1.0; */
    /* } else { */
    /*   y[i] = 0.0; */
    /* } */
  }
  fprintf(stdout, "\n");

  fprintf(stdout, "\n----------\n"
                    "VERA RATES\n"
                    "----------\n\n");

  /* VERA vera = create_VERA_from_file("in/vera.def"); */

  size_t n_chl    = p->n_chl;
  size_t n_car    = p->n_car;

  bool hybrid = true;
  std::vector<size_t> n_s_cars(n_car, 0);
  size_t n_total = n_chl + 1;
  if (hybrid) {
    for (unsigned i = 0; i < n_car; i++) {
      n_s_cars[i]  = pow(cars[i].n_vib + 1, cars[i].n_normal);
      n_total     += (n_s_cars[i] + 1);
    }
  } else {
    for (unsigned i = 0; i < n_car; i++) {
      n_s_cars[i] = cars[i].n_total;
      n_total    += n_s_cars[i];
    }
  }

  double beta     = 1.439 / protocol.T;
  double **car_decays = (double **)calloc(n_car, sizeof(double));
  for (unsigned i = 0; i < n_car; i++) {
    car_decays[i] = (double *)calloc(n_s_cars[i], sizeof(double));
  }

  if (car_decays == NULL) {
    fprintf(stderr, "car_decays allocation failed!\n");
    exit(EXIT_FAILURE);
  }

  std::vector<std::vector<double>> k_chl_car;
  if (hybrid) {
    k_chl_car = k_i_xa_hybrid(cars, n_chl,
      n_car, p->tau, eig, eigvals,
      Jij, normed_ai, normed_fi,
      VERA_absorption, beta, car_decays);
  } else {
    /* k_chl_car = k_i_xa(vera, n_chl, */
    /*     n_car, p->tau, eig, eigvals, Jij, normed_ai, normed_fi, */
    /*     VERA_absorption, beta); */
  }
  double k_sum = 0.;

  bool print_decays = true;
  if (print_decays) {
    fprintf(stdout, "Calculated S1 decay rates (ps^{-1}):");
    for (unsigned i = 0; i < n_car; i++) {
      for (unsigned j = 0; j < n_s_cars[i]; j++) {
        std::vector<size_t> subs = ind2sub(j, cars[i].get_pop_extents());
        fprintf(stdout, "%2u (%1lu, %1lu, %1lu) = %10.6e\n",
                i, subs[0] + 1, subs[1], subs[2], car_decays[i][j]);
      }
    }
  }

  bool print_i_xa = true;
  if (print_i_xa) {
    for (unsigned k = 0; k < n_car; k++) {
      for (unsigned j = 0; j < k_chl_car[k].size(); j = j + 2) {
        std::vector<size_t> subs = ind2sub(j, {n_chl, n_s_cars[k], 2});
        std::vector<size_t> xa = ind2sub(subs[2],
            cars[k].get_pop_extents());
        fprintf(stdout, "%4u (%1lu)<->(%1lu)(%1lu %1lu %1lu)"
            " %10.6e %10.6e\n", 
            j, subs[0], subs[1], xa[0], xa[1], xa[2],
            k_chl_car[k][j], k_chl_car[k][j + 1]);
        k_sum += k_chl_car[k][j];
      }
    }
    fprintf(stdout, "sum of rates = %12.8e\n", k_sum);
  }

  double **k_tot = (double **)calloc(n_total, sizeof(double *));
  for (unsigned i = 0; i < n_total; i++) {
    k_tot[i] = (double *)calloc(n_total, sizeof(double));
  }

  /* not finished - should give option to just read in the total rates */
  bool read_car_rates = false;
  if (read_car_rates) {
    fp = fopen("in/car_rates.dat", "rb");

  }

  /* if (hybrid) { */
    k_tot = hybrid_transfer(n_chl, n_car, cars, gamma, Jij,
        k_chl_car, kij, car_decays);
  /* } else { */
    /* k_tot = total_rates(n_chl, vera.intra_rates(), n_car, n_s_car, gamma, Jij, */
    /*         k_chl_car, odep.Tij); */
  /* } */

  strcpy(fn, p->pop_file);
  status = generate_filename(sizeof(fn), fn,
           "populations", "k_tot");
  if(status != 0) {
  } else {
    fp = fopen(fn, "w");
    print_matrix(fp, NULL, n_total, k_tot);

    cl = fclose(fp);
    if (cl != 0) {
    }
  }

  fprintf(stdout, "\n-------------------------\n"
                  "STEADY STATE FLUORESCENCE\n"
                  "-------------------------\n\n");

  void *params = &odep;
  double *p0 = guess(population_guess, boltz, musq, max, p->n_chl);
  double *p_i = (double *)calloc(p->n_chl, sizeof(double));

  bool steady_state = false;
  if (steady_state) {
    p_i = steady_state_populations(p0, params, p->n_chl);
  } else {
    p_i = hybrid_boltz(p->n_chl, n_car, beta, eigvals, cars); 
  }
  fprintf(stdout, "HYBRID BOLTZ\n");
  double boltz_sum = 0.0;
  for (i = 0; i < p->N; i++) {
    if (i < p->n_chl) {
      fprintf(stdout, "%2u %8.5f\n", i, p_i[i]);
    }
    boltz_sum += boltz[i] * musq[i];
  }

  /* write out the normalised steady-state populations */
  strcpy(fn, p->pop_file);
  status = generate_filename(sizeof(fn), fn,
           "populations", "ss_populations");
  if(status != 0) {
    fprintf(stdout, "Cannot find string 'populations' in "
        "pop_file needed to print out steady-state populations. "
        "Something got renamed accidentally?\n");
    exit(EXIT_FAILURE);
  } else {
    fp = fopen(fn, "w");

    fprintf(stdout, "i\t p_i^eq(norm)\t boltz\t\t boltz*|μ^2|\n");
    fprintf(fp, "# i\t p_i^eq(norm)\t boltz\t\t boltz*|μ^2|\n");
    for (i = 0; i < p->n_chl; i++) {
      fprintf(stdout, "%2d\t%+12.8e\t%+12.8e\t%+12.8e\n",
          i, p_i[i], boltz[i],
          (boltz[i] * musq[i]) / boltz_sum);
      fprintf(fp, "%2d\t%+12.8e\t%+12.8e\t%+12.8e\n",
          i, p_i[i], boltz[i],
          (boltz[i] * musq[i]) / boltz_sum);
    }
    cl = fclose(fp);
    if (cl != 0) {
        fprintf(stdout, "Failed to close steady state"
            " population file, error no. %d.\n", cl);
        exit(EXIT_FAILURE);
    }
  }
  free(boltz); 

  /* fluorescence spectrum */
  /* need to zero the running integral before FFT! */
  free(integral);
  integral = (double *)calloc(p->tau, sizeof(double));
  for (i = 0; i < p->n_chl; i++) {
    for (unsigned int j = 0; j < p->tau; j++) {
      Ftv[j] = p_i[i] * Ft(eigvals[i],
              gi_array[i][j][0], gi_array[i][j][1],
              lambda[i], (double)j * TOFS,
              1. / (1000000. * gamma[i]));
      in[j][0] = REAL(Ftv[j]);
      in[j][1] = IMAG(Ftv[j]);
    }

    fftw_execute(plan); 
    for (unsigned int j = 0; j < p->tau; j++) {
      chiw[i][j]   = out[j][0] * musq[i];
      integral[j] += out[j][0] * musq[i];
    }

  }
  free(Ftv);

  fftw_free(out); fftw_free(in);
  fftw_destroy_plan(plan);
  fftw_free(ft_out); fftw_free(ft_in);
  fftw_destroy_plan(ft_plan);

  fprintf(stdout, "\nWriting F(w) file\n");
  fp = fopen(p->fw_file, "w");
  for (i = 0; i < p->tau; i++) {
    /* unpack the ordering used by FFTW */
    kd = i * 2. * M_PI / (p->tau);
    fprintf(fp, "%18.10f %18.10f\n", kd / TOFS, 
        (pow(kd, 3.) * integral[i] * (1./ sqrt(p->tau))));
  }
  cl = fclose(fp);
  if (cl != 0) {
      fprintf(stdout, "Failed to close F(w) output file %d.\n", cl);
      exit(EXIT_FAILURE);
  }

  fprintf(stdout, "\nWriting χ^{bar}_i(w) files\n");
  /* replace with chi_bar_ this time */
  for (i = 0; i < p->n_chl; i++) {
    strcpy(fn, p->gi_files[i]);
    status = generate_filename(sizeof(fn), fn, "g_", "chi_bar_");
    if(status != 0) {
      fprintf(stdout, "Cannot find string 'g_' in "
          "g_i filename no. %d, needed to print out chi_i. "
          "Something got renamed accidentally?\n", i);
      break;
    }
    fp = fopen(fn, "w");
    for (j = 0; j < p->tau; j++) {
      /* unpack the ordering used by FFTW */
      kd = j * 2. * M_PI / (p->tau);
      fprintf(fp, "%18.10f %18.10f\n", kd / TOFS, 
          kd * chiw[i][j] * (1./ sqrt(p->tau)));
    }
    cl = fclose(fp);
    if (cl != 0) {
        fprintf(stdout, "Failed to close χ_i(w) "
            "output file no. %d, error no. %d.\n", i, cl);
        exit(EXIT_FAILURE);
    }
  }

  free(p_i);

  fprintf(stdout, "\n---------------\n"
                    "TRANSFER MATRIX\n"
                    "---------------\n\n");
  
  double **Tij_vr, **Tij_vr_inv, **Tij_wr;
  Tij_vr = (double **)calloc(n_total, sizeof(double*));
  Tij_vr_inv = (double **)calloc(n_total, sizeof(double*));
  Tij_wr = (double **)calloc(n_total, sizeof(double*));
  for (i = 0; i < n_total; i++) {
    Tij_vr[i]     = (double *)calloc(n_total, sizeof(double));
    Tij_vr_inv[i] = (double *)calloc(n_total, sizeof(double));
    Tij_wr[i]     = (double *)calloc(n_total, sizeof(double));
  }
  decompose_transfer_matrix(n_total, k_tot, Tij_vr, Tij_vr_inv, Tij_wr);

  strcpy(fn, p->fw_file);
  status = generate_filename(sizeof(fn), fn, "fw", "k_tot");
  if (status == 0) {
    fp = fopen(fn, "w");
    print_matrix(fp, NULL, n_total, k_tot);
    cl = fclose(fp);
    if (cl != 0) {
        fprintf(stdout, "Failed to close k_tot "
            "output file %s, error no. %d.\n", fn, cl);
        exit(EXIT_FAILURE);
    }
  } else {
    fprintf(stdout, "Filename could not be generated - not writing\n");
  }
  status = generate_filename(sizeof(fn), fn, "k_tot", "tij_vr");
  if (status == 0) {
    fp = fopen(fn, "w");
    print_matrix(fp, NULL, n_total, Tij_vr);
    cl = fclose(fp);
    if (cl != 0) {
        fprintf(stdout, "Failed to close tij_vr "
            "output file %s, error no. %d.\n", fn, cl);
        exit(EXIT_FAILURE);
    }
  } else {
    fprintf(stdout, "Filename could not be generated - not writing\n");
  }
  status = generate_filename(sizeof(fn), fn, "tij_vr", "tij_wr");
  if (status == 0) {
    fp = fopen(fn, "w");
    print_matrix(fp, NULL, n_total, Tij_wr);
    cl = fclose(fp);
    if (cl != 0) {
        fprintf(stdout, "Failed to close tij_wr "
            "output file %s, error no. %d.\n", fn, cl);
        exit(EXIT_FAILURE);
    }
  } else {
    fprintf(stdout, "Filename could not be generated - not writing\n");
  }
  status = generate_filename(sizeof(fn), fn, "tij_wr", "tij_vr_inv");
  if (status == 0) {
    fp = fopen(fn, "w");
    print_matrix(fp, NULL, n_total, Tij_vr_inv);
    cl = fclose(fp);
    if (cl != 0) {
        fprintf(stdout, "Failed to close tij_vr_inv "
            "output file %s, error no. %d.\n", fn, cl);
        exit(EXIT_FAILURE);
    }
  } else {
    fprintf(stdout, "Filename could not be generated - not writing\n");
  }

  /* change P(0) to just be on the chls */
  free(p0);

  p0 = (double *)calloc(n_total, sizeof(double));
  double *p0_chl = (double *)calloc(p->n_chl, sizeof(double));
  p0_chl = guess(population_guess,
      bcs(p->n_chl, eigvals, protocol.T), 
      musq, max, p->n_chl);
  /* p0_chl = ; */
  double p0_sum = 0.;
  for (unsigned i = 0; i < p->n_chl; i++) {
    /* p0[i] = 1. / p->n_chl; */
    p0[i + 1] = p0_chl[i];
    p0_sum += p0_chl[i];
    fprintf(stdout, "%2d %10.6f\n", i, p0[i + 1]);
  }
  fprintf(stdout, "sum = %10.6f\n", p0_sum);

  free(p0_chl);

  if (calculate_CD) {
    double **com = (double **)calloc(p->n_chl, sizeof(double *));
    for (unsigned i = 0; i < p->n_chl; i++) {
      com[i] = (double *)calloc(3, sizeof(double));
    }
    com = read_mu(p->com_file, p->n_chl);

    double *cdw = (double *)calloc(p->tau, sizeof(double));
    cd_calc(n_chl, p->tau, chiw, mu, eig, com, eigvals, cdw);

    for (unsigned i = 0; i < p->n_chl; i++) {
      free(com[i]);
    }
    free(com);
    strcpy(fn, p->fw_file);
    status = generate_filename(sizeof(fn), fn, "fw", "cd");
    if(status != 0) {
      fprintf(stdout, "CD filename not constructed, status %d\n", i);
    } else {
      fp = fopen(fn, "w");
      fprintf(stdout, "\nWriting CD file\n\n");
      print_vector(fp, NULL, p->tau, cdw);
      cl = fclose(fp);
      if (cl != 0) {
          fprintf(stdout, "Failed to close CD "
              "output file, error no. %d.\n", cl);
          exit(EXIT_FAILURE);
      }
    }
    free(cdw);
  }

  fprintf(stdout, "\n-------------------\n"
                    "POPULATION DYNAMICS\n"
                    "-------------------\n\n");

  /* note - this could be sped up by making life_yet the condition
   * of a while loop, that way we break out once we reach 1/e on the
   * chlorophylls. also write_pop_to_file currently does nothing - I
   * changed the averaging script so I can delete unwanted files, so
   * it doesn't really matter whether we save the populations here */

  double ti = 0.0; double dt = 1.0; double thresh = 1e-10;
  double *pt = (double *)calloc(n_total, sizeof(double));
  /* use previous step to check convergence */
  double *pt_prev = (double *)calloc(n_total, sizeof(double));
  unsigned int MAX_ITER = 5000; /* 4 ns */
  unsigned int print_pop = 0;
  status = 0;
  double sum = 0.;
  bool life_yet = false;
  /* bool write_pop_to_file = false; */

  FILE *gp;
  strcpy(fn, p->pop_file);
  status = generate_filename(sizeof(fn), fn, "pop", "gs_pop");
  if (status == 0) {
    gp = fopen(fn, "w");
  } else {
    fprintf(stdout, "gs_population file could not be created");
  }

  fprintf(stdout, "\nWriting populations\n");
  fp = fopen(p->pop_file, "w");

  double total_sum = 0.;
  for (i = 0; i < MAX_ITER; i++) {
    total_sum = 0.;
    sum = 0.;
    pt_prev[j] = pt[j];

    ti = (i * dt);
    population(n_total, ti, pt, Tij_vr, Tij_vr_inv, Tij_wr, p0);

    fprintf(fp, "%6.3f ", ti);
    fprintf(gp, "%6.3f ", ti);

    if (print_pop) {
      fprintf(stdout, "ti = %6.3f ", ti);
    }

    for (j = 0; j < n_total; j++) {
      fprintf(fp, "%+12.8e ", pt[j]);

      total_sum += pt[j];
      if (j > 0 && j <= n_chl) {
        sum += pt[j];
      }

      /* if (j == 0) { // redfield gs */
      /*   sum += pt[j]; */
      /* } */
      /* if (hybrid) { */
      /*   if (j == n_chl + 1) { */
      /*     sum += pt[j]; */
      /*   } */
      /*   if (j == n_s_car + n_chl + 2) { */
      /*     sum += pt[j]; */
      /*   } */
      /* } else { */
      /*   if (j > n_chl && j < (n_s_car + n_chl + 1)) { // 620 */
      /*     std::vector<size_t> subs = ind2sub(j - (n_chl + 1), */
      /*                                vera.get_pop_extents()); */
      /*     if (subs[0] == 0) { */
      /*       sum += pt[j]; */
      /*     } */
      /*   } */
      /*   if (j >= (n_s_car + n_chl + 1)) { // 621 */
      /*     std::vector<size_t> subs = ind2sub(j - (n_s_car + n_chl + 1), */
      /*                                vera.get_pop_extents()); */
      /*     if (subs[0] == 0) { */
      /*       sum += pt[j]; */
      /*     } */
      /*   } */
      /* } */

      if (print_pop) {
        fprintf(stdout, "%+12.8e ", pt[j]);
      }

    }

    if (abs(total_sum - 1.) > 1E-8) {
      fprintf(stdout, "step %4u, total sum = %10.6e, loss = %10.6e\n",
          i, total_sum, 1. - total_sum);
    }

    fprintf(gp, "%+12.8e\n", sum);

    if (!life_yet) {
      if (sum < (exp(-1.))) {
        status = generate_filename(sizeof(fn), fn, "gs_populations", "tau");
        if (status == 0) {
          FILE *hp = fopen(fn, "w");
          fprintf(hp, "%+12.8e\n", float(i));
          fprintf(stdout, "<tau> = %+12.8e\n", float(i));
          cl = fclose(hp);
          if (cl != 0) {
              fprintf(stdout, "Failed to close tau "
                  "output file %s, error no. %d.\n", fn, cl);
              exit(EXIT_FAILURE);
          }
        } else {
          fprintf(stdout, "tau file could not be created");
        }
        status = generate_filename(sizeof(fn), fn, "tau", "pop_at_tau");
        if (status == 0) {
          FILE *hp = fopen(fn, "w");
          for (unsigned j = 0; j < n_total; j++) {
            fprintf(hp, "%+12.8e ", pt[j]);
          }
          cl = fclose(hp);
          if (cl != 0) {
              fprintf(stdout, "Failed to close pop_at_tau "
                  "output file %s, error no. %d.\n", fn, cl);
              exit(EXIT_FAILURE);
          }
        } else {
          fprintf(stdout, "pop_at_tau file could not be created");
        }
        life_yet = true;
      }
    }

    fprintf(fp, "\n");
    if(print_pop) {
      fprintf(stdout, ". chl sum = %12.8e\n", sum);
    }

    if (i > 0) {
      unsigned short int converged = pop_converge(pt, pt_prev, n_total, thresh);
      if (converged) {
        fprintf(stdout, "all populations converged to within %f; "
            "ending integration at interation %d\n", thresh, i);
        break;
      }
    }

  }

  status = generate_filename(sizeof(fn), fn, "pop_at_tau", "final_pop");
  if (status == 0) {
    FILE *hp = fopen(fn, "w");
    fprintf(hp, "%+12.8e ", total_sum);
    cl = fclose(hp);
    if (cl != 0) {
        fprintf(stdout, "Failed to close final_pop "
            "output file %s, error no. %d.\n", fn, cl);
        exit(EXIT_FAILURE);
    }
  } else {
    fprintf(stdout, "pop_at_tau file could not be created");
  }

  fclose(fp);
  fclose(gp);

  /* fprintf(stdout, "\n------------\n" */
  /*                   "FORSTER TEST\n" */
  /*                   "------------\n\n"); */

  /* prototype for an array of them */
  /* chromophore *cs[p->N]; */
  /* for (i = 0; i < p->N; i++) { */
  /*   cs[i] = create_chromophore(p->tau); */
  /* } */

  /* double k_610_620 = forster_test(); */
  /* fprintf(stdout, "610-620 Forster rate = %8.6e => %12.8f ps\n", */
  /*         k_610_620, 1. / k_610_620); */


  /* deallocations of 2d stuff */
  for (i = 0; i < n_total; i++) {
    /* free_chromophore(cs[i]); */
    free(Tij_vr[i]);
    free(Tij_vr_inv[i]);
    free(Tij_wr[i]);
    free(k_tot[i]);
  }
  for (i = 0; i < p->n_chl; i++) {
    free(lineshape_files[i]);
    free(eig[i]);
    free(mu[i]);
    free(wij[i]);
    free(kij[i]);
    free(chiw[i]);
    free(pump[i]);
    free(odep.Tij[i]);
    free(normed_ai[i]);
    free(normed_fi[i]);
    fftw_free(gi_array[i]);
  }
  free(mu);
  free(eig);
  free(wij);
  free(kij); 
  free(k_tot); 
  free(odep.Tij);
  free(pt);
  free(pt_prev);
  free(normed_ai);
  free(normed_fi);

  for (i = 0; i < p->N; i++) {
    free(Jij[i]);
  }
  free(Jij);

  fftw_free(gi_array);
  /* fftw_cleanup(); */

  free(p0); free(Tij_vr); free(Tij_vr_inv); free(Tij_wr);
  free(line); free(lineshape_files);
  free(integral);
  free(eigvals); free(gamma); free(lambda);
  free(line_params);
  free(car_decays);
  free(f);
  free(rates);
  free(chiw); free(pump); free(chiw_ints);
  free(p);

  exit(EXIT_SUCCESS);
}
