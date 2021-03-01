#include <fftw3.h>
#include <stdio.h>
#include <cmath>
#include <complex>
#include <gsl/gsl_eigen.h>
#include "input.h"
#include "steady_state.h"
#include "forster.h"
#include "vera.h"

int
main(int argc, char** argv)
{
  FILE *fp;
  unsigned int i, j;
  int status;
  char *line, **lineshape_files;
  double kd;
  fftw_complex **gi_array;
  double *eigvals, *ei, *gamma, *rates, *musq, *chi_p,
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

  /* specifies form of incident light for source term in P_i eqns */
  pulse pump_properties = {
    .type=DELTA,
    .target_state = 2,
    .amplitude = 100.,
    .centre=15000.,
    .width=300.,
    .t_peak = 0.1,
    .duration = 70.0E-3
  };

  pulse VERA_pump = {
    .type = GAUSSIAN,
    .target_state = 2,
    .amplitude = 100.,
    .centre = 20000.,
    .width = 10.,
    .t_peak = 0.1,
    .duration = 70.0E-3,
  };

  pulse VERA_absorption = {
    .type = LORENTZIAN,
    .target_state = 0,
    .amplitude = 0.,
    .centre = 0.,
    .width = 0.,
    .t_peak = 0.1,
    .duration = 70.0E-3,
  };

  /* form of P_i(0) - check steady_state.h for details */
  ss_init population_guess = MUSQ;

  /* malloc 1d stuff, read them in */
  eigvals     = (double *)calloc(p->N, sizeof(double));
  gamma       = (double *)calloc(p->N, sizeof(double));
  rates       = (double *)calloc(p->N, sizeof(double));
  musq        = (double *)calloc(p->N, sizeof(double));
  lambda      = (double *)calloc(p->N, sizeof(double));
  line_params = (Parameters *)malloc(p->N * sizeof(Parameters));
  chiw_ints   = (double *)calloc(p->N, sizeof(double));
  ww          = (double *)calloc(p->tau, sizeof(double));
  integral    = (double *)calloc(p->tau, sizeof(double));
  in          = (fftw_complex *)fftw_malloc(p->tau * sizeof(fftw_complex));
  out         = (fftw_complex *)fftw_malloc(p->tau * sizeof(fftw_complex));
  ft_in          = (fftw_complex *)fftw_malloc(p->tau * sizeof(fftw_complex));
  ft_out         = (fftw_complex *)fftw_malloc(p->tau * sizeof(fftw_complex));

  /* malloc 2d stuff */
  lineshape_files = (char **)malloc(p->N * sizeof(char*));
  mu              = (double **)calloc(p->N, sizeof(double*));
  gi_array        = (fftw_complex**)
                    fftw_malloc(p->N * sizeof(fftw_complex*));
  eig             = (double **)calloc(p->N, sizeof(double*));
  wij             = (double **)calloc(p->N, sizeof(double*));
  kij             = (double **)calloc(p->N, sizeof(double*));
  Jij             = (double **)calloc(p->N, sizeof(double*));
  chiw            = (double **)calloc(p->N, sizeof(double*));
  pump            = (double **)calloc(p->N, sizeof(double*));
                      
  for (i = 0; i < p-> N; i++) {
    lineshape_files[i] = (char *)malloc(200 * sizeof(char));
    gi_array[i]        = (fftw_complex*)
                          fftw_malloc(p->N * sizeof(fftw_complex));
    mu[i]              = (double *)calloc(3, sizeof(double));
    eig[i]             = (double *)calloc(p->N, sizeof(double));
    wij[i]             = (double *)calloc(p->N, sizeof(double));
    kij[i]             = (double *)calloc(p->N, sizeof(double));
    Jij[i]             = (double *)calloc(p->N, sizeof(double));
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
  /* this'll need adding to the input struct */
  char jij_file[200] = "out/7_PROD_1/1000/J_ij.out";

  /* read 2d stuff */
  gi_array = read_gi(p->gi_files, p->N, p->tau);
  eig      = read_eigvecs(p->eigvecs_file, p->N);
  Jij      = read_eigvecs(jij_file, p->N);
  mu       = read_mu(p->mu_file, p->N);
  ww       = incident(pump_properties, p->tau);

  fp = fopen(argv[3], "r"); /* read in list of lineshape files here */
  line = (char *)malloc(200 * sizeof(char));
  for (i = 0; i < p->N; i++) {
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

    for (j = 0; j < p->N; j++) {
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

  kij = rate_calc(p->N, eig, wij, line_params);
  check_detailed_balance(p->N, protocol.T, 1e-5, kij, wij);
  rates = relaxation_rates(p->N, gamma);
  char fn[200];
  strcpy(fn, p->fw_file);
  status = generate_filename(sizeof(fn), fn, "fw", "rates");
  fp = fopen(fn, "w");
  print_matrix(fp, NULL, p->N, kij);
  cl = fclose(fp);
  if (cl != 0) {
    fprintf(stdout, "Failed to close list of lineshape files %d.\n", cl);
    exit(EXIT_FAILURE);
  }


  plan = fftw_plan_dft_1d(p->tau, 
  	 in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  ft_plan = fftw_plan_dft_1d(p->tau, 
  	 ft_in, ft_out, FFTW_BACKWARD, FFTW_ESTIMATE);

  double **normed_ai = (double **)calloc(p->N, sizeof(double*));
  double **normed_fi = (double **)calloc(p->N, sizeof(double*));
  double ai_sum = 0., fi_sum = 0.;
  for (unsigned i = 0; i < p->N; i++) {
    normed_ai[i] = (double *)calloc(p->tau, sizeof(double));
    normed_fi[i] = (double *)calloc(p->tau, sizeof(double));
  }

  COMPLEX *Atv = (COMPLEX *)calloc(p->tau, sizeof(COMPLEX));
  COMPLEX *Ftv = (COMPLEX *)calloc(p->tau, sizeof(COMPLEX));
  for (i = 0; i < p->N; i++) {
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
      ai_sum += out[j][0];
      fi_sum += ft_out[j][0];
    }
    for (unsigned int j = 0; j < p->tau; j++) {
      normed_ai[i][j] = out[j][0] / ai_sum;
      normed_fi[i][j] = ft_out[j][0] / fi_sum;
    }


    chi_p = pump[i];
    chiw_ints[i] = trapezoid(chi_p, p->tau);

  }
  free(Atv);
  free(ww);

  fprintf(stdout, "\nWriting A(w) file\n");
  fp = fopen(p->aw_file, "w");
  for (i = 0; i < p->tau; i++) {
    /* unpack the ordering used by FFTW */
    kd = i * 2. * M_PI / (p->tau);
    fprintf(fp, "%18.10f %18.10f\n", kd / TOFS, 
        kd * integral[i] * (1./ sqrt(p->tau)) * 6.4);
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
  for (i = 0; i < p->N; i++) {
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
          kd * chiw[i][j] * (1./ sqrt(p->tau)) * 6.4);
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
  for (i = 0; i < p->N; i++) {
    fprintf(stdout, "%7d %10.6e %10.6e\n", i + 601, musq[i], chiw_ints[i]);
    if (musq[i] > musq_max) {
      max = i;
      musq_max = musq[i];
      musq_sum += musq[i];
      chiw_sum += chiw_ints[i];
    }
  }

  ode_params odep;
  odep.N = p->N;
  odep.kij = kij;
  odep.rates = rates;
  odep.chiw = chiw_ints;
  odep.Tij = transfer_matrix(p->N, rates, kij);
  double *f = (double *)calloc(p->N, sizeof(double));
  double *y = (double *)calloc(p->N, sizeof(double));
  /* check convergence */
  double *yprev = (double *)calloc(p->N, sizeof(double));
  double *boltz = (double *)calloc(p->N, sizeof(double));

  boltz = bcs(p->N, eigvals, protocol.T);
  fprintf(stdout, "\n-----------------\nBOLTZMANN WEIGHTS\n"
                  "-----------------\n\n");
  fprintf(stdout, "Pigment    p_i   |μ^2|*p_i\n");
  for (i = 0; i < p->N; i++) {
    fprintf(stdout, "%7d %8.6f %8.6f\n", i + 601, boltz[i],
        boltz[i] * musq[i]);
    /* possible initial values for transient absorption? */
    if (i == max) {
      y[i] = 1.0;
    } else {
      y[i] = 0.0;
    }
  }
  fprintf(stdout, "\n");

  fprintf(stdout, "\n-------------------------\n"
                  "STEADY STATE FLUORESCENCE\n"
                  "-------------------------\n\n");

  void *params = &odep;
  double *p0 = guess(population_guess, boltz, musq, max, p->N);

  double *p_i_equib = steady_state_populations(p0, params, p->N);
  double boltz_sum = 0.0;
  for (i = 0; i < p->N; i++) {
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
    for (i = 0; i < p->N; i++) {
      fprintf(stdout, "%2d\t%+12.8e\t%+12.8e\t%+12.8e\n",
          i, p_i_equib[i], boltz[i],
          (boltz[i] * musq[i]) / boltz_sum);
      fprintf(fp, "%2d\t%+12.8e\t%+12.8e\t%+12.8e\n",
          i, p_i_equib[i], boltz[i],
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
  for (i = 0; i < p->N; i++) {
    for (unsigned int j = 0; j < p->tau; j++) {
      Ftv[j] = p_i_equib[i] * Ft(eigvals[i],
              gi_array[i][j][0], gi_array[i][j][1],
              lambda[i], (double)j * TOFS,
              1. / (1000000. * gamma[i]));
      in[j][0] = REAL(Ftv[j]);
      in[j][1] = IMAG(Ftv[j]);
      /* in[j][1] = p_i_equib[i] * *Ft(eigvals[i], */
      /*         gi_array[i][j][0], gi_array[i][j][1], */
      /*         lambda[i], (double)j * TOFS, */
      /*         1. / (1000000. * gamma[i]))[1]; */
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
        (pow(kd, 3.) * integral[i] * (1./ sqrt(p->tau))) * 6.4);
  }
  cl = fclose(fp);
  if (cl != 0) {
      fprintf(stdout, "Failed to close F(w) output file %d.\n", cl);
      exit(EXIT_FAILURE);
  }

  fprintf(stdout, "\nWriting χ^{bar}_i(w) files\n");
  /* replace with chi_bar_ this time */
  for (i = 0; i < p->N; i++) {
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
          kd * chiw[i][j] * (1./ sqrt(p->tau)) * 6.4);
    }
    cl = fclose(fp);
    if (cl != 0) {
        fprintf(stdout, "Failed to close χ_i(w) "
            "output file no. %d, error no. %d.\n", i, cl);
        exit(EXIT_FAILURE);
    }
  }

  free(p_i_equib);

  fprintf(stdout, "\n----------\n"
                    "VERA RATES\n"
                    "----------\n\n");

  VERA vera = create_VERA_from_file("in/vera_params.dat");

  size_t n_chl    = 14; /* number of chlorophylls */
  size_t n_car    = 2;
  size_t n_s_car  = 48;
  size_t n_total  = n_chl + 1 + (n_car * n_s_car);
  double beta     = 1. / protocol.T;
  std::vector<double> k_chl_car = k_i_xa(vera, n_chl,
      n_car, p->tau, eig, eigvals, Jij, normed_ai, normed_fi,
      VERA_absorption, beta);
  double k_sum = 0.;

  for (unsigned i = 0; i < k_chl_car.size(); i = i + 2) {
    std::vector<size_t> subs = ind2sub(i, {14, 2, 48, 2});
    std::vector<size_t> xa = ind2sub(subs[2],
        vera.get_pop_extents());
    fprintf(stdout, "%4u (%1u)<->(%1u)(%1u %1u %1u)"
        " %10.6e %10.6e\n", 
        i, subs[0], subs[1], xa[0], xa[1], xa[2],
        k_chl_car[i], k_chl_car[i + 1]);
    k_sum += k_chl_car[i];

  }
  fprintf(stdout, "sum of rates = %12.8e\n", k_sum);

  double **k_tot = (double **)calloc(n_total, sizeof(double *));
  for (unsigned i = 0; i < n_total; i++) {
    double *k_tot = (double *)calloc(n_total, sizeof(double));
  }

  k_tot = total_rates(n_chl, vera, n_car, n_s_car, gamma, Jij,
          k_chl_car, odep.Tij);

  /* for (unsigned i = 0; i < n_total; i++) { */
  /*   for (unsigned j = 0; j < n_total; j++) { */
  /*     fprintf(stdout, "%3u %3u %10.6e\n", i, j, k_tot[i][j]); */
  /*   } */
  /* } */
  
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
  for (unsigned i = 1; i < n_chl + 1; i++) {
    p0[i] = 1./n_chl;
  }

  /* fprintf(stdout, "\n-------------------\n" */
  /*                   "EXCITATION LIFETIME\n" */
  /*                   "-------------------\n\n"); */

  /* okay - for this we need to delete the ground state lines */
  /* n_no_gs = n_total - (n_car + 1); */
  /* double **t_tau = (double **)calloc(n_no_gs, sizeof(double*)); */
  /* double **t_tau_vr = (double **)calloc(n_no_gs, sizeof(double*)); */
  /* double **t_tau_vr_inv = (double **)calloc(n_no_gs, sizeof(double*)); */
  /* double **t_tau_wr = (double **)calloc(n_no_gs, sizeof(double*)); */
  /* double **k_no_gs = (double **)calloc(n_no_gs, sizeof(double*)); */
  /* double *p_no_gs = (double *)calloc(n_no_gs, sizeof(double)); */
  /* for (unsigned i = 0; i < n_no_gs; i++) { */
  /*   t_tau[i] = (double *)calloc(n_no_gs, sizeof(double)); */
  /*   t_tau_vr[i] = (double *)calloc(n_no_gs, sizeof(double)); */
  /*   t_tau_vr_inv[i] = (double *)calloc(n_no_gs, sizeof(double)); */
  /*   t_tau_wr[i] = (double *)calloc(n_no_gs, sizeof(double)); */
  /*   k_no_gs[i] = (double *)calloc(n_no_gs, sizeof(double)); */
  /* } */

  /* for (unsigned i = 0; i < n_total; i++) { */
  /*   if (i != 0 && i != 15 && i != 63) { */
  /*     for (unsigned j = 0; j < n_no_gs; j++) { */
  /*       k_no_gs[i][j] */ 
  /*     } */
  /*   } */
  /* } */
  /* decompose_transfer_matrix(n_total, k_tot, Tij_vr, Tij_vr_inv, Tij_wr); */

  /* double excite = mean_excitation_lifetime(n_total, Tij_vr, */
  /*                                          Tij_vr_inv, */
  /*                                          Tij_wr, p0); */

  /* fprintf(stdout, "\n<τ> (ps) = %12.8e\n", excite); */

  /* status = generate_filename(sizeof(fn), fn, "tij_vr_inv", "tau"); */
  /* status = generate_filename(sizeof(fn), fn, "fw", "tau"); */
  /* if (status == 0) { */
  /*   fp = fopen(fn, "w"); */
  /*   fprintf(fp, "%12.8e", excite); */
  /*   cl = fclose(fp); */
  /*   if (cl != 0) { */
  /*       fprintf(stdout, "Failed to close <τ> " */
  /*           "output file %s, error no. %d.\n", fn, cl); */
  /*       exit(EXIT_FAILURE); */
  /*   } */
  /* } else { */
  /*   fprintf(stdout, "Filename could not be generated - not writing\n"); */
  /* } */

  fprintf(stdout, "\n-------------------\n"
                    "POPULATION DYNAMICS\n"
                    "-------------------\n\n");

  double ti = 0.0; double dt = 1.0; double thresh = 1e-10;
  double *pt = (double *)calloc(n_total, sizeof(double));
  /* use previous step to check convergence */
  double *pt_prev = (double *)calloc(n_total, sizeof(double));
  unsigned int MAX_ITER = 4000; /* 5 ns */
  unsigned int print_pop = 0;
  status = 0;
  double sum = 0.;

  fprintf(stdout, "\nWriting populations\n");
  fp = fopen(p->pop_file, "w");

  for (i = 0; i < MAX_ITER; i++) {
    sum = 0.;
    pt_prev[j] = pt[j];

    ti = (i * dt);
    population(n_total, ti, pt, Tij_vr, Tij_vr_inv, Tij_wr, p0);

    fprintf(fp, "%6.3f ", ti);
    if (print_pop) {
      fprintf(stdout, "ti = %6.3f ", ti);
    }
    for (j = 0; j < n_total; j++) {
      sum += pt[j];
      fprintf(fp, "%+12.8e ", pt[j]);
      if (print_pop) {
        fprintf(stdout, "%+12.8e ", pt[j]);
      }
    }
    fprintf(fp, "\n");
    if(print_pop) {
      fprintf(stdout, ". sum = %12.8e\n", sum);
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
  fclose(fp);

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
  for (i = 0; i < p->N; i++) {
    /* free_chromophore(cs[i]); */
    fftw_free(gi_array[i]);
    free(Tij_vr[i]);
    free(Tij_vr_inv[i]);
    free(Tij_wr[i]);
    free(lineshape_files[i]);
    free(eig[i]);
    free(mu[i]);
    free(wij[i]);
    free(kij[i]);
    free(k_tot[i]);
    free(Jij[i]);
    free(chiw[i]);
    free(pump[i]);
    free(odep.Tij[i]);
    free(normed_ai[i]);
    free(normed_fi[i]);
  }
  /* free(gi_array); */
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
  /* free(ei); */

  fftw_free(gi_array);
  fftw_cleanup();

  free(p0); free(Tij_vr); free(Tij_vr_inv); free(Tij_wr);
  free(line); free(lineshape_files); free(integral);
  free(eigvals); free(gamma); free(lambda);
  free(line_params);
  free(y); free(f); free(yprev); free(rates);
  free(Jij); free(chiw); free(pump); free(chiw_ints);
  free(p);

  exit(EXIT_SUCCESS);
}
