#include "input.h"
#include "steady_state.h"
#include "forster.h"
#include <fftw3.h>
#include <stdio.h>
#include <gsl/gsl_eigen.h>

int
main(int argc, char** argv)
{
  FILE *fp;
  unsigned int i, j;
  int status;
  char *line, **lineshape_files;
  double kd;
  double complex **gi_array;
  double *eigvals, *gamma, *rates, *musq, *chi_p,
         *lambda, *integral, *chiw_ints, *ww,
         **wij, **kij, **Jij, **mu, **eig, **chiw, **pump;
  Parameters *line_params;
  fftw_complex *out, *in;
  fftw_plan plan;

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
  pulse pump_properties = { .type=DELTA, .centre=15000., .width=300. };
  /* form of P_i(0) - check steady_state.h for details */
  ss_init population_guess = MUSQ;

  /* malloc 1d stuff, read them in */
  eigvals     = calloc(p->N, sizeof(double));
  gamma       = calloc(p->N, sizeof(double));
  rates       = calloc(p->N, sizeof(double));
  musq        = calloc(p->N, sizeof(double));
  lambda      = calloc(p->N, sizeof(double));
  line_params = malloc(p->N * sizeof(Parameters));
  chiw_ints   = calloc(p->N, sizeof(double));
  ww          = calloc(p->tau, sizeof(double));
  integral    = calloc(p->tau, sizeof(double));
  in          = (fftw_complex*) fftw_malloc(
                 sizeof(fftw_complex) * p->tau);
  out         = (fftw_complex*) fftw_malloc(
                 sizeof(fftw_complex) * p->tau);

  /* malloc 2d stuff */
  lineshape_files = malloc(p->N * sizeof(char*));
  mu              = calloc(p->N, sizeof(double*));
  gi_array        = calloc(p->N, sizeof(double complex*));
  eig             = calloc(p->N, sizeof(double*));
  wij             = calloc(p->N, sizeof(double*));
  kij             = calloc(p->N, sizeof(double*));
  Jij             = calloc(p->N, sizeof(double*));
  chiw            = calloc(p->N, sizeof(double*));
  pump            = calloc(p->N, sizeof(double*));

  for (i = 0; i < p->N; i++) {
    lineshape_files[i] = malloc(200 * sizeof(char));
    gi_array[i]        = calloc(p->tau, sizeof(double complex));
    mu[i]              = calloc(3, sizeof(double));
    eig[i]             = calloc(p->N, sizeof(double));
    wij[i]             = calloc(p->N, sizeof(double));
    kij[i]             = calloc(p->N, sizeof(double));
    Jij[i]             = calloc(p->N, sizeof(double));
    chiw[i]            = calloc(p->tau, sizeof(double));
    pump[i]            = calloc(p->tau, sizeof(double));
  }

  /* read function reads a set of doubles: here we read in
   * the rates (gamma), reorganisation energies (lambda) 
   * and exciton energies (eigvals - they're the eigenvalues
   * of the diagonalised Hamiltonian */
  gamma   = read(p->gamma_file, p->N);
  lambda  = read(p->lambda_file, p->N);
  eigvals = read(p->eigvals_file, p->N);

  /* read 2d stuff */
  gi_array = read_gi(p->gi_files, p->N, p->tau);
  eig      = read_eigvecs(p->eigvecs_file, p->N);
  mu       = read_mu(p->mu_file, p->N);
  ww       = incident(pump_properties, p->tau);

  fp = fopen(argv[3], "r"); /* read in list of lineshape files here */
  line = malloc(200 * sizeof(char));
  for (i = 0; i < p->N; i++) {
    /* now load in the parameters for each ligand */
    fgets(line, 200, fp);
    line[strcspn(line, "\n")] = 0;
    strcpy(lineshape_files[i], line);
    line_params[i] = get_parameters(lineshape_files[i],
                     protocol.chl_ansatz);
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

  for (i = 0; i < p->N; i++) {
    musq[i] = pow(mu[i][0], 2.) + pow(mu[i][1], 2.) + pow(mu[i][2], 2.);

    for (unsigned int j = 0; j < p->tau; j++) {
      in[j] = At(eigvals[i], creal(gi_array[i][j]), cimag(gi_array[i][j]),
                 (double)j * TOFS,
                 1. / (1000000. * gamma[i]));
    }

    fftw_execute(plan); 
    for (unsigned int j = 0; j < p->tau; j++) {
      chiw[i][j] = creal(out[j]) * musq[i];
      pump[i][j] = chiw[i][j] * ww[j];
      integral[j] += creal(out[j]) * musq[i];
    }

    chi_p = pump[i];
    chiw_ints[i] = trapezoid(chi_p, p->tau);

  }

  fprintf(stdout, "\nWriting A(w) file\n");
  fp = fopen(p->aw_file, "w");
  for (i = 0; i < p->tau; i++) {
    /* unpack the ordering used by FFTW */
    kd = i * 2. * M_PI / (p->tau);
    fprintf(fp, "%18.10f %18.10f\n", kd / TOFS, 
        kd * creal(integral[i]) * (1./ sqrt(p->tau)) * 6.4);
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
          kd * creal(chiw[i][j]) * (1./ sqrt(p->tau)) * 6.4);
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
  double *f = calloc(p->N, sizeof(double));
  double *y = calloc(p->N, sizeof(double));
  /* check convergence */
  double *yprev = calloc(p->N, sizeof(double));
  double *boltz = calloc(p->N, sizeof(double));

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

  /* fluorescence spectrum */
  /* need to zero the running integral before FFT! */
  free(integral);
  integral = calloc(p->tau, sizeof(double));
  for (i = 0; i < p->N; i++) {
    for (unsigned int j = 0; j < p->tau; j++) {
      in[j] = p_i_equib[i] * Ft(eigvals[i],
              creal(gi_array[i][j]), cimag(gi_array[i][j]),
              lambda[i], (double)j * TOFS,
              1. / (1000000. * gamma[i]));
    }

    fftw_execute(plan); 
    for (unsigned int j = 0; j < p->tau; j++) {
      chiw[i][j]   = creal(out[j]) * musq[i];
      integral[j] += creal(out[j]) * musq[i];
    }

  }

  fprintf(stdout, "\nWriting F(w) file\n");
  fp = fopen(p->fw_file, "w");
  for (i = 0; i < p->tau; i++) {
    /* unpack the ordering used by FFTW */
    kd = i * 2. * M_PI / (p->tau);
    fprintf(fp, "%18.10f %18.10f\n", kd / TOFS, 
        (pow(kd, 3.) * creal(integral[i]) * (1./ sqrt(p->tau))) * 6.4);
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
          kd * creal(chiw[i][j]) * (1./ sqrt(p->tau)) * 6.4);
    }
    cl = fclose(fp);
    if (cl != 0) {
        fprintf(stdout, "Failed to close χ_i(w) "
            "output file no. %d, error no. %d.\n", i, cl);
        exit(EXIT_FAILURE);
    }
  }

  free(p_i_equib);

  fprintf(stdout, "\n-------------------\n"
                    "EXCITATION LIFETIME\n"
                    "-------------------\n\n");

  double **Tij_vr, **Tij_vr_inv, **Tij_wr;
  Tij_vr = calloc(p->N, sizeof(double*));
  Tij_vr_inv = calloc(p->N, sizeof(double*));
  Tij_wr = calloc(p->N, sizeof(double*));
  for (i = 0; i < p->N; i++) {
    Tij_vr[i] = calloc(p->N, sizeof(double));
    Tij_vr_inv[i] = calloc(p->N, sizeof(double));
    Tij_wr[i] = calloc(p->N, sizeof(double));
  }
  decompose_transfer_matrix(p->N, odep.Tij, Tij_vr, Tij_vr_inv, Tij_wr);

  /* p0 is still our guess from earlier */
  print_vector(stdout, "P(0)", p->N, p0);
  double excite = mean_excitation_lifetime(p->N, Tij_vr,
                                           Tij_vr_inv,
                                           Tij_wr, p0);
  fprintf(stdout, "\n<τ> (ps) = %12.8e\n", excite);

  strcpy(fn, p->fw_file);
  status = generate_filename(sizeof(fn), fn, "fw", "tau");
  if (status == 0) {
    fp = fopen(fn, "w");
    fprintf(fp, "%12.8e", excite);
    cl = fclose(fp);
    if (cl != 0) {
        fprintf(stdout, "Failed to close <τ> "
            "output file %s, error no. %d.\n", fn, cl);
        exit(EXIT_FAILURE);
    }
  } else {
    fprintf(stdout, "Filename could not be generated - not writing\n");
  }

  fprintf(stdout, "\n-------------------\n"
                    "POPULATION DYNAMICS\n"
                    "-------------------\n\n");

  double ti = 0.0; double dt = 1.0; double thresh = 1e-10;
  double *pt = calloc(p->N, sizeof(double));
  /* use previous step to check convergence */
  double *pt_prev = calloc(p->N, sizeof(double));
  unsigned int MAX_ITER = 20000; /* 20 ns */
  unsigned int print_pop = 0;
  status = 0;
  double sum = 0.;

  fprintf(stdout, "\nWriting populations\n");
  fp = fopen(p->pop_file, "w");

  for (i = 0; i < MAX_ITER; i++) {
    sum = 0.;
    pt_prev[j] = pt[j];

    ti = (i * dt);
    population(p->N, ti, pt, Tij_vr, Tij_vr_inv, Tij_wr, p0);

    fprintf(fp, "%6.3f ", ti);
    if (print_pop) {
      fprintf(stdout, "ti = %6.3f ", ti);
    }
    for (j = 0; j < p->N; j++) {
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
      unsigned short int converged = pop_converge(pt, pt_prev, p->N, thresh);
      if (converged) {
        fprintf(stdout, "all populations converged to within %f; "
            "ending integration at interation %d\n", thresh, i);
        break;
      }
    }

  }
  fclose(fp);

  fprintf(stdout, "\n------------\n"
                    "FORSTER TEST\n"
                    "------------\n\n");

  fprintf(stdout, "610-620 Forster rate = %8.6e\n", forster_test());

  free(p0); free(pt); free(Tij_vr); free(Tij_vr_inv); free(Tij_wr);
  free(line); free(lineshape_files); free(integral);
  free(gi_array); free(eigvals); free(gamma); free(lambda); free(mu);
  free(eig); free(wij); free(kij); free(p); free(line_params);
  free(in); free(out); free(y); free(f); free(boltz); free(yprev);
  free(ww);
  exit(EXIT_SUCCESS);
}
