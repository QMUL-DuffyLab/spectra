#include "hybrid_vera.h"

std::vector<std::vector<double>>
k_i_xa_hybrid(std::vector<VERA> x, unsigned n_chl, unsigned n_car, unsigned tau,
               double **eig, double *eigvals,
               double **Jij, double **normed_ai, double **normed_fi,
               pulse v_abs, double beta, double **car_decays)
{
  std::vector<size_t> vib_totals(n_car, 0.);
  for (unsigned i = 0; i < n_car; i++) {
    vib_totals[i] = pow(x[i].n_vib + 1, x[i].n_normal);
  }
  double ji = 0., ji_work = 0.;
  /* double *abs = (double *)calloc(tau, sizeof(double)); */
  double *fi_ad  = (double *)calloc(tau, sizeof(double));
  double *ai_fd  = (double *)calloc(tau, sizeof(double));
  std::vector<std::vector<double>> k_i_xa(n_car);
  unsigned short print_ji = 0;
  unsigned short print_delta_fc = 0;
  std::vector<std::vector<size_t>> pop_extents;
  for (unsigned i = 0; i < n_car; i++) {
    pop_extents.push_back(x[i].get_pop_extents()); 
  }
  std::vector<std::vector<double>> car_rates;
  for (unsigned i = 0; i < n_car; i++) {
    car_rates.push_back(x[i].intra_rates()); 
  }
  std::vector<size_t> n_s_totals(n_car, 0);
  for (unsigned i = 0; i < n_car; i++) {
    n_s_totals[i] = vib_totals[i] * x[i].n_elec;
  }
  bool print_decay_details = false;
  bool print_details = false;
  bool output_lineshapes = false;

  for (unsigned carotenoid = n_chl;
      carotenoid < n_chl + n_car; carotenoid++) {
    for (unsigned chl_index = 0; chl_index < n_chl; chl_index++) {
      /* assumes the carotenoids are at the end of the arrays;
       * fortran code puts them there so shouldn't be an issue */

      ji = 0.;
      for (unsigned m = 0; m < n_chl; m++) {
        ji += eig[m][chl_index] * Jij[m][carotenoid];
      }

      if (print_ji) {
        fprintf(stdout, "%5u %12.8e\n", chl_index, ji);
      }

      /* outer loop runs from vib_total to 2 * vib_total because
       * we're looping over the first excited electronic state (S1).
       * the upward rate is just the rate from the true g/s to this
       * state; then the inner loop is over the whole set of ground
       * states, to get the effective decay rate. */
      size_t car_index = carotenoid - n_chl;
      for (unsigned ii = vib_totals[car_index];
          ii < 2 * vib_totals[car_index]; ii++) {
        std::vector<size_t> a = ind2sub(ii, pop_extents[car_index]);
        /* we're only considering the transition from the true g/s
         * up to (1, {a}) */
        double e_xa = x[car_index].get_w_elec(a[0]);
        double chl_car = 0.;
        double car_chl = 0.;
        double fc = 1.;

        for (unsigned alpha = 0; alpha < x[car_index].n_normal; alpha++) {
          e_xa += x[car_index].get_w_normal(alpha) * a[alpha + 1];
          fc *= x[car_index].get_fc({0, 1, alpha, a[alpha + 1], 0});
        }
        ji_work = ji * fc;

        /* e_xa is also delta_xy_ba in this case since y_b = {0} */
        /* the S0-S1 width is the first element which is all we need */
        v_abs.width = x[car_index].get_widths(0);
        /* now we have to calculate the rate between every delta_xy_ba */
        v_abs.centre = e_xa;
        double *abs = incident(v_abs, tau);

        for (unsigned step = 0; step < tau; step++) {
          fi_ad[step] = normed_fi[chl_index][step] * abs[step];
          ai_fd[step] = normed_ai[chl_index][step] * abs[step];
        }
        
        chl_car = pow(ji_work, 2.)
                * trapezoid(fi_ad, 2. * M_PI / (TOFS * tau), tau);
        car_chl = pow(ji_work, 2.)
                * trapezoid(ai_fd, 2. * M_PI / (TOFS * tau), tau);

        if (print_delta_fc) {
          fprintf(stdout, "%2u %1u %2u %10.6e %10.6e %10.6e "
              "%10.6e %10.6e %10.6e %10.6e %10.6e %10.6e\n",
              chl_index, carotenoid, ii, e_xa, fc, pow(ji_work, 2.), 
              trapezoid(normed_fi[chl_index], 2. * M_PI / (TOFS * tau), tau),
              trapezoid(normed_ai[chl_index], 2. * M_PI / (TOFS * tau), tau),
              trapezoid(fi_ad, 2. * M_PI / (TOFS * tau), tau),
              trapezoid(ai_fd, 2. * M_PI / (TOFS * tau), tau),
              chl_car,
              car_chl
              );
        }


        /* only calculate this once for each carotenoid lol */
        if (chl_index == 0) {
          for (unsigned kk = 0; kk < vib_totals[car_index]; kk++) {

            /* NB: indexing!!! */
            car_decays[car_index][ii - vib_totals[car_index]] += 
              car_rates[car_index][sub2ind({kk, ii},
                  {n_s_totals[car_index], n_s_totals[car_index]})];
            if (print_decay_details) {
              fprintf(stdout, "%2u %2u %10.6e\n",
                  ii, kk, car_rates[car_index][sub2ind({kk, ii},
                  {n_s_totals[car_index], n_s_totals[car_index]})]);
            }

          } // kk

        }

        if (output_lineshapes) {
          /* both carotenoids have same lineshapes */
          /* this is a nightmare lol strings in C!!!
           * snprintf will add a null char which we don't want,
           * so define another char[2] and memcpy the formatted 
           * number to that without the null char. then gen filename.
           * set the number back to 01 each time by redefining fn,
           * otherwise we'd have to keep track of the previous number
           * and search for that. i guess alternatively we could
           * hardcode the index of the underscore and overwrite
           * two characters after that. but who cares. it's ugly either way
           *
           * also note that this'll fail if the hardcoded directory doesn't
           * exist. this function doesn't have access to the filenames atm
           * so i did it like this, please do not judge lol
           * */
          char snum[3], car_num[2], num[2];
          if (carotenoid == n_chl) {
            char fn[200] = "out/monomer/1/1/1/car_01.dat\0";
            int status = snprintf(snum, 3, "%02lu", ii - vib_totals[car_index]);
            memcpy(car_num, snum, 2);
            status = generate_filename(sizeof(fn), fn, "01", snum);
            if (status != 0) {
              fprintf(stdout, "filename generation for outputting"
                  "car lineshape %d didn't work\n", ii);
              break;
            } else {
              FILE *fp = fopen(fn, "w");
              fprintf(fp, "# chl_car = %10.6e\n",
                      chl_car * CM_PER_PS * PI * 2);
              fprintf(fp, "# car_chl = %10.6e\n",
                      car_chl * CM_PER_PS * PI * 2);

              for (unsigned step = 0; step < tau; step++) {
                fprintf(fp, "%10.6e\n", abs[step]);
              }
              int cl = fclose(fp);
              if (cl != 0) {
                fprintf(stdout, "Failed to close car lineshape"
                    "output file no. %d, error no. %d.\n", ii, cl);
                exit(EXIT_FAILURE);
              }
            }
          }
          if (ii == vib_totals[car_index]) {
            char fn[200] = "out/monomer/1/1/1/chl_01.dat\0";
            int status = snprintf(snum, 3, "%02u", chl_index);
            memcpy(num, snum, 2);
            /* fprintf(stdout, "num = %s", num); */
            status = generate_filename(sizeof(fn), fn, "01", snum);
            /* fprintf(stdout, "%s", fn); */
            if (status != 0) {
              fprintf(stdout, "filename generation for outputting"
                  "chl lineshape %d didn't work\n", chl_index);
              break;
            } else {
              FILE *fp = fopen(fn, "w");
              for (unsigned step = 0; step < tau; step++) {
                fprintf(fp, "%10.6e %10.6e\n",
                    normed_ai[chl_index][step],
                    normed_fi[chl_index][step]);
              }
              int cl = fclose(fp);
              if (cl != 0) {
                fprintf(stdout, "Failed to close chl lineshape"
                    "output file no. %d, error no. %d.\n", chl_index, cl);
                exit(EXIT_FAILURE);
              }
            }
          }
        } // output_lineshapes

        chl_car *= CM_PER_PS * 2. * M_PI;
        car_chl *= CM_PER_PS * 2. * M_PI;
        if (eigvals[chl_index] > e_xa) {
          car_chl *= exp(-beta * (eigvals[chl_index] - e_xa));
        } else {
          chl_car *= exp(-beta * (e_xa - eigvals[chl_index]));
        }
        k_i_xa[car_index].push_back(chl_car);
        k_i_xa[car_index].push_back(car_chl);

        if (print_details) {
          fprintf(stdout, "chl = %2u, car = %2u,"
          " state = (%2lu, %2lu, %2lu): "
          "e_chl = %10.6e, e_car = %10.6e, "
          "FC = %10.6e, chl->car = %10.6e,"
          " car->chl = %10.6e\n", chl_index, carotenoid,
          a[0], a[1], a[2],
          eigvals[chl_index], e_xa,
          fc, chl_car, car_chl);
        }

        free(abs);

      } // ii
    }
  } // chl_index
  /* free(abs); */
  free(fi_ad);
  free(ai_fd);

  return k_i_xa;
}

double**
hybrid_transfer(unsigned n_chl, unsigned n_car, std::vector<VERA> x,
    double *gamma, double **Jij, std::vector<std::vector<double>> k_i_delta,
    double **redfield_rates, double **car_decays)
{
  /* n_vib_tot is the total number of vibrational levels
   * per electronic state; we only need this here because
   * we're not concerned with hot ground states/S2 */
  std::vector<size_t> n_vib_tot(n_car, 0);
  size_t size = n_chl + 1;
  for (unsigned i = 0; i < n_car; i++) {
    n_vib_tot[i] = pow(x[i].n_vib + 1, x[i].n_normal);
    size += n_vib_tot[i] + 1;
  }

  double **k_tot = (double **)calloc(size, sizeof(double*));
  for (unsigned i = 0; i < size; i++) {
    k_tot[i] = (double *)calloc(size, sizeof(double));
  }

  std::vector<std::vector<size_t>> car_extents;
  for (unsigned i = 0; i < n_car; i++) {
    car_extents.push_back(x[i].get_pop_extents());
  }

  /* NB: this assumes n_vib_tot is the same for each carotenoid!
   * this could be fixed pretty easily - make k_i_xa_hybrid return
   * a vector of vectors, one for each carotenoid, then make the 
   * below a vector of vectors in the same way
   */
  /* std::vector<size_t> chl_car_extents = {n_chl, n_car, n_vib_tot[0], 2}; */
  std::vector<std::vector<size_t>> chl_car_extents;
  for (unsigned i = 0; i < n_car; i++) {
    chl_car_extents.push_back({n_chl, n_vib_tot[i], 2});
  }

  std::vector<std::vector<double>> car_rates;
  for (unsigned i = 0; i < n_car; i++) {
    car_rates.push_back(x[i].intra_rates());
  }

  std::vector<size_t> car_gs_indices(n_car, 0);
  car_gs_indices[0] = n_chl + 1;
  for (unsigned k = 1; k < n_car; k++) {
    car_gs_indices[k] = car_gs_indices[k - 1] + n_vib_tot[k - 1] + 1;
  }

  std::vector<bool> i_in_car(n_car, false);
  std::vector<bool> j_in_car(n_car, false);

  unsigned rgs = 0;
  for (unsigned i = 0; i < size; i++) {
    for (unsigned j = 0; j < size; j++) {
      for (unsigned k = 0; k < n_car; k++) {
        /* need to reset each time */
        i_in_car[k] = false;
        j_in_car[k] = false;
      }

      for (unsigned k = n_car; k--;) {
        /* this looks weird but bear with me:
         * we want it to be true if i / j correspond to an actual
         * vibronic state on the carotenoid. if they correspond to the
         * ad-hoc ground state index, we don't want these to be true;
         * we deal with that case separately */
        if (i == car_gs_indices[k]) break;
        /* also: cast both to int bc otherwise you get wraparound */
        if ((int)i > (int)car_gs_indices[k]) {
          i_in_car[k] = true;
          break;
        }
      }

      for (unsigned k = n_car; k--;) {
        if (j == car_gs_indices[k]) break;
        if ((int)j > (int)car_gs_indices[k]) {
          j_in_car[k] = true;
          break;
        }
      }

      bool i_rgs  = (i == rgs);
      bool j_rgs  = (j == rgs);
      bool i_chls = (i > rgs && i <= n_chl);
      bool j_chls = (j > rgs && j <= n_chl);

      if (i == j) continue; /* deal with this at the end */

      if (i_rgs) {
        if (j_rgs) {
          continue;
        }
        if (j_chls) {
          k_tot[i][j] += (1. / (1000 * gamma[j - 1]));
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
          k_tot[j][i] = redfield_rates[i - 1][j - 1];
        }
      } // i_chls

      for (unsigned k = 0; k < n_car; k++) {
        if (i_chls) {
          if (j_in_car[k]) {
            k_tot[j][i] = k_i_delta[k][sub2ind({i - 1,
                j - (car_gs_indices[k] + 1), 0}, chl_car_extents[k])];

            /* fprintf(stdout, "i_chls, j_in_car: "); */
            /* fprintf(stdout, "%2u %2u %2lu %2u %10.6e %10.6e\n", */
            /*     j, i, j - (car_gs_indices[k] + 1), i - 1, */
            /*     k_i_delta[k][sub2ind({i - 1, */
            /*     j - (car_gs_indices[k] + 1), 0}, chl_car_extents[k])], */
            /*     k_tot[j][i]); */
          }
        }

        if (i == car_gs_indices[k]) {
          if (j_in_car[k]) {
            /* fprintf(stdout, "i==car_gs, j_in_car: "); */
            /* fprintf(stdout, "%2u %2u %2u %2lu %10.6e\n", */
            /*     i, j, i - 1, j - (car_gs_indices[k] + 1), */
            /*     car_decays[k][j - (car_gs_indices[k] + 1)]); */
            k_tot[i][j] += car_decays[k][j - (car_gs_indices[k] + 1)];
          }
        }

        if (i_in_car[k]) {
          if (j_chls) {

            k_tot[j][i] = k_i_delta[k][sub2ind({j - 1,
                i - (car_gs_indices[k] + 1), 1}, chl_car_extents[k])];

            /* fprintf(stdout, "i_in_car, j_chls: "); */
            /* fprintf(stdout, "%2u %2u %2u %2lu %10.6e %10.6e\n", */
            /*     j, i, j - 1, i - (car_gs_indices[k] + 1), */
            /*     k_i_delta[k][sub2ind({j - 1, */
            /*     i - (car_gs_indices[k] + 1), 1}, chl_car_extents[k])], */
            /*     k_tot[j][i]); */
          }
          if (j == car_gs_indices[k]) {
            /* fprintf(stdout, "i_in_car, j==car_gs: "); */
            /* fprintf(stdout, "%2u %2u %2u %2lu %10.6e\n", */
            /*     i, j, i - 1, j - (car_gs_indices[k] + 1), */
            /*     car_decays[k][j - (car_gs_indices[k] + 1)]); */
            k_tot[j][i] = car_decays[k][i - (car_gs_indices[k] + 1)];
          }

          if (j_in_car[k]) {
            std::vector<size_t> i_subs = ind2sub(i - (car_gs_indices[k] + 1),
                car_extents[k]);
            std::vector<size_t> j_subs = ind2sub(j - (car_gs_indices[k] + 1),
                car_extents[k]);
            /* S1! */
            i_subs[0]++; j_subs[0]++;
            size_t i_index = sub2ind(i_subs, car_extents[k]);
            size_t j_index = sub2ind(j_subs, car_extents[k]);
            k_tot[i][j] += car_rates[k][sub2ind({i_index, j_index},
                {x[k].n_total, x[k].n_total})];
            /* fprintf(stdout, "i_in_car, j_in_car: %2u %2u %10.6e %10.6e\n", */
            /*     i, j, car_rates[k][sub2ind({i_index, j_index}, */
            /*     {x[k].n_total, x[k].n_total})], k_tot[j][i]); */
          }
        }
      }

    } // j
  }

  /* now deal with the diagonals */
  for (unsigned i = 0; i < size; i++) {

    bool car_gs = false;
    /* this could easily be achieved by starting the loop
     * at 1 but i've kept it this way for clarity */
    if (i == rgs) continue;

    for (unsigned k = 0; k < n_car; k++) {
      if (i == car_gs_indices[k]) car_gs = true;;
    }

    if (car_gs) continue;

    for (unsigned j = 0; j < size; j++) {
      if (i == j) {
        continue;
      } else {
        k_tot[i][i] -= k_tot[j][i];
      }
    }
  }

  return k_tot;

}

double*
hybrid_boltz(unsigned n_chl, unsigned n_car, double beta,
             double *eigvals, std::vector<VERA> x)
{
  /* get the boltzmann factor for every state and normalise, 
   * then return just the chlorophyll boltzmann factors. we do this
   * because the carotenoid doesn't fluoresce but does have some
   * population on it. */
  std::vector<size_t> n_vib_tot(n_car, 0);
  std::vector<size_t> car_starts(n_car, 0);
  size_t size = n_chl;
  car_starts[0] = n_chl;
  for (unsigned i = 0; i < n_car; i++) {
    n_vib_tot[i] = pow(x[i].n_vib + 1, x[i].n_normal);
    size += n_vib_tot[i];
    if (i > 0) {
      car_starts[i] = car_starts[i - 1] + n_vib_tot[i - 1];
    }
  }
  double *work = (double *)calloc(size, sizeof(double));
  double *boltz = (double *)calloc(n_chl, sizeof(double));
  std::vector<std::vector<size_t>> extents;
  for (unsigned i = 0; i < n_car; i++) {
    extents.push_back(x[i].get_pop_extents());
  }
  double work_sum = 0.;

  std::vector<bool> i_in_car(n_car, false);
  for (unsigned i = 0; i < size; i++) {
    if (i < n_chl) {
      work[i] = exp(-beta * eigvals[i]); // dimensions!
      work_sum += work[i];
    } else {
      for (unsigned k = 0; k < n_car; k++) {
        i_in_car[k] = false;
      }
      for (unsigned k = n_car; k--;) {
        if ((int)i >= (int)car_starts[k]) {
          i_in_car[k] = true;
          break;
        }
      }
      for (unsigned k = 0; k < n_car; k++) {
        if (i_in_car[k]) {
          std::vector<size_t> subs = ind2sub(i - car_starts[k], extents[k]);
          subs[0]++;
          double e_xa = x[k].get_w_elec(subs[0]);
          for (unsigned alpha = 0; alpha < x[k].n_normal; alpha++) {
            e_xa += x[k].get_w_normal(alpha) * subs[alpha + 1];
          }

          work[i] = exp(-beta * e_xa);
          work_sum += work[i];
          fprintf(stdout, "%2u %8.5e %8.5e\n", i, e_xa, work_sum);
        }
      }

    }
  }

  for (unsigned i = 0; i < size; i++) {
    work[i] /= work_sum;
  }
  for (unsigned i = 0; i < n_chl; i++) {
    boltz[i] = work[i];
  }
  free(work);
  return boltz;

}
