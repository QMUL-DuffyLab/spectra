#include "hybrid_vera.h"

std::vector<double>
k_i_xa_hybrid(VERA x, unsigned n_chl, unsigned n_car, unsigned tau,
               double **eig, double *eigvals,
               double **Jij, double **normed_ai, double **normed_fi,
               pulse v_abs, double beta, double *car_decays)
{
  size_t vib_total = pow(x.n_vib + 1, x.n_normal);
  double ji = 0., ji_work = 0.;
  double *abs = (double *)calloc(tau, sizeof(double));
  double *fi_ad  = (double *)calloc(tau, sizeof(double));
  double *ai_fd  = (double *)calloc(tau, sizeof(double));
  std::vector<double> k_i_xa;
  unsigned short print_ji = 0;
  unsigned short print_delta_fc = 0;
  std::vector<size_t> pop_extents = x.get_pop_extents();
  std::vector<double> car_rates = x.intra_rates();
  size_t n_s_total = vib_total * x.n_elec;
  bool print_decay_details = false;
  bool print_details = false;
  bool output_lineshapes = false;

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

      /* outer loop runs from vib_total to 2 * vib_total because
       * we're looping over the first excited electronic state (S1).
       * the upward rate is just the rate from the true g/s to this
       * state; then the inner loop is over the whole set of ground
       * states, to get the effective decay rate. */
      for (unsigned ii = vib_total; ii < 2 * vib_total; ii++) {
        std::vector<size_t> a = ind2sub(ii, pop_extents);
        /* we're only considering the transition from the true g/s
         * up to (1, {a}) */
        double e_xa = x.get_w_elec(a[0]);
        double chl_car = 0.;
        double car_chl = 0.;
        double fc_sq = 1.;

        for (unsigned alpha = 0; alpha < x.n_normal; alpha++) {
          e_xa += x.get_w_normal(alpha) * a[alpha + 1];
          fc_sq *= pow(x.get_fc({0, 1, alpha,
                a[alpha + 1], 0}), 2.);
        }
        ji_work = ji * fc_sq;

        /* e_xa is also delta_xy_ba in this case since y_b = {0} */
        /* first element is the S0-S1 width which is all we need */
        v_abs.width = x.get_widths(0);
        /* now we have to calculate the rate between every delta_xy_ba */
        v_abs.centre = e_xa;
        abs = incident(v_abs, tau);
        if (output_lineshapes) {
          /* both carotenoids have same lineshapes */
          if (carotenoid == n_chl) {
            /* this is a nightmare lol strings in C!!!
             * snprintf will add a null char which we don't want,
             * so define another char[2] and memcpy the formatted 
             * number to that without the null char. then gen filename.
             * set the number back to 01 each time by redefining fn,
             * otherwise we'd have to keep track of the previous number
             * and search for that. i guess alternatively we could
             * hardcode the index of the underscore and overwrite
             * two characters after that. but who cares. it's ugly either way
             * */
            char fn[200] = "out/NLLZ/7_PROD_2/2/1000/car_01.dat\0";
            char snum[3], car_num[2], num[2];
            int status = snprintf(snum, 3, "%02u", ii);
            memcpy(car_num, snum, 2);
            status = generate_filename(sizeof(fn), fn, "01", snum);
            if (status != 0) {
              fprintf(stdout, "filename generation for outputting"
                  "car lineshape %d didn't work\n", ii);
              break;
            } else {
              FILE *fp = fopen(fn, "w");
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
            if (ii == vib_total) {
              char fn[200] = "out/NLLZ/7_PROD_2/2/1000/chl_01.dat\0";
              status = snprintf(snum, 3, "%02u", chl_index);
              memcpy(num, snum, 2);
              fprintf(stdout, "num = %s", num);
              status = generate_filename(sizeof(fn), fn, "01", snum);
              fprintf(stdout, "%s", fn);
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

          } // car = n_chl
        }

        for (unsigned step = 0; step < tau; step++) {
          fi_ad[step] = normed_fi[chl_index][step] * abs[step];
          ai_fd[step] = normed_ai[chl_index][step] * abs[step];
        }
        
        if (print_delta_fc) {
          fprintf(stdout, "%5u %12.8e %12.8e %12.8e\n", ii,
              e_xa, fc_sq, ji_work);
        }

        chl_car += pow(ji_work, 2.) * trapezoid(fi_ad, tau);
        car_chl += pow(ji_work, 2.) * trapezoid(ai_fd, tau);

        /* only calculate this once lol */
        if (chl_index == 0 && carotenoid == n_chl) {
          for (unsigned kk = 0; kk < vib_total; kk++) {

            /* NB: indexing!!! */
            car_decays[ii - vib_total] += car_rates[sub2ind({kk, ii},
                            {n_s_total, n_s_total})];
            if (print_decay_details) {
              fprintf(stdout, "%2u %2u %10.6e\n",
                  ii, kk, car_rates[sub2ind({kk, ii},
                  {n_s_total, n_s_total})]);
            }

          } // kk
        }

        chl_car *= CM_PER_PS * 2 * PI;
        car_chl *= CM_PER_PS * 2 * PI;
        if (eigvals[chl_index] > e_xa) {
          car_chl *= exp(-beta * (eigvals[chl_index] - e_xa));
        } else {
          chl_car *= exp(-beta * (e_xa - eigvals[chl_index]));
        }
        k_i_xa.push_back(chl_car);
        k_i_xa.push_back(car_chl);

        if (print_details) {
          fprintf(stdout, "chl = %2u, car = %2u,"
          " state = (%2lu, %2lu, %2lu): e_chl = %10.6e,"
          " e_car = %10.6e, fc^2 = %10.6e, forward rate = %10.6e,"
          " backward_rate = %10.6e\n", chl_index, carotenoid,
          a[0], a[1], a[2], eigvals[chl_index], e_xa, fc_sq, chl_car, car_chl);
        }

      } // ii
    }
  } // chl_index
  free(abs);
  free(fi_ad);
  free(ai_fd);

  return k_i_xa;
}

double**
car_transfer(VERA x, double *decays)
{
  std::vector<size_t> extents = x.get_pop_extents();
  std::vector<size_t> subscripts(extents.size(), 0.),
                      j_subs(extents.size(), 0.);

  unsigned s1_total = pow(x.n_vib + 1, x.n_normal) + 1;
  unsigned n_total = x.n_total;
  std::vector<double> rates = x.intra_rates();

  double **t = (double **)calloc(n_total, sizeof(double*));
  for (unsigned i = 0; i < n_total; i++) {
    t[i] = (double *)calloc(n_total, sizeof(double));
  }

  /* indices - ground state is index 0 */
  for (unsigned i = 1; i <= s1_total; i++) {
    subscripts = ind2sub(i - 1, extents);  
    /* this will return 0, a_1, a_2, so add one to the electronic index */
    subscripts[0]++;
    size_t ind = sub2ind(subscripts, extents);

    for (unsigned j = 1; j <= s1_total; j++) {
      if (i == j) {
        t[i - 1][j - 1] -= decays[i - 1];
        t[0][j - 1] += decays[i - 1];
      } else {
        j_subs = ind2sub(j - 1, extents);

        j_subs[0]++;
        size_t j_ind = sub2ind(j_subs, extents);
        t[j - 1][i - 1] = rates[sub2ind({ind, j_ind}, {n_total, n_total})];
      }

    }

  }
  return t;
}

double**
hybrid_transfer(unsigned n_chl, unsigned n_car, VERA *x,
    double *gamma, double **Jij, std::vector<double> k_i_delta,
    double **redfield_rates, double *car_decays)
{
  /* n_vib_tot is the total number of vibrational levels
   * per electronic state; we only need this here because
   * we're not concerned with hot ground states/S2 */
  unsigned n_vib_tot = pow(x->n_vib + 1, x->n_normal);
  unsigned size = n_chl + 1 + n_car * (n_vib_tot + 1);
  double **k_tot = (double **)calloc(size, sizeof(double*));
  for (unsigned i = 0; i < size; i++) {
    k_tot[i] = (double *)calloc(size, sizeof(double));
  }

  std::vector<size_t> car_extents = x->get_pop_extents();
  std::vector<size_t> chl_car_extents = {n_chl, n_car, n_vib_tot, 2};
  std::vector<double> car_rates = x->intra_rates();

  for (unsigned i = 0; i < size; i++) {
    for (unsigned j = 0; j < size; j++) {
      unsigned rgs = 0;
      unsigned gs_620 = n_chl + 1;
      unsigned gs_621 = n_vib_tot + 1 + gs_620;
      bool i_rgs  = (i == rgs);
      bool j_rgs  = (j == rgs);
      bool i_chls = (i > rgs && i <= n_chl);
      bool j_chls = (j > rgs && j <= n_chl);
      bool i_620_gs  = (i == gs_620);
      bool j_620_gs  = (j == gs_620);
      bool i_620  = (i > gs_620 && i < gs_621);
      bool j_620  = (j > gs_620 && j < gs_621);
      bool i_621_gs  = (i == gs_621);
      bool j_621_gs  = (j == gs_621);
      bool i_621  = (i > gs_621);
      bool j_621  = (j > gs_621);

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
            for (unsigned k = 0; k < n_vib_tot; k++) {
              k_tot[i][j] -= k_i_delta[sub2ind({i - 1, 0, k, 0},
                             chl_car_extents)];
              k_tot[i][j] -= k_i_delta[sub2ind({i - 1, 1, k, 0},
                             chl_car_extents)];
            }
          }
        }
        if (j_620) {
          k_tot[j][i] = k_i_delta[sub2ind({i - 1, 0, j - (gs_620 + 1), 0},
                        chl_car_extents)];
        }
        if (j_621) {
          k_tot[j][i] = k_i_delta[sub2ind({i - 1, 1,
              j - (gs_621 + 1), 0}, chl_car_extents)];
        }
      } // i_chls

      if (i_620_gs) {
        if (j_rgs) {
        }
        if (j_chls) {
        }
        if (j_620) {
          k_tot[i][j] += car_decays[j - (gs_620 + 1)];
        }
        if (j_621) {
        }
      } // i_620_gs

      if (i_620) {
        if (j_rgs) {
        }
        if (j_chls) {
          k_tot[j][i] += k_i_delta[sub2ind({j - 1, 0,
              i - (gs_620 + 1), 1}, chl_car_extents)];
        }
        if (j_620) {
          if (i == j) {
            k_tot[i][j] -= car_decays[j - (gs_620 + 1)];
            for (unsigned k = 0; k < n_chl; k++) {
              k_tot[i][j] -= k_i_delta[sub2ind({k, 0, j - (gs_620 + 1), 1},
                            chl_car_extents)];
            }
          }
          /* this needs confirming, but:
           * the total rate matrix, VERA.intra_rates(), contains
           * both the IVR and IC rates. but on a given electronic
           * state, only IVR is nonzero, and across electronic
           * states, only IC is nonzero. hence, if we ensure we're 
           * on the same electronic state, intra_rates()[i][j]
           * should just be the IVR rate.
           */
          std::vector<size_t> i_subs = ind2sub(i - (gs_620 + 1), car_extents);
          std::vector<size_t> j_subs = ind2sub(j - (gs_620 + 1), car_extents);
          /* S1! */
          i_subs[0]++; j_subs[0]++;
          size_t i_index = sub2ind(i_subs, car_extents);
          size_t j_index = sub2ind(j_subs, car_extents);
          k_tot[i][j] = car_rates[sub2ind({i_index, j_index},
              {x->n_total, x->n_total})];
        }
        if (j_621) {
        }
      }

      /* this can be folded into the 620 bit i guess; these blocks
       * will be identical, it's only the chl-car blocks which differ */

      if (i_621_gs) {
        if (j_rgs) {
        }
        if (j_chls) {
        }
        if (j_620) {
        }
        if (j_621) {
          k_tot[i][j] += car_decays[j - (gs_621 + 1)];
        }
      } // i_620_gs
      if (i_621) {
        if (j_rgs) {
        }
        if (j_chls) {
          k_tot[j][i] += k_i_delta[sub2ind({j - 1, 1,
              i - (gs_621 + 1), 1}, chl_car_extents)];
        }
        if (j_620) {
        }
        if (j_621) {
          if (i == j) {
            k_tot[i][j] -= car_decays[j - (gs_621 + 1)];
            for (unsigned k = 0; k < n_chl; k++) {
              k_tot[i][j] -= k_i_delta[sub2ind({k, 1,
                  j - (gs_621 + 1), 1}, chl_car_extents)];
            }
          }
          std::vector<size_t> i_subs = ind2sub(i - (gs_621 + 1),
                                       car_extents);
          std::vector<size_t> j_subs = ind2sub(j - (gs_621 + 1),
                                       car_extents);
          i_subs[0]++; j_subs[0]++;
          size_t i_index = sub2ind(i_subs, car_extents);
          size_t j_index = sub2ind(j_subs, car_extents);
          k_tot[i][j] = car_rates[sub2ind({i_index, j_index},
              {x->n_total, x->n_total})];
        }
      }

    } // j
  }
  return k_tot;

}

double*
hybrid_boltz(unsigned n_chl, unsigned n_car, double beta,
             double *eigvals, VERA *x)
{
  /* get the boltzmann factor for every state and normalise, 
   * then return just the chlorophyll boltzmann factors. we do this
   * because the carotenoid doesn't fluoresce but does have some
   * population on it. note that, as everywhere else in the code,
   * this assumes all the carotenoids are identical */
  unsigned n_vib_tot = pow(x->n_vib + 1, x->n_normal);
  double *work = (double *)calloc(n_chl + n_car * n_vib_tot, sizeof(double));
  double *boltz = (double *)calloc(n_chl, sizeof(double));
  std::vector<size_t> extents = x->get_pop_extents();
  double work_sum = 0.;

  for (unsigned i = 0; i < n_chl + n_vib_tot; i++) {
    if (i < n_chl) {
      work[i] = exp(-beta * eigvals[i]); // dimensions!
      work_sum += work[i];
    } else {
      std::vector<size_t> subs = ind2sub(i - n_chl, extents);
      subs[0]++;
      double e_xa = x->get_w_elec(subs[0]);
      for (unsigned alpha = 0; alpha < x->n_normal; alpha++) {
        e_xa += x->get_w_normal(alpha) * subs[alpha + 1];
      }
      work[i] = exp(-beta * e_xa);
      for (unsigned j = 0; j < n_car; j++) {
        work[i + (j * n_vib_tot)] = exp(-beta * e_xa);
      }
      work_sum += n_car * work[i];
    }
  }
  for (unsigned i = 0; i < n_chl + n_car * n_vib_tot; i++) {
    work[i] /= work_sum;
  }
  for (unsigned i = 0; i < n_chl; i++) {
    boltz[i] = work[i];
  }
  free(work);
  return boltz;

}
