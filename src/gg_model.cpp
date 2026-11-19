// Gibbs sampler for the Gaussian precision factor model of
// Chandra, Mueller and Sarkar (2026), Bayesian Analysis 21(1), 75-103,
// doi:10.1214/24-BA1461.
//
// Original implementation: Noirrit Kiran Chandra (PrecFact 1.0, 2024).
// This version keeps the algorithm and removes the system dependencies of the original
// (GSL, mlpack, RcppDist, RcppEnsmallen, OpenMP) so that the package builds with only
// Rcpp and RcppArmadillo on Windows, macOS and Linux. All random numbers come from R's
// generator, so set.seed() makes runs reproducible.
//
// The generalized inverse Gaussian generator (rgig and helpers) is adapted from the
// GIGrvg package by Josef Leydold and Wolfgang Hormann (GPL >= 2).

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// ---------------------------------------------------------------------------------------
// Small utilities
// ---------------------------------------------------------------------------------------

// [[Rcpp::export]]
double log_sum_exp(const arma::vec &x) {
  double maxVal = x.max();
  double sum_exp = sum(exp(x - maxVal));
  return std::log(sum_exp) + maxVal;
}

// log of the marginal density (up to the constant -n/2 log(2 pi)) of an n-vector v with
// v | delta2 ~ N(0, I / delta2) and delta2 ~ Gamma(as, rate = bs); vnorm = ||v||^2.
// lgamma(as + n/2) - lgamma(as) replaces gsl_sf_lnpoch(as, n/2) of the original.
inline double log_marginal_pdf(unsigned n, double as, double bs, double vnorm) {
  const double h = ((double) n) / 2.0;
  return R::lgammafn(as + h) - R::lgammafn(as) + as * std::log(bs) -
         (as + h) * std::log(bs + vnorm / 2.0);
}

// Inverse Gaussian with mean mu and shape 1 (Michael, Schucany and Haas 1976).
double rig(double mu) {
  double y = R::norm_rand();
  y *= y;
  const double mu2 = mu * mu;
  const double quad = 4 * mu * y + mu2 * y * y;
  const double x_div_mu = (2 + mu * y - std::sqrt(quad)) / 2;
  const double x = mu * x_div_mu;
  if (x <= 0) return 1e-5;
  const double u = std::log(R::unif_rand());
  if (u <= -std::log1p(x_div_mu)) return x;
  return mu / x_div_mu;
}

// ---------------------------------------------------------------------------------------
// Generalized inverse Gaussian generator, adapted from GIGrvg (Leydold and Hormann).
// Density proportional to x^(lambda-1) exp(-(chi/x + psi x)/2).
// ---------------------------------------------------------------------------------------

static double gig_mode(double lambda, double omega) {
  if (lambda >= 1.)
    return (std::sqrt((lambda - 1.) * (lambda - 1.) + omega * omega) + (lambda - 1.)) / omega;
  return omega / (std::sqrt((1. - lambda) * (1. - lambda) + omega * omega) + (1. - lambda));
}

static double rgig_ROU_shift_alt(double lambda, double lambda_old, double omega, double alpha) {
  const double t = 0.5 * (lambda - 1.);
  const double s = 0.25 * omega;
  const double xm = gig_mode(lambda, omega);
  const double nc = t * std::log(xm) - s * (xm + 1. / xm);
  const double a = -(2. * (lambda + 1.) / omega + xm);
  const double b = (2. * (lambda - 1.) * xm / omega - 1.);
  const double c = xm;
  const double p = b - a * a / 3.;
  const double q = (2. * a * a * a) / 27. - (a * b) / 3. + c;
  const double fi = std::acos(-q / (2. * std::sqrt(-(p * p * p) / 27.)));
  const double fak = 2. * std::sqrt(-p / 3.);
  const double y1 = fak * std::cos(fi / 3.) - a / 3.;
  const double y2 = fak * std::cos(fi / 3. + 4. / 3. * M_PI) - a / 3.;
  const double uplus = (y1 - xm) * std::exp(t * std::log(y1) - s * (y1 + 1. / y1) - nc);
  const double uminus = (y2 - xm) * std::exp(t * std::log(y2) - s * (y2 + 1. / y2) - nc);
  double U, V, X;
  do {
    U = uminus + R::unif_rand() * (uplus - uminus);
    V = R::unif_rand();
    X = U / V + xm;
  } while ((X <= 0.) || (std::log(V) > (t * std::log(X) - s * (X + 1. / X) - nc)));
  return (lambda_old < 0.) ? (alpha / X) : (alpha * X);
}

static double rgig_newapproach1(double lambda, double lambda_old, double omega, double alpha) {
  if (lambda >= 1. || omega > 1.) Rcpp::stop("invalid parameters for the GIG generator");
  double A[3];
  const double xm = gig_mode(lambda, omega);
  const double x0 = omega / (1. - lambda);
  const double k0 = std::exp((lambda - 1.) * std::log(xm) - 0.5 * omega * (xm + 1. / xm));
  double k1, k2;
  A[0] = k0 * x0;
  if (x0 >= 2. / omega) {
    k1 = 0.;
    A[1] = 0.;
    k2 = std::pow(x0, lambda - 1.);
    A[2] = k2 * 2. * std::exp(-omega * x0 / 2.) / omega;
  } else {
    k1 = std::exp(-omega);
    A[1] = (lambda == 0.) ? k1 * std::log(2. / (omega * omega))
                          : k1 / lambda * (std::pow(2. / omega, lambda) - std::pow(x0, lambda));
    k2 = std::pow(2 / omega, lambda - 1.);
    A[2] = k2 * 2 * std::exp(-1.) / omega;
  }
  const double Atot = A[0] + A[1] + A[2];
  double U, V, X, hx, a;
  while (true) {
    V = Atot * R::unif_rand();
    do {
      if (V <= A[0]) {
        X = x0 * V / A[0];
        hx = k0;
        break;
      }
      V -= A[0];
      if (V <= A[1]) {
        if (lambda == 0.) {
          X = omega * std::exp(std::exp(omega) * V);
          hx = k1 / X;
        } else {
          X = std::pow(std::pow(x0, lambda) + (lambda / k1 * V), 1. / lambda);
          hx = k1 * std::pow(X, lambda - 1.);
        }
        break;
      }
      V -= A[1];
      a = (x0 > 2. / omega) ? x0 : 2. / omega;
      X = -2. / omega * std::log(std::exp(-omega / 2. * a) - omega / (2. * k2) * V);
      hx = k2 * std::exp(-omega / 2. * X);
    } while (0);
    U = R::unif_rand() * hx;
    if (std::log(U) <= (lambda - 1.) * std::log(X) - omega / 2. * (X + 1. / X))
      return (lambda_old < 0.) ? (alpha / X) : (alpha * X);
  }
}

static double rgig_ROU_noshift(double lambda, double lambda_old, double omega, double alpha) {
  const double t = 0.5 * (lambda - 1.);
  const double s = 0.25 * omega;
  const double xm = gig_mode(lambda, omega);
  const double nc = t * std::log(xm) - s * (xm + 1. / xm);
  const double ym = ((lambda + 1.) + std::sqrt((lambda + 1.) * (lambda + 1.) + omega * omega)) / omega;
  const double um = std::exp(0.5 * (lambda + 1.) * std::log(ym) - s * (ym + 1. / ym) - nc);
  double U, V, X;
  do {
    U = um * R::unif_rand();
    V = R::unif_rand();
    X = U / V;
  } while (std::log(V) > (t * std::log(X) - s * (X + 1. / X) - nc));
  return (lambda_old < 0.) ? (alpha / X) : (alpha * X);
}

#define ZTOL (DBL_EPSILON * 10.0)

double rgig(double lambda, double chi, double psi) {
  if (!(R_FINITE(lambda) && R_FINITE(chi) && R_FINITE(psi)) || (chi < 0. || psi < 0) ||
      (chi == 0. && lambda <= 0.) || (psi == 0. && lambda >= 0.)) {
    Rcpp::stop("invalid parameters for the GIG distribution: lambda = %f, chi = %f, psi = %f",
               lambda, chi, psi);
  }
  if (chi < ZTOL) {
    if (lambda > 0.0) return R::rgamma(lambda, 2.0 / psi);
    return 1.0 / R::rgamma(-lambda, 2.0 / psi);
  }
  if (psi < ZTOL) {
    if (lambda > 0.0) return 1.0 / R::rgamma(lambda, 2.0 / chi);
    return R::rgamma(-lambda, 2.0 / chi);
  }
  const double lambda_old = lambda;
  if (lambda < 0.) lambda = -lambda;
  const double alpha = std::sqrt(chi / psi);
  const double omega = std::sqrt(psi * chi);
  if (lambda > 2. || omega > 3.) return rgig_ROU_shift_alt(lambda, lambda_old, omega, alpha);
  if (lambda >= 1. - 2.25 * omega * omega || omega > 0.2)
    return rgig_ROU_noshift(lambda, lambda_old, omega, alpha);
  if (lambda >= 0. && omega > 0.) return rgig_newapproach1(lambda, lambda_old, omega, alpha);
  Rcpp::stop("parameters must satisfy lambda >= 0 and omega > 0");
  return 0.0;
}

// ---------------------------------------------------------------------------------------
// Fast inverse and log-determinant of Omega = Lambda Lambda' + diag(ps)
// ---------------------------------------------------------------------------------------

// [[Rcpp::export]]
arma::mat fast_fact_inv(const arma::mat &Lambda, const arma::vec &ps) {
  arma::mat Lmsg = Lambda.each_col() / ps;
  arma::mat Lmsg2 = symmatu(Lmsg.t() * Lambda);
  Lmsg2.diag() += 1;
  arma::mat S = chol(Lmsg2, "lower");
  arma::mat exp_mat2 = solve(trimatl(S), Lmsg.t());
  return diagmat(1 / ps) - symmatu(exp_mat2.t() * exp_mat2);
}

// [[Rcpp::export]]
double fast_fact_det(const arma::mat &Lambda, const arma::vec &ps) {
  arma::vec s;
  arma::vec sq_ps = sqrt(ps);
  arma::mat Lambda_orig = Lambda.each_col() / sq_ps;
  svd(s, Lambda_orig);
  return sum(log1p(square(s))) + sum(log(ps));
}

// ---------------------------------------------------------------------------------------
// Gibbs steps
// ---------------------------------------------------------------------------------------

// Step 4 of the sampler: Dirichlet-Laplace(a, b) scales. T_phi holds tau * phi.
static void update_loading_params(const arma::mat &Lambda, arma::mat &T_phi, arma::mat &Plam,
                                  const double a, const double b) {
  arma::mat abs_theta = abs(Lambda);
  arma::mat psi_inv = (T_phi / abs_theta);
  psi_inv.transform([](double val) { return rig(val); });
  T_phi = abs_theta;
  T_phi.transform([&a, &b](double val) {
    double tmp = rgig(a - 1, 2 * val, 2 * b);
    return std::max(tmp, 9e-6);
  });
  if (T_phi.has_nan()) Rcpp::stop("NaN in the Dirichlet-Laplace scales");
  Plam = (psi_inv / square(T_phi));
  if (Plam.has_nan()) Rcpp::stop("NaN in the prior precisions of the loadings");
}

// Step 2 of the sampler: rows of Lambda given the latent u and v.
static void update_phi(arma::mat &phimat, const arma::mat &Plam, const arma::mat &u,
                       const arma::mat &v, const arma::vec &ssqv) {
  const unsigned p = phimat.n_rows;
  arma::mat v_phi = (v * phimat);
  for (unsigned j = 0; j < p; ++j) {
    arma::mat tmp = v.col(j) * phimat.row(j);
    arma::mat util = u - v_phi + tmp;
    rowvec w = (util.t() * v.col(j)).t();
    rowvec d = 1 / sqrt(Plam.row(j) + ssqv(j));
    phimat.row(j) = d % (d % w + rowvec(phimat.n_cols, fill::randn));
    v_phi -= (tmp - v.col(j) * phimat.row(j));
  }
}

// Step 5: Dirichlet process concentration (Escobar and West 1995).
static double update_dir_prec(double alpha, double a, double b, unsigned n, unsigned k) {
  const double phi = R::rbeta(alpha + 1, (double) n);
  const double gam1 = a + k, gam2 = b - std::log(phi);
  const double pi = (gam1 - 1) / (n * gam2);
  return (std::log(R::unif_rand()) < std::log(pi) - std::log1p(pi))
             ? R::rgamma(gam1, 1 / gam2)
             : R::rgamma(gam1 - 1, 1 / gam2);
}

// Two-cluster 1-D k-means used only to initialise the Dirichlet process allocation.
// Replaces the mlpack::KMeans call of the original, which also used min(ceil(log p), 2)
// clusters on the column sums of squares. Starts from a split at the median.
static uvec kmeans_1d_init(const arma::vec &x, unsigned k) {
  const unsigned p = x.n_elem;
  uvec lab(p, fill::zeros);
  if (k < 2 || p < 2) return lab;
  const double med = median(x);
  for (unsigned j = 0; j < p; ++j) lab(j) = (x(j) > med) ? 1 : 0;
  for (int it = 0; it < 100; ++it) {
    uvec i0 = find(lab == 0), i1 = find(lab == 1);
    if (i0.n_elem == 0 || i1.n_elem == 0) break;
    const double m0 = mean(x(i0)), m1 = mean(x(i1));
    uvec newlab(p);
    for (unsigned j = 0; j < p; ++j)
      newlab(j) = (std::abs(x(j) - m1) < std::abs(x(j) - m0)) ? 1 : 0;
    if (all(newlab == lab)) break;
    lab = newlab;
  }
  if (all(lab == 1)) lab.zeros();
  return lab;
}

// [[Rcpp::export]]
Rcpp::List prec_est_cpp(const double as, const double bs, const double a, const double b,
                        const int nrun, const int thin, const arma::mat &lambda_init,
                        const arma::mat &Y, const bool alldiff, const double a_dir,
                        const double b_dir, const bool verbose) {
  mat phimat = lambda_init;
  if (Y.n_cols != phimat.n_rows)
    Rcpp::stop("Dimension mismatch between Y and the initial value of lambda");
  unsigned j;
  const unsigned n = Y.n_rows, p = Y.n_cols, k = phimat.n_cols;

  // Dirichlet process on the residual precisions
  double dir_prec = 1.0;
  const unsigned nmix = p;
  double prob_empty = std::log(dir_prec);
  uvec nj_val(nmix, fill::zeros), del;
  arma::vec log_nj_val(nmix, fill::zeros), probs(nmix), log_probs(nmix);
  std::vector<uvec> inds_eq_j(nmix);
  arma::vec ps(p, fill::ones);

  arma::mat T_phi = randg<arma::mat>(p, k, distr_param(a, 0.5));
  T_phi.transform([](double val) { return std::max(val, 9e-17); });
  arma::mat Plam = 1 / (1 + square(T_phi));

  const unsigned n_mc = std::floor((double) nrun / thin);
  arma::mat ps_mat(n_mc, p);
  arma::vec n_clusters(n_mc, fill::zeros);
  cube Phi_mc(p, k, n_mc, fill::zeros);

  const int report_every = std::max(1, nrun / 10);
  for (int i = 0; i < nrun; ++i) {
    if (i % 50 == 0) Rcpp::checkUserInterrupt();
    if (verbose && (i + 1) % report_every == 0)
      Rcpp::Rcout << "PrecFact: iteration " << (i + 1) << " of " << nrun << std::endl;

    // Step 1: u and v
    arma::mat Lmsg = phimat.each_col() / ps;
    arma::mat Lmsg2 = symmatu(Lmsg.t() * phimat);
    Lmsg2.diag() += 1.0;
    arma::mat noise(n, k, fill::randn);
    arma::mat chol_precmat_u = chol(Lmsg2);
    arma::mat u = noise * trimatu(chol_precmat_u);
    arma::mat tmpmt = solve(trimatl(chol_precmat_u.t()), Lmsg.t());
    arma::mat v = Y + noise * tmpmt;
    arma::vec ssqv = sum(square(v), 0).t();

    // Step 2: loadings
    update_phi(phimat, Plam, u, v, ssqv);

    // Step 4: Dirichlet-Laplace parameters
    update_loading_params(phimat, T_phi, Plam, a, b);

    // Step 3: residual precisions delta_j^2
    unsigned counttttt = p;
    if (alldiff) {
      ps = bs + 0.5 * ssqv;
      ps.transform([&as, &n](double val) { return R::rgamma(as + 0.5 * n, 1 / val); });
    } else {
      if (!i) {
        unsigned n_init = std::min((unsigned) std::ceil(std::log((double) p)), 2u);
        del = kmeans_1d_init(ssqv, n_init);
        for (j = 0; j < nmix; ++j) {
          inds_eq_j[j] = find(del == j);
          nj_val(j) = inds_eq_j[j].n_elem;
          log_nj_val(j) = std::log((double) nj_val(j));
        }
      }
      arma::vec v_atoms(nmix, fill::zeros);
      for (j = 0; j < nmix; ++j)
        if (nj_val(j)) v_atoms(j) = sum(ssqv(inds_eq_j[j]));

      for (unsigned jj = 0; jj < p; ++jj) {
        double dens, cluster_prob;
        const double pdf_empty = log_marginal_pdf(n, as, bs, ssqv(jj));
        uvec non_empty_clusts = find(nj_val > 0);
        const uword current_ind = del(jj);
        if (nj_val(current_ind) > 1) {
          uvec empty_inds = find(nj_val == 0, 1, "first");
          non_empty_clusts.resize(non_empty_clusts.n_elem + 1);
          non_empty_clusts.tail(1) = empty_inds(0);
        }
        const unsigned nclust = non_empty_clusts.n_elem;
        log_probs.set_size(nclust);
        for (unsigned j1 = 0; j1 < nclust; ++j1) {
          j = non_empty_clusts(j1);
          if (j != current_ind) {
            if (nj_val(j) < 1) {
              dens = pdf_empty;
              cluster_prob = prob_empty;
            } else {
              dens = log_marginal_pdf(n, as + ((double) (n * nj_val(j))) / 2,
                                      bs + v_atoms(j) / 2, ssqv(jj));
              cluster_prob = log_nj_val(j);
            }
            log_probs(j1) = dens + cluster_prob;
          }
        }
        j = current_ind;
        if (nj_val(j) == 1) {
          dens = pdf_empty;
          cluster_prob = prob_empty;
        } else {
          dens = log_marginal_pdf(n, as + ((double) (n * (nj_val(j) - 1))) / 2,
                                  bs + (v_atoms(j) - ssqv(jj)) / 2, ssqv(jj));
          cluster_prob = std::log((double) (nj_val(j) - 1));
        }
        uvec indddd = find(non_empty_clusts == j, 1, "first");
        log_probs(indddd(0)) = dens + cluster_prob;

        const double log_DEN = log_sum_exp(log_probs);
        probs = normalise(exp(log_probs - log_DEN), 1);
        if (std::abs(sum(probs) - 1.0) > 1e-5) {
          log_probs -= max(log_probs);
          probs = normalise(exp(log_probs), 1);
        }
        if (sum(probs) == 0)
          Rcpp::stop("allocation probabilities sum to zero at variable %d, iteration %d", jj + 1, i + 1);

        // draw the new cluster
        const double uu = R::unif_rand();
        double cum = 0.0;
        unsigned pick = nclust - 1;
        for (unsigned j1 = 0; j1 < nclust; ++j1) {
          cum += probs(j1);
          if (uu <= cum) {
            pick = j1;
            break;
          }
        }
        del(jj) = non_empty_clusts(pick);

        if (del(jj) != current_ind) {
          j = current_ind;
          --nj_val(j);
          if (nj_val(j)) {
            log_nj_val(j) = std::log((double) nj_val(j));
            v_atoms(j) -= ssqv(jj);
          }
          j = del(jj);
          log_nj_val(j) = std::log((double) (++nj_val(j)));
          if (nj_val(j) == 1)
            v_atoms(j) = ssqv(jj);
          else
            v_atoms(j) += ssqv(jj);
        }
      }
      counttttt = 0;
      for (j = 0; j < nmix; ++j) {
        if (nj_val(j)) {
          ++counttttt;
          inds_eq_j[j] = find(del == j);
          const double atom =
              R::rgamma(as + 0.5 * (n * nj_val(j)), 1 / (bs + .5 * v_atoms(j)));
          (ps(inds_eq_j[j])).fill(atom);
        } else
          inds_eq_j[j].reset();
        if (inds_eq_j[j].n_elem != nj_val(j)) Rcpp::stop("cluster occupancy mismatch");
      }
      // Step 5
      dir_prec = update_dir_prec(dir_prec, a_dir, b_dir, p, counttttt);
      prob_empty = std::log(dir_prec);
    }

    if ((i + 1) % thin == 0) {
      const unsigned s = (i + 1) / thin - 1;
      ps_mat.row(s) = ps.t();
      Phi_mc.slice(s) = phimat;
      n_clusters(s) = counttttt;
    }
  }

  return Rcpp::List::create(Rcpp::Named("Lambda") = Phi_mc,
                            Rcpp::Named("residuals") = ps_mat,
                            Rcpp::Named("n_clusters") = n_clusters);
}

// Posterior summaries over the saved draws, computed one draw at a time so that the
// d x d x S array of the original post-processing is never formed.
// Returns the posterior mean of Omega, the posterior mean of the partial correlation
// matrix, and, for each epsilon in the grid, the posterior probability that
// |partial correlation| > epsilon.
// [[Rcpp::export]]
Rcpp::List posterior_summaries_cpp(const arma::cube &Lambda, const arma::mat &residuals,
                                   const arma::uvec &keep, const arma::vec &eps_grid) {
  const unsigned p = Lambda.n_rows, S = keep.n_elem, G = eps_grid.n_elem;
  arma::mat Omega_mean(p, p, fill::zeros), Pcor_mean(p, p, fill::zeros);
  arma::mat Pcor_sq(p, p, fill::zeros);
  arma::cube counts(p, p, G, fill::zeros);
  arma::vec logdet(S);
  for (unsigned s = 0; s < S; ++s) {
    if (s % 20 == 0) Rcpp::checkUserInterrupt();
    const unsigned idx = keep(s);
    arma::mat L = Lambda.slice(idx);
    arma::vec ps = residuals.row(idx).t();
    logdet(s) = fast_fact_det(L, ps);
    arma::mat Om = L * L.t();
    Om.diag() += ps;
    Omega_mean += Om;
    arma::vec isd = 1 / sqrt(Om.diag());
    arma::mat R = Om.each_col() % isd;
    R.each_row() %= isd.t();
    arma::mat Pc = -R;
    Pc.diag().ones();
    Pcor_mean += Pc;
    Pcor_sq += square(Pc);
    arma::mat A = abs(R);
    for (unsigned g = 0; g < G; ++g) counts.slice(g) += conv_to<mat>::from(A > eps_grid(g));
  }
  Omega_mean /= S;
  Pcor_mean /= S;
  Pcor_sq /= S;
  arma::mat Pcor_sd = sqrt(clamp(Pcor_sq - square(Pcor_mean), 0.0, datum::inf));
  counts /= S;
  return Rcpp::List::create(Rcpp::Named("Omega") = Omega_mean,
                            Rcpp::Named("pcor") = Pcor_mean,
                            Rcpp::Named("pcor_sd") = Pcor_sd,
                            Rcpp::Named("postprob") = counts,
                            Rcpp::Named("logdet") = logdet);
}
