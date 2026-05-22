// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// ---------------------------------------------------------------------------
// Diagonal preconditioning conventions used throughout this file
// ---------------------------------------------------------------------------
// USER-FACING CONVENTION:
//   The `precond` argument to *every* sampler is the diagonal of an estimate
//   of the posterior COVARIANCE (i.e. target marginal variances).  Larger
//   entries -> looser posterior in that coordinate.  This is the natural
//   quantity to estimate from a pilot run via apply(samples, 2, var).
//
// Stored as an arma::vec `v` (length p), v_j > 0.
//
//   * RWM:   proposal  w' = w + h * sqrt(v) % randn(p),
//            so proposal covariance is h^2 * diag(v).
//
//   * HMC:   the mass matrix is M = diag(v)^{-1} (i.e. mass = 1/variance,
//            so heavy directions move slowly, light directions move fast).
//            Equivalently, in code we use d = 1/v and:
//              momentum p ~ N(0, M):           p = sqrt(d) % randn(p)
//                                              = randn(p) / sqrt(v).
//              position step qdot = M^{-1} p = p % v.
//              kinetic energy   K = 0.5 * p^T M^{-1} p = 0.5 * sum(p^2 * v).
//            This is the standard HMC convention (Neal 2011, eq. 4.20+):
//            picking M ~ Cov(posterior)^{-1} makes the trajectory cover
//            ellipsoids matched to the posterior scale.
//
//   * Proximal/Guo HMC: same mass-matrix HMC; the proximal step is kept
//     in the original (unscaled) parameterisation so it stays a cheap
//     elementwise soft-threshold on w.
//
// All "% / sqrt / square" operations on `v` are O(p), so diagonal
// preconditioning adds essentially no cost.
// ---------------------------------------------------------------------------


// Gradient of f for Durmus
// [[Rcpp::export]]
arma::vec gradf_dur(const arma::vec& w,
                    const arma::vec& y,
                    const arma::mat& B,
                    double nu,
                    double sigma)
{
  arma::vec r = y - B * w;
  
  arma::vec weights =
    r / (nu * sigma * sigma + r % r);
  
  arma::vec grad =
    - (nu + 1.0) * B.t() * weights;
    
    return grad;
}

// log target
// [[Rcpp::export]]
double log_pi(const arma::vec& w,
              const arma::vec& y,
              const arma::mat& B,
              double nu,
              double alpha,
              double sigma)
{
  arma::vec r = y - B * w;
  
  double log_term = 0.5 * (nu + 1.0) *
    sum(log(1.0 + square(r) / (nu * sigma * sigma)));
  
  // Optional normalization term
  log_term += y.n_elem * std::log(sigma);
  
  double penalty = alpha * norm(w, 1);
  
  return -(log_term + penalty);
}

// soft thresholding function
arma::vec soft_threshold(const arma::vec& z, double threshold) {
  return arma::sign(z) % arma::max(arma::abs(z) - threshold,
                    arma::zeros(z.n_elem));
}

// proximal mapping for phmc
arma::vec prox_phmc(const arma::vec& w, double lambda, double alpha) {
  return soft_threshold(w, lambda * alpha);
}

// gradient of the partial proximal methods
// [[Rcpp::export]]
arma::vec grad_logpiLamg(const arma::mat& B,
                         const arma::vec& y,
                         const arma::vec& w,
                         double lambda,
                         double alpha,
                         double nu,
                         double sigma)
{
  arma::vec grad_f = gradf_dur(w, y, B, nu, sigma);
  
  arma::vec w_prox = prox_phmc(w, lambda, alpha);
  
  return -(grad_f + (w - w_prox) / lambda);
}

// gradient for the guo method
// [[Rcpp::export]]
arma::vec grad_logpiLam_guo(const arma::vec& w, double lambda, double alpha) {
  arma::vec w_prox = prox_phmc(w, lambda, alpha);
  return - (w - w_prox);
}

// ---------------------------------------------------------------------------
// Default diagonal preconditioner (covariance scale)
// ---------------------------------------------------------------------------
// Per the convention above, `precond` is the posterior-VARIANCE diagonal.
// A cheap default is the inverse of the expected Fisher diagonal:
//          v_j  =  (nu+3) * sigma^2 / ( (nu+1) * ||B_{.j}||^2 ).
// (i.e. 1 / H_{jj}, where H_{jj} ~ ((nu+1)/(nu+3)) * ||B_{.j}||^2 / sigma^2
// is the diagonal of the expected Hessian of the smooth part of the
// negative log-target).  We cap large values to a ceiling to guard
// against (near-)zero columns.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
arma::vec default_precond(const arma::mat& B,
                          double nu,
                          double sigma,
                          double ceil_val = 1e8)
{
  // sum of squares of each column = diag(B^T B)
  arma::vec col_ss = arma::sum(arma::square(B), 0).t();
  double scale = (nu + 3.0) * sigma * sigma / (nu + 1.0);
  // Avoid division by zero on degenerate columns
  arma::vec safe_col_ss = col_ss;
  safe_col_ss.elem(arma::find(safe_col_ss < 1.0 / ceil_val)).fill(1.0 / ceil_val);
  arma::vec v = scale / safe_col_ss;
  v.elem(arma::find(v > ceil_val)).fill(ceil_val);
  return v;
}

// Helper: resolve the preconditioner.  If `precond` is empty (length 0),
// fall back to the default diagonal preconditioner; otherwise validate and
// return as-is.  This keeps the export signatures clean while making
// preconditioning the *default* behaviour.
static inline arma::vec resolve_precond(const arma::vec& precond,
                                        const arma::mat& B,
                                        double nu,
                                        double sigma,
                                        arma::uword p)
{
  if (precond.n_elem == 0) {
    return default_precond(B, nu, sigma);
  }
  if (precond.n_elem != p) {
    Rcpp::stop("Length of `precond` must equal the dimension of `start`.");
  }
  if (arma::any(precond <= 0.0)) {
    Rcpp::stop("`precond` must be strictly positive on every entry.");
  }
  return precond;
}


// ---------------------------------------------------------------------------
// Proximal HMC (partial proximal) with diagonal mass matrix
// ---------------------------------------------------------------------------
// `precond` is the posterior-variance diagonal v.  Mass matrix M = diag(1/v).
// Leapfrog with mass matrix:
//   p_{1/2} = p_0 + (eps/2) * grad_log_pi(q_0)
//   q_1     = q_0 + eps * M^{-1} p_{1/2}   = q_0 + eps * (p_{1/2} % v)
//   p_1     = p_{1/2} + (eps/2) * grad_log_pi(q_1)
// Momentum draw:   p ~ N(0, M)  =>  p = randn(p) / sqrt(v).
// Kinetic energy:  K(p) = 0.5 * p^T M^{-1} p = 0.5 * sum(p^2 * v).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List phmc_cpp(const arma::mat& B, const arma::vec& y,
              double lambda, double alpha, double sigma,
              int iter, double eps_hmc, int L, double nu,
              arma::vec start,
              Rcpp::Nullable<arma::vec> precond = R_NilValue,
              bool blather = true)
{
  int p = start.n_elem;
  
  arma::vec v = resolve_precond(
    precond.isNotNull() ? Rcpp::as<arma::vec>(precond) : arma::vec(),
    B, nu, sigma, p
  );
  arma::vec sqrt_inv_v = 1.0 / arma::sqrt(v);   // for momentum draw: p = randn / sqrt(v)
  // (We use `v` directly in the position update and kinetic-energy formulas.)
  
  arma::mat samples(iter, p);
  arma::vec samp = start;
  samples.row(0) = samp.t();
  
  // Standard-normal draws for momentum; scaled to N(0, M) below
  arma::mat mom_mat = arma::randn(iter, p);
  
  int accept = 0;
  
  for (int i = 1; i < iter; ++i) {
    arma::vec p_prop    = mom_mat.row(i).t() % sqrt_inv_v;   // N(0, diag(1/v)) = N(0, M)
    arma::vec p_current = p_prop + 0.5 * eps_hmc *
      grad_logpiLamg(B, y, samp, lambda, alpha, nu, sigma);
    arma::vec q_current = samp;
    
    int L_i = (R::runif(0, 1) < 0.05) ? 1 : L;
    for (int j = 0; j < L_i; ++j) {
      samp = samp + eps_hmc * (p_current % v);               // M^{-1} p = p % v
      arma::vec grad = grad_logpiLamg(B, y, samp, lambda, alpha, nu, sigma);
      if (j != L_i - 1)
        p_current += eps_hmc * grad;
    }
    
    p_current += 0.5 * eps_hmc *
      grad_logpiLamg(B, y, samp, lambda, alpha, nu, sigma);
    p_current = -p_current;
    
    double U_curr = -log_pi(q_current, y, B, nu, alpha, sigma);
    double U_prop = -log_pi(samp,      y, B, nu, alpha, sigma);
    // K(p) = 0.5 * p^T M^{-1} p = 0.5 * sum(p^2 * v)
    double K_curr = 0.5 * arma::dot(p_prop,    p_prop    % v);
    double K_prop = 0.5 * arma::dot(p_current, p_current % v);
    
    double log_acc = U_curr - U_prop + K_curr - K_prop;
    
    if (std::log(R::runif(0, 1)) <= log_acc) {
      samples.row(i) = samp.t();
      accept++;
    } else {
      samples.row(i) = q_current.t();
      samp = q_current;
    }
    
    if (blather) {
      if (i % std::max(1, iter / 10) == 0) {
        Rcpp::Rcout << "Iteration " << i
                    << ", Acceptance rate so far: "
                    << (double) accept / i << std::endl;
      }
    }
  }
  
  return List::create(Named("samples")     = samples,
                      Named("accept_rate") = (double) accept / iter,
                      Named("precond")     = v);
}

// ---------------------------------------------------------------------------
// Guo HMC with diagonal mass matrix (same convention as phmc_cpp)
// `precond` is the posterior-variance diagonal v; mass M = diag(1/v).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List guohmc_cpp(const arma::mat& B, const arma::vec& y,
                double lambda, double alpha, double sigma, int iter,
                double eps_hmc, int L, double nu,
                arma::vec start,
                Rcpp::Nullable<arma::vec> precond = R_NilValue,
                bool blather = true)
{
  int p = start.n_elem;
  
  arma::vec v = resolve_precond(
    precond.isNotNull() ? Rcpp::as<arma::vec>(precond) : arma::vec(),
    B, nu, sigma, p
  );
  arma::vec sqrt_inv_v = 1.0 / arma::sqrt(v);
  
  arma::mat samples(iter, p);
  arma::vec samp = start;
  samples.row(0) = samp.t();
  
  arma::mat mom_mat = arma::randn(iter, p);
  int accept = 0;
  
  for (int i = 1; i < iter; ++i) {
    arma::vec p_prop    = mom_mat.row(i).t() % sqrt_inv_v;   // N(0, M)
    arma::vec p_current = p_prop + 0.5 * eps_hmc *
      grad_logpiLam_guo(samp, lambda, alpha);
    arma::vec q_current = samp;
    
    int L_i = (R::runif(0, 1) < 0.05) ? 1 : L;
    for (int j = 0; j < L_i; ++j) {
      samp = samp + eps_hmc * (p_current % v);               // M^{-1} p
      arma::vec grad = grad_logpiLam_guo(samp, lambda, alpha);
      if (j != L_i - 1)
        p_current += eps_hmc * grad;
    }
    
    p_current += 0.5 * eps_hmc * grad_logpiLam_guo(samp, lambda, alpha);
    p_current = -p_current;
    
    double U_curr = -log_pi(q_current, y, B, nu, alpha, sigma);
    double U_prop = -log_pi(samp,      y, B, nu, alpha, sigma);
    double K_curr = 0.5 * arma::dot(p_prop,    p_prop    % v);
    double K_prop = 0.5 * arma::dot(p_current, p_current % v);
    
    double log_acc = U_curr - U_prop + K_curr - K_prop;
    
    if (std::log(R::runif(0, 1)) <= log_acc) {
      samples.row(i) = samp.t();
      accept++;
    } else {
      samples.row(i) = q_current.t();
      samp = q_current;
    }
    
    if (blather) {
      if (i % std::max(1, iter / 10) == 0) {
        Rcpp::Rcout << "Iteration " << i
                    << ", Acceptance rate so far: "
                    << (double) accept / i << std::endl;
      }
    }
  }
  
  return List::create(Named("samples")     = samples,
                      Named("accept_rate") = (double) accept / iter,
                      Named("precond")     = v);
}


// ---------------------------------------------------------------------------
// Preconditioned Random-Walk Metropolis
// `precond` is the posterior-variance diagonal v.
// Proposal:  w' = w + h * sqrt(v) % randn(p),  i.e. N(w, h^2 * diag(v)).
// Since the proposal is symmetric, the acceptance ratio only uses log_pi.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List rwm_cpp(const arma::mat& B, const arma::vec& y, int iter,
             double h, arma::vec start,
             double alpha, double sigma, double nu,
             Rcpp::Nullable<arma::vec> precond = R_NilValue,
             bool blather = true)
{
  int p = start.n_elem;
  
  arma::vec v = resolve_precond(
    precond.isNotNull() ? Rcpp::as<arma::vec>(precond) : arma::vec(),
    B, nu, sigma, p
  );
  arma::vec sqrt_v = arma::sqrt(v);
  
  arma::mat samples(iter, p);
  arma::vec samp = start;
  samples.row(0) = samp.t();
  int accept = 0;
  
  // Cache current log-target to avoid recomputing it every iteration
  double logpi_curr = log_pi(samp, y, B, nu, alpha, sigma);
  
  for (int i = 1; i < iter; ++i) {
    arma::vec prop = samp + h * (sqrt_v % arma::randn<arma::vec>(p));
    double logpi_prop = log_pi(prop, y, B, nu, alpha, sigma);
    double logr = logpi_prop - logpi_curr;
    
    if (std::log(R::runif(0, 1)) <= logr) {
      samples.row(i) = prop.t();
      samp = prop;
      logpi_curr = logpi_prop;
      accept++;
    } else {
      samples.row(i) = samp.t();
    }
    
    if (blather) {
      if (i % std::max(1, iter / 10) == 0) {
        Rcpp::Rcout << "Iteration " << i
                    << ", Acceptance rate so far: "
                    << (double) accept / i << std::endl;
      }
    }
  }
  
  return List::create(Named("samples")     = samples,
                      Named("accept_rate") = (double) accept / iter,
                      Named("precond")     = v);
}


// [[Rcpp::export]]
arma::vec map_estimate(const arma::mat& B, const arma::vec& y,
                       double alpha, double nu, double sigma,
                       arma::vec start,
                       int max_iter = 2000,
                       double tol = 1e-8) {
  int p = start.n_elem;
  arma::vec w = start;
  arma::vec y_k = w;
  arma::vec w_prev = w;
  double t_km1 = 1.0;
  
  // Negative log-likelihood (the smooth part we're minimizing)
  auto f_smooth = [&](const arma::vec& v) {
    arma::vec r = y - B * v;
    
    return 0.5 * (nu + 1.0) *
      arma::sum(
        arma::log(
          1.0 +
            arma::square(r) /
              (nu * sigma * sigma)
        )
      );
  };
  
  // Start with a small step, grow/shrink with backtracking
  double step = 1.0 / (((nu + 1.0) / (8.0 * nu * sigma * sigma)) *
                       std::pow(arma::norm(B, 2), 2.0));
  
  for (int k = 0; k < max_iter; ++k) {
    arma::vec g = gradf_dur(y_k, y, B, nu, sigma);   // gradient of NEG log-lik
    double f_yk = f_smooth(y_k);
    
    // Backtracking line search: ensure majorization
    arma::vec w_new;
    for (int bt = 0; bt < 50; ++bt) {
      arma::vec z = y_k - step * g;
      w_new = soft_threshold(z, alpha * step);
      arma::vec diff = w_new - y_k;
      double lhs = f_smooth(w_new);
      double rhs = f_yk + arma::dot(g, diff)
        + 0.5 / step * arma::dot(diff, diff);
      if (lhs <= rhs + 1e-12) break;
      step *= 0.5;
    }
    
    // FISTA momentum
    double t_k = 0.5 * (1.0 + std::sqrt(1.0 + 4.0 * t_km1 * t_km1));
    y_k = w_new + ((t_km1 - 1.0) / t_k) * (w_new - w);
    
    if (arma::norm(w_new - w, 2) < tol * std::max(1.0, arma::norm(w, 2))) {
      w = w_new; break;
    }
    w = w_new;
    t_km1 = t_k;
    
    // Optional: grow step occasionally to avoid getting too small
    if (k % 20 == 0) step *= 1.1;
  }
  return w;
}