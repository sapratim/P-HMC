// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

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
  
  double log_term = 0.5 * (nu + 1.0) *sum(log(1.0 +square(r) / (nu * sigma * sigma)));
  
  // Optional normalization term
  log_term += y.n_elem * std::log(sigma);
  
  double penalty = alpha * norm(w, 1);
  
  return -(log_term + penalty);
}

// soft thresholding function
arma::vec soft_threshold(const arma::vec& z, double threshold) {
  return arma::sign(z) % arma::max(arma::abs(z) - threshold, arma::zeros(z.n_elem));
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
  //arma::vec eta = X * beta;
  arma::vec w_prox = prox_phmc(w, lambda, alpha);
  return - (w - w_prox);
}

// Proximal HMC (partial proximal)
// [[Rcpp::export]]
List phmc_cpp(const arma::mat& B, const arma::vec& y, double lambda, double alpha, double sigma,
              int iter, double eps_hmc, int L, double nu, arma::vec start, bool blather = true) {
  int p = start.n_elem;
  arma::mat samples(iter, p);
  arma::vec samp = start;
  samples.row(0) = samp.t();
  arma::mat mom_mat = arma::randn(iter, p);
  int accept = 0;
  
  for (int i = 1; i < iter; ++i) {
    arma::vec p_prop = mom_mat.row(i).t();
    arma::vec p_current = p_prop + 0.5 * eps_hmc * grad_logpiLamg(B, y, samp, lambda, alpha, nu, sigma);
    arma::vec q_current = samp;
    
    for (int j = 0; j < L; ++j) {
      samp = samp + eps_hmc * p_current;
      arma::vec grad = grad_logpiLamg(B, y, samp, lambda, alpha, nu, sigma);
      if (j != L - 1)
        p_current += eps_hmc * grad;
    }
    
    p_current += 0.5 * eps_hmc * grad_logpiLamg(B, y, samp, lambda, alpha, nu, sigma);
    p_current = -p_current;
    
    double U_curr = -log_pi(q_current, y, B, nu, alpha, sigma);
    double U_prop = -log_pi(samp, y, B, nu, alpha, sigma);
    double K_curr = 0.5 * dot(p_prop, p_prop);
    double K_prop = 0.5 * dot(p_current, p_current);
    
    double log_acc = U_curr - U_prop + K_curr - K_prop;
    
    if (std::log(R::runif(0, 1)) <= log_acc) {
      samples.row(i) = samp.t();
      accept++;
    } else {
      samples.row(i) = q_current.t();
      samp = q_current;
    }
    
    if(blather)
    {
      if (i % std::max(1, iter / 10) == 0) {
        Rcpp::Rcout << "Iteration " << i
                    << ", Acceptance rate so far: "
                    << (double) accept / i << std::endl;
      }
    }
  }
  
  return List::create(Named("samples") = samples,
                      Named("accept_rate") = (double) accept / iter);
}

// Guo HMC 
// [[Rcpp::export]]
List guohmc_cpp(const arma::mat& B, const arma::vec& y, double lambda, double alpha, double sigma, int iter,
                double eps_hmc, int L, double nu, arma::vec start, bool blather = true) {
  int p = start.n_elem;
  arma::mat samples(iter, p);
  arma::vec samp = start;
  samples.row(0) = samp.t();
  arma::mat mom_mat = arma::randn(iter, p);
  int accept = 0;
  
  for (int i = 1; i < iter; ++i) {
    arma::vec p_prop = mom_mat.row(i).t();
    arma::vec p_current = p_prop + 0.5 * eps_hmc * grad_logpiLam_guo(samp, lambda, alpha);
    arma::vec q_current = samp;
    
    for (int j = 0; j < L; ++j) {
      samp = samp + eps_hmc * p_current;
      arma::vec grad = grad_logpiLam_guo(samp, lambda, alpha);
      if (j != L - 1)
        p_current += eps_hmc * grad;
    }
    
    p_current += 0.5 * eps_hmc * grad_logpiLam_guo(samp, lambda, alpha);
    p_current = -p_current;
    
    double U_curr = -log_pi(q_current, y, B, nu, alpha, sigma);
    double U_prop = -log_pi(samp, y, B, nu, alpha, sigma);
    double K_curr = 0.5 * dot(p_prop, p_prop);
    double K_prop = 0.5 * dot(p_current, p_current);
    
    double log_acc = U_curr - U_prop + K_curr - K_prop;
    
    if (std::log(R::runif(0, 1)) <= log_acc) {
      samples.row(i) = samp.t();
      accept++;
    } else {
      samples.row(i) = q_current.t();
      samp = q_current;
    }
    
    if(blather)
    {
      if (i % std::max(1, iter / 10) == 0) {
        Rcpp::Rcout << "Iteration " << i
                    << ", Acceptance rate so far: "
                    << (double) accept / i << std::endl;
      }
    }
  }
  
  return List::create(Named("samples") = samples,
                      Named("accept_rate") = (double) accept / iter);
}


// [[Rcpp::export]]
List rwm_cpp(const arma::mat& B, const arma::vec& y, int iter,double h, arma::vec start,
                                double alpha, double sigma, double nu, bool blather = true) {
  int p = start.n_elem;
  arma::mat samples(iter, p);
  arma::vec samp = start;
  samples.row(0) = samp.t();
  int accept = 0;
  
  for (int i = 1; i < iter; ++i) {
    arma::vec prop = samp + h * arma::randn<arma::vec>(p);
    double logr = log_pi(prop, y, B, nu, alpha, sigma) - log_pi(samp, y, B, nu, alpha, sigma);
    
    if (std::log(R::runif(0, 1)) <= logr) {
      samples.row(i) = prop.t();
      samp = prop;
      accept++;
    } else {
      samples.row(i) = samp.t();
    }
    
    if(blather)
    {
      if (i % std::max(1, iter / 10) == 0) {
        Rcpp::Rcout << "Iteration " << i
                    << ", Acceptance rate so far: "
                    << (double) accept / i << std::endl;
      }
    }
  }
  
  return List::create(Named("samples") = samples,
                      Named("accept_rate") = (double) accept / iter);
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
  double step = 1.0 /(((nu + 1.0) /(8.0 * nu * sigma * sigma))*std::pow(arma::norm(B, 2), 2.0));
  
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