// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// // [[Rcpp::export]]
// arma::vec calc_prox_fista(const arma::vec& x,         // Prox center
//                           const arma::mat& X,         // n x p
//                           const arma::vec& y,         // binary labels
//                           double lambda = 1.0,        // prox weight
//                           double alpha = 1.0,         // L1 penalty weight
//                           double step_size = 1e-3,
//                           int max_iter = 500,
//                           double tol = 1e-6) {
//   
//   int n = X.n_rows;
//   int p = X.n_cols;
//   
//   arma::vec z = x;          // initialize z_0 = x
//   arma::vec z_prev = z;     // z_{k-1}
//   arma::vec y_k = z;        // extrapolated point
//   arma::vec grad(p);
//   arma::vec z_new = z;
//   
//   double t_k = 1.0;
//   double t_km1 = 1.0;
//   
//   int iter;
//   
//   for (iter = 0; iter < max_iter; ++iter) {
//     // Logistic gradient at y_k
//     arma::vec Xy = X * y_k;
//     arma::vec sigmoid = 1.0 / (1.0 + exp(-Xy));
//     grad = X.t() * (sigmoid - y) / n;
//     
//     // Gradient step on smooth part: includes z - x term
//     arma::vec temp = y_k - step_size * (y_k - x + lambda * grad);
//     
//     // Proximal step: soft-thresholding
//     for (int j = 0; j < p; ++j) {
//       double thresh = lambda * alpha * step_size;
//       double val = temp[j];
//       if (val > thresh)
//         z_new[j] = val - thresh;
//       else if (val < -thresh)
//         z_new[j] = val + thresh;
//       else
//         z_new[j] = 0.0;
//     }
//     
//     // Convergence check
//     if (norm(z_new - z, 2) < tol)
//       break;
//     
//     // Update momentum
//     t_k = (1 + std::sqrt(1 + 4 * t_km1 * t_km1)) / 2.0;
//     double momentum = (t_km1 - 1) / t_k;
//     
//     y_k = z_new + momentum * (z_new - z);
//     
//     // Update iterates
//     z_prev = z;
//     z = z_new;
//     t_km1 = t_k;
//   }
//   
//   return z_new;
// }

// [[Rcpp::export]]
// Gradient of f for Durmus
arma::vec gradf_dur(const arma::vec& w, const arma::vec& y,
                    const arma::mat& B, double nu)
{
  // residual vector
  arma::vec r = y - B * w;
  
  // element-wise: r_i / (nu + r_i^2)
  arma::vec weights = r / (nu + r % r);
  
  // gradient: - B^T * weights
  arma::vec grad = - (nu + 1) *B.t() * weights;
  
  return grad;
}
// arma::vec gradf_dur(const arma::vec& z, const arma::mat& x, const arma::vec& y) {
//   arma::vec xz = x * z;
//   arma::vec sigmoid = 1 / (1 + exp(-xz));
//   return x.t() * (sigmoid - y);
// }


// log target
// [[Rcpp::export]]
double log_pi(const arma::vec& w, const arma::vec& y, const arma::mat& B,
               double nu, double alpha)
{
  // residuals
  arma::vec r = y - B * w;
  
  // first term:
  // sum log(1 + r_i^2 / (nu * sigma^2))
  double log_term = 0.5 * (nu + 1) * sum(log(1.0 + square(r) / nu ));
  
  // L1 penalty
  double penalty = alpha * norm(w, 1);
  
  // negative objective
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
arma::vec grad_logpiLamg(const arma::mat& B, const arma::vec& y,
                         const arma::vec& w, double lambda, double alpha, double nu) {
  //arma::vec eta = X * beta;
  arma::vec grad_f = gradf_dur(w, y, B, nu); 
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
List phmc_cpp(const arma::mat& B, const arma::vec& y, double lambda, double alpha,
              int iter, double eps_hmc, int L, double nu, arma::vec start, bool blather = true) {
  int p = start.n_elem;
  arma::mat samples(iter, p);
  arma::vec samp = start;
  samples.row(0) = samp.t();
  arma::mat mom_mat = arma::randn(iter, p);
  int accept = 0;
  
  for (int i = 1; i < iter; ++i) {
    arma::vec p_prop = mom_mat.row(i).t();
    arma::vec p_current = p_prop + 0.5 * eps_hmc * grad_logpiLamg(B, y, samp, lambda, alpha, nu);
    arma::vec q_current = samp;
    
    for (int j = 0; j < L; ++j) {
      samp = samp + eps_hmc * p_current;
      arma::vec grad = grad_logpiLamg(B, y, samp, lambda, alpha, nu);
      if (j != L - 1)
        p_current += eps_hmc * grad;
    }
    
    p_current += 0.5 * eps_hmc * grad_logpiLamg(B, y, samp, lambda, alpha, nu);
    p_current = -p_current;
    
    double U_curr = -log_pi(q_current, y, B, nu, alpha);
    double U_prop = -log_pi(samp, y, B, nu, alpha);
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
List guohmc_cpp(const arma::mat& B, const arma::vec& y, double lambda, double alpha,int iter,
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
    
    double U_curr = -log_pi(q_current, y, B, nu, alpha);
    double U_prop = -log_pi(samp, y, B, nu, alpha);
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
                                double alpha, double nu, bool blather = true) {
  int p = start.n_elem;
  arma::mat samples(iter, p);
  arma::vec samp = start;
  samples.row(0) = samp.t();
  int accept = 0;
  
  for (int i = 1; i < iter; ++i) {
    arma::vec prop = samp + h * arma::randn<arma::vec>(p);
    double logr = log_pi(prop, y, B, nu, alpha) - log_pi(samp, y, B, nu, alpha);
    
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

