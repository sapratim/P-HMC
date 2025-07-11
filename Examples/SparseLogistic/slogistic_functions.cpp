// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// not used
Rcpp::List calc_prox(const arma::vec& x,         // Prox center
                     const arma::mat& X,         // n x p design matrix
                     const arma::vec& y,         // n binary labels
                     double lambda = 1.0,        // overall proximal weight
                     double alpha = 1.0,         // L1 weight inside penalty
                     double step_size = 1e-3, 
                     int max_iter = 500,
                     double tol = 1e-6) {
  
  int n = X.n_rows;
  int p = X.n_cols;
  arma::vec z = x;              // initialize at prox center
  arma::vec z_new = z;
  arma::vec grad(p);
  int iter;
  
  for (iter = 0; iter < max_iter; ++iter) {
    // Gradient of logistic loss
    arma::vec Xz = X * z;
    arma::vec sigmoid = 1.0 / (1.0 + exp(-Xz));
    grad = X.t() * (sigmoid - y) / n;
    
    // Gradient step on smooth part: prox center + logistic loss
    arma::vec temp = z - step_size * (z - x + lambda * grad);
    
    // Prox step on L1: soft-thresholding
    for (int j = 0; j < p; ++j) {
      double thresh = lambda * alpha * step_size;
      double val = temp[j];
      if (val > thresh)
        z_new[j] = val - thresh;
      else if (val < -thresh)
        z_new[j] = val + thresh;
      else
        z_new[j] = 0.0;
    }
    
    // Convergence check
    if (norm(z_new - z, 2) < tol)
      break;
    
    z = z_new;
  }
  
  return List::create(
    Named("z") = z_new,
    Named("n_iter") = iter + 1  // actual number of iterations performed
  );
}


#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec calc_prox_fista(const arma::vec& x,         // Prox center
                           const arma::mat& X,         // n x p
                           const arma::vec& y,         // binary labels
                           double lambda = 1.0,        // prox weight
                           double alpha = 1.0,         // L1 penalty weight
                           double step_size = 1e-3,
                           int max_iter = 500,
                           double tol = 1e-6) {
  
  int n = X.n_rows;
  int p = X.n_cols;
  
  arma::vec z = x;          // initialize z_0 = x
  arma::vec z_prev = z;     // z_{k-1}
  arma::vec y_k = z;        // extrapolated point
  arma::vec grad(p);
  arma::vec z_new = z;
  
  double t_k = 1.0;
  double t_km1 = 1.0;
  
  int iter;
  
  for (iter = 0; iter < max_iter; ++iter) {
    // Logistic gradient at y_k
    arma::vec Xy = X * y_k;
    arma::vec sigmoid = 1.0 / (1.0 + exp(-Xy));
    grad = X.t() * (sigmoid - y) / n;
    
    // Gradient step on smooth part: includes z - x term
    arma::vec temp = y_k - step_size * (y_k - x + lambda * grad);
    
    // Proximal step: soft-thresholding
    for (int j = 0; j < p; ++j) {
      double thresh = lambda * alpha * step_size;
      double val = temp[j];
      if (val > thresh)
        z_new[j] = val - thresh;
      else if (val < -thresh)
        z_new[j] = val + thresh;
      else
        z_new[j] = 0.0;
    }
    
    // Convergence check
    if (norm(z_new - z, 2) < tol)
      break;
    
    // Update momentum
    t_k = (1 + std::sqrt(1 + 4 * t_km1 * t_km1)) / 2.0;
    double momentum = (t_km1 - 1) / t_k;
    
    y_k = z_new + momentum * (z_new - z);
    
    // Update iterates
    z_prev = z;
    z = z_new;
    t_km1 = t_k;
  }
  
  return z_new;
}

// [[Rcpp::export]]
// Gradient of f for Durmus
arma::vec gradf_dur(const arma::vec& z, const arma::mat& x, const arma::vec& y) {
  arma::vec xz = x * z;
  arma::vec sigmoid = 1 / (1 + exp(-xz));
  return x.t() * (sigmoid - y);
}


// log target

double log_pi(const arma::mat& X, const arma::vec& y, const arma::vec& beta, double alpha) {
  arma::vec eta = X * beta;
  arma::vec loglik = y % log1p(exp(-eta)) + (1 - y) % (eta + log1p(exp(-eta)));
  double penalty = alpha * sum(abs(beta));
  return -sum(loglik) - penalty;
}

// soft thresholding function

arma::vec soft_threshold(const arma::vec& z, double threshold) {
  return arma::sign(z) % arma::max(arma::abs(z) - threshold, arma::zeros(z.n_elem));
}
// arma::vec soft_threshold(const arma::vec& u, double pen) {
//   arma::vec out = u;
//   for (unsigned int i = 0; i < u.n_elem; ++i) {
//     if (u[i] > pen)
//       out[i] = u[i] - pen;
//     else if (u[i] < -pen)
//       out[i] = u[i] + pen;
//     else
//       out[i] = 0.0;
//   }
//   return out;
// }


// calculate gradient of the ns-hmc method

arma::vec grad_logpiLam(const arma::mat& X, const arma::vec& y,
                        const arma::vec& beta, double lambda,
                        double alpha, double step_size = 0.1,
                        int max_iter = 100, double tol = 1e-6) {

  arma::vec beta_prox = calc_prox_fista(beta, X, y, lambda, alpha, step_size = .1);
  // for (int iter = 0; iter < max_iter; ++iter) {
  //   arma::vec eta = X * z;
  //   arma::vec grad = X.t() * (1.0 / (1.0 + exp(-eta)) - y) / X.n_rows;
  //   arma::vec temp = z - step_size * (grad);
  //   z_new = soft_threshold(temp, thresh);
  //   if (norm(z_new - z, 2) < tol)
  //     break;
  //   z = z_new;
  // }
  return -(beta - beta_prox) / lambda;
}

// proximal mapping for phmc

arma::vec prox_phmc(const arma::vec& beta, double lambda, double alpha) {
  return soft_threshold(beta, lambda * alpha);
}

// gradient of the partial proximal methods

arma::vec grad_logpiLamg(const arma::mat& X, const arma::vec& y,
                         const arma::vec& beta, double lambda, double alpha) {
  //arma::vec eta = X * beta;
  arma::vec grad_f = gradf_dur(beta, X, y); // X.t() * (1.0 / (1.0 + exp(-eta)) - y) / X.n_rows;
  arma::vec beta_prox = prox_phmc(beta, lambda, alpha);
  return -(grad_f + (beta - beta_prox) / lambda);
}




// Chaari's non-smooth method
// [[Rcpp::export]]
List nshmc_cpp(const arma::mat& X, const arma::vec& y, double lambda, double alpha,
               int iter, double eps_hmc, int L, arma::vec start) {
  int p = start.n_elem;
  arma::mat samples(iter, p);
  arma::vec samp = start;
  samples.row(0) = samp.t();
  arma::mat mom_mat = arma::randn(iter, p);
  int accept = 0;
  
  for (int i = 1; i < iter; ++i) {
    arma::vec p_prop = mom_mat.row(i).t();
    arma::vec p_current = p_prop - 0.5 * eps_hmc * grad_logpiLam(X, y, samp, lambda, alpha);
    arma::vec q_current = samp;
    
    for (int j = 0; j < L; ++j) {
      samp = samp + eps_hmc * p_current;
      arma::vec grad = grad_logpiLam(X, y, samp, lambda, alpha);
      if (j != L - 1)
        p_current -= eps_hmc * grad;
    }
    
    p_current -= 0.5 * eps_hmc * grad_logpiLam(X, y, samp, lambda, alpha);
    p_current = -p_current;
    
    double U_curr = -log_pi(X, y, q_current, alpha);
    double U_prop = -log_pi(X, y, samp, alpha);
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
    
    if (i % std::max(1, iter / 10) == 0) {
      Rcpp::Rcout << "Iteration " << i
                  << ", Acceptance rate so far: "
                  << (double) accept / i << std::endl;
    }
  }
  
  return List::create(Named("samples") = samples,
                      Named("accept_rate") = (double) accept / iter);
}





// Proximal HMC (paritial proximal)
// [[Rcpp::export]]
List phmc_cpp(const arma::mat& X, const arma::vec& y, double lambda, double alpha,
              int iter, double eps_hmc, int L, arma::vec start) {
  int p = start.n_elem;
  arma::mat samples(iter, p);
  arma::vec samp = start;
  samples.row(0) = samp.t();
  arma::mat mom_mat = arma::randn(iter, p);
  int accept = 0;
  
  for (int i = 1; i < iter; ++i) {
    arma::vec p_prop = mom_mat.row(i).t();
    arma::vec p_current = p_prop + 0.5 * eps_hmc * grad_logpiLamg(X, y, samp, lambda, alpha);
    arma::vec q_current = samp;
    
    for (int j = 0; j < L; ++j) {
      samp = samp + eps_hmc * p_current;
      arma::vec grad = grad_logpiLamg(X, y, samp, lambda, alpha);
      if (j != L - 1)
        p_current += eps_hmc * grad;
    }
    
    p_current += 0.5 * eps_hmc * grad_logpiLamg(X, y, samp, lambda, alpha);
    p_current = -p_current;
    
    double U_curr = -log_pi(X, y, q_current, alpha);
    double U_prop = -log_pi(X, y, samp, alpha);
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
    
    if (i % std::max(1, iter / 10) == 0) {
      Rcpp::Rcout << "Iteration " << i
                  << ", Acceptance rate so far: "
                  << (double) accept / i << std::endl;
    }
  }
  
  return List::create(Named("samples") = samples,
                      Named("accept_rate") = (double) accept / iter);
}



// [[Rcpp::export]]
List rwm_cpp(const arma::mat& X, const arma::vec& y, int iter,
             double h, arma::vec start, double alpha) {
  int p = start.n_elem;
  arma::mat samples(iter, p);
  arma::vec samp = start;
  samples.row(0) = samp.t();
  int accept = 0;
  
  for (int i = 1; i < iter; ++i) {
    arma::vec prop = samp + h * arma::randn<arma::vec>(p);
    double logr = log_pi(X, y, prop, alpha) - log_pi(X, y, samp, alpha);
    
    if (std::log(R::runif(0, 1)) <= logr) {
      samples.row(i) = prop.t();
      samp = prop;
      accept++;
    } else {
      samples.row(i) = samp.t();
    }
    
    if (i % std::max(1, iter / 10) == 0) {
      Rcpp::Rcout << "Iteration " << i
                  << ", Acceptance rate so far: "
                  << (double) accept / i << std::endl;
    }
  }
  
  return List::create(Named("samples") = samples,
                      Named("accept_rate") = (double) accept / iter);
}

