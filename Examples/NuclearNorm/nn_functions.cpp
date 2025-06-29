// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// Compute nuclear norm of a square matrix from a vector
double nucl_norm(const arma::vec& vect) {
  int nvar = vect.n_elem;
  int n = sqrt(nvar);
  arma::mat A = reshape(vect, n, n);
  arma::vec s;
  svd(s, A);
  return sum(s);
}


double log_pi(const arma::vec& x, const arma::vec& y, double sigma2, double alpha) {
  double n_norm = nucl_norm(x);
  double dens_val = alpha * n_norm + sum(square(y - x)) / (2.0 * sigma2);
  return -dens_val;
}


double log_pilambda(const arma::vec& eta, const arma::vec& x, double lambda,
                    const arma::vec& y, double sigma2, double alpha) {
  double n_norm = nucl_norm(eta);
  double dens_val = alpha * n_norm + sum(square(eta - x)) / (2.0 * lambda)
    + sum(square(y - eta)) / (2.0 * sigma2);
  return -dens_val;
}

// Soft-thresholding
arma::vec softthreshold(const arma::vec& u, double lambda) {
  arma::vec result = abs(u) - lambda;
  result.elem(find(result < 0)).zeros();
  return result % sign(u);
}


arma::vec prox_func(const arma::vec& x, double lambda,
                    const arma::vec& y, double sigma2, double alpha) {
  arma::vec num = lambda * y + sigma2 * x;
  double denom = lambda + sigma2;
  int nvar = x.n_elem;
  int n = sqrt(nvar);
  arma::mat mat = reshape(num / denom, n, n);
  arma::mat U, V;
  arma::vec s;
  svd(U, s, V, mat);
  arma::vec s_thresh = softthreshold(s, (alpha * sigma2 * lambda) / denom);
  arma::mat out = U * diagmat(s_thresh) * V.t();
  return vectorise(out);
}

arma::vec dur_prox_func(const arma::vec& x, double lambda, double alpha) {
  int nvar = x.n_elem;
  int n = sqrt(nvar);
  arma::mat mat = reshape(x, n, n);
  arma::mat U, V;
  arma::vec s;
  svd(U, s, V, mat);
  arma::vec s_thresh = s - lambda * alpha;
  s_thresh.elem(find(s_thresh < 0)).zeros();
  arma::mat out = U * diagmat(s_thresh) * V.t();
  return vectorise(out);
}


arma::vec grad_logpiLam(const arma::vec& x, double lambda, const arma::vec& y,
                        double sigma2, double alpha) {
  arma::vec x_prox = prox_func(x, lambda, y, sigma2, alpha);
  return -(x - x_prox) / lambda;
}


arma::vec grad_log_durpiLam(const arma::vec& x, double lambda, const arma::vec& y,
                            double sigma2, double alpha) {
  arma::vec x_prox = dur_prox_func(x, lambda, alpha);
  arma::vec term2 = -(x - x_prox) / lambda;
  arma::vec term1 = (y - x) / sigma2;
  return term1 + term2;
}

// [[Rcpp::export]]
List durhmc_cpp(const arma::vec& y, double alpha, double lambda, double sigma2,
            int iter, double eps_hmc, int L, arma::vec start) {
  
  int nvar = y.n_elem;
  arma::mat samp_hmc(iter, nvar, fill::zeros);
  arma::mat mom_mat = randn<arma::mat>(iter, nvar);
  
  arma::vec samp = start;
  samp_hmc.row(0) = samp.t();
  
  int accept = 0;
  
  for (int i = 1; i < iter; ++i) {
    arma::vec p_prop = mom_mat.row(i).t();
    arma::vec U_samp = -grad_log_durpiLam(samp, lambda, y, sigma2, alpha);
    arma::vec p_current = p_prop - 0.5 * eps_hmc * U_samp;
    arma::vec q_current = samp;
    
    for (int j = 0; j < L; ++j) {
      samp += eps_hmc * p_current;
      U_samp = -grad_log_durpiLam(samp, lambda, y, sigma2, alpha);
      if (j != L - 1) {
        p_current -= eps_hmc * U_samp;
      }
    }
    
    p_current -= 0.5 * eps_hmc * U_samp;
    p_current = -p_current;
    
    double U_curr = sum(square(y - q_current)) / (2.0 * sigma2) + alpha * nucl_norm(q_current);
    double U_prop = sum(square(y - samp)) / (2.0 * sigma2) + alpha * nucl_norm(samp);
    double K_curr = dot(p_prop, p_prop) / 2.0;
    double K_prop = dot(p_current, p_current) / 2.0;
    
    double log_acc_prob = U_curr - U_prop + K_curr - K_prop;
    
    if (std::log(R::runif(0, 1)) <= log_acc_prob) {
      samp_hmc.row(i) = samp.t();
      accept++;
    } else {
      samp_hmc.row(i) = q_current.t();
      samp = q_current;
    }
    
    if (i % (iter / 10) == 0) {
      Rcout << "Iter: " << i << "  Accept Rate: " << static_cast<double>(accept) / i << std::endl;
    }
  }
  
  double acc_rate = static_cast<double>(accept) / iter;
  Rcout << "Final Acceptance Rate: " << acc_rate << std::endl;
  return List::create(Named("samples") = samp_hmc, Named("acceptance_rate") = acc_rate);
}

// [[Rcpp::export]]
List pxhmc_cpp(const arma::vec& y, double alpha, double lambda, double sigma2,
           int iter, double eps_hmc, int L, arma::vec start) {
  
  int nvar = y.n_elem;
  arma::mat samp_hmc(iter, nvar, fill::zeros);
  arma::mat mom_mat = randn<arma::mat>(iter, nvar);
  
  arma::vec samp = start;
  samp_hmc.row(0) = samp.t();
  
  int accept = 0;
  
  for (int i = 1; i < iter; ++i) {
    arma::vec p_prop = mom_mat.row(i).t();
    arma::vec U_samp = -grad_logpiLam(samp, lambda, y, sigma2, alpha);
    arma::vec p_current = p_prop - 0.5 * eps_hmc * U_samp;
    arma::vec q_current = samp;
    
    for (int j = 0; j < L; ++j) {
      samp += eps_hmc * p_current;
      U_samp = -grad_logpiLam(samp, lambda, y, sigma2, alpha);
      if (j != L - 1) {
        p_current -= eps_hmc * U_samp;
      }
    }
    
    p_current -= 0.5 * eps_hmc * U_samp;
    p_current = -p_current;
    
    double U_curr = sum(square(y - q_current)) / (2.0 * sigma2) + alpha * nucl_norm(q_current);
    double U_prop = sum(square(y - samp)) / (2.0 * sigma2) + alpha * nucl_norm(samp);
    double K_curr = dot(p_prop, p_prop) / 2.0;
    double K_prop = dot(p_current, p_current) / 2.0;
    
    double log_acc_prob = U_curr - U_prop + K_curr - K_prop;
    
    if (std::log(R::runif(0, 1)) <= log_acc_prob) {
      samp_hmc.row(i) = samp.t();
      accept++;
    } else {
      samp_hmc.row(i) = q_current.t();
      samp = q_current;
    }
    
    if (i % (iter / 10) == 0) {
      Rcout << "Iter: " << i << "  Accept Rate: " << static_cast<double>(accept) / i << std::endl;
    }
  }
  
  double acc_rate = static_cast<double>(accept) / iter;
  Rcout << "Final Acceptance Rate: " << acc_rate << std::endl;
  return List::create(Named("samples") = samp_hmc, Named("acceptance_rate") = acc_rate);
}
