// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Soft-thresholding function
arma::vec softthreshold(const arma::vec& z, double threshold) {
  return arma::sign(z) % arma::max(arma::abs(z) - threshold, arma::zeros(z.n_elem));
}

// Objective function f
double f(const arma::vec& z, const arma::mat& x, const arma::vec& y, 
         const arma::vec& beta_point, double lamb) {
  arma::vec xz = x * z;
  return accu(log(1 + exp(xz)) - y % xz) + accu(square(beta_point - z)) / (2 * lamb);
}

// Gradient of f
arma::vec gradf(const arma::vec& z, const arma::mat& x, const arma::vec& y, 
                const arma::vec& beta_point, double lamb) {
  arma::vec xz = x * z;
  arma::vec sigmoid = 1 / (1 + exp(-xz));
  return x.t() * (sigmoid - y) + (beta_point - z) / lamb;
}

// [[Rcpp::export]]
// Gradient of f for Durmus
arma::vec gradf_dur(const arma::vec& z, const arma::mat& x, const arma::vec& y) {
  arma::vec xz = x * z;
  arma::vec sigmoid = 1 / (1 + exp(-xz));
  return x.t() * (sigmoid - y);
}

// Function g
double g(const arma::vec& z, double alpha) {
  return alpha * accu(abs(z));
}

// Proximal operator for g
arma::vec proxg(const arma::vec& z, double tau_fasta, double alpha) {
  return softthreshold(z, alpha * tau_fasta);
}

// [[Rcpp::export]]
List rcpp_fasta(const arma::mat& x, const arma::vec& y, 
                arma::vec x0, const arma::vec& beta_point, 
                double alpha, double lamb, double tau1, 
                int max_iters = 200, int w = 10, bool backtrack = true, 
                bool recordIterates = false, double stepsizeShrink = 0.1, 
                double eps_n = 1e-15) {
  
  int n = x0.n_elem;
  arma::vec residual(max_iters, fill::zeros);
  arma::vec normalizedResid(max_iters, fill::zeros);
  arma::vec taus(max_iters, fill::zeros);
  arma::vec fVals(max_iters, fill::zeros);
  arma::vec objective(max_iters + 1, fill::zeros);
  arma::mat iterates;
  
  if (recordIterates) {
    iterates = arma::mat(n, max_iters + 1, fill::zeros);
    iterates.col(0) = x0;
  }
  
  int totalBacktracks = 0;
  double maxResidual = -std::numeric_limits<double>::infinity();
  double minObjectiveValue = std::numeric_limits<double>::infinity();
  
  arma::vec x1 = x0;
  arma::vec d1 = x1;
  double f1 = f(d1, x, y, beta_point, lamb);
  fVals(0) = f1;
  arma::vec gradf1 = gradf(d1, x, y, beta_point, lamb);
  objective(0) = f1 + g(x0, alpha);
  arma::vec bestObjectiveIterate = x1;
  
  for (int i = 0; i < max_iters; i++) {
    arma::vec x0 = x1;
    arma::vec gradf0 = gradf1;
    double tau0 = tau1;
    
    arma::vec x1hat = x0 - tau0 * gradf0;
    x1 = proxg(x1hat, tau0, alpha);
    arma::vec Dx = x1 - x0;
    d1 = x1;
    f1 = f(d1, x, y, beta_point, lamb);
    
    if (backtrack) {
      double M = fVals.subvec(std::max(0, i - w), std::max(0, i - 1)).max();
      int backtrackCount = 0;
      bool prop = (f1 - 1e-12 > M + dot(Dx, gradf0) + 0.5 * norm(Dx, 2) * norm(Dx, 2) / tau0) &&
        (backtrackCount < 20);
      
      while (prop) {
        tau0 *= stepsizeShrink;
        x1hat = x0 - tau0 * gradf0;
        x1 = proxg(x1hat, tau0, alpha);
        d1 = x1;
        f1 = f(d1, x, y, beta_point, lamb);
        Dx = x1 - x0;
        backtrackCount++;
        prop = (f1 - 1e-12 > M + dot(Dx, gradf0) + 0.5 * norm(Dx, 2) * norm(Dx, 2) / tau0) &&
          (backtrackCount < 20);
      }
      totalBacktracks += backtrackCount;
    }
    
    taus(i) = tau0;
    residual(i) = norm(Dx, 2) / tau0;
    maxResidual = std::max(maxResidual, residual(i));
    double normalizer = std::max(norm(gradf0, 2), norm(x1 - x1hat, 2) / tau0) + eps_n;
    normalizedResid(i) = residual(i) / normalizer;
    fVals(i) = f1;
    objective(i + 1) = f1 + g(x1, alpha);
    
    if (objective(i + 1) < minObjectiveValue) {
      bestObjectiveIterate = x1;
      minObjectiveValue = objective(i + 1);
    }
    
    gradf1 = gradf(d1, x, y, beta_point, lamb);
    arma::vec Dg = gradf1 + (x1hat - x0) / tau0;
    double dotprod = dot(Dx, Dg);
    double tau_s = norm(Dx, 2) * norm(Dx, 2) / dotprod;
    double tau_m = dotprod / (norm(Dg, 2) * norm(Dg, 2));
    tau_m = std::max(tau_m, 0.0);
    
    if (std::abs(dotprod) < 1e-15) break;
    
    if (2 * tau_m > tau_s) {
      tau1 = tau_m;
    } else {
      tau1 = tau_s - 0.5 * tau_m;
    }
    
    if (tau1 <= 0 || std::isinf(tau1) || std::isnan(tau1)) {
      tau1 = tau0 * 1.5;
    }
    
    if (recordIterates) {
      iterates.col(i + 1) = x1;
    }
  }
  
  if (recordIterates) {
    iterates = iterates.cols(0, max_iters);
  }
  
  return List::create(
    Named("x") = bestObjectiveIterate,
    Named("objective") = objective.subvec(0, max_iters),
    Named("fVals") = fVals,
    Named("totalBacktracks") = totalBacktracks,
    Named("residual") = residual,
    Named("taus") = taus,
    Named("iterates") = iterates
  );
}
