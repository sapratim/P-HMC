// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
List optimized_algorithm(
    Function f, Function gradf, Function g, Function proxg,
    arma::vec x0, double tau1, int max_iters = 100, int w = 10,
    bool backtrack = true, bool recordIterates = false, double stepsizeShrink = 0.5,
    double eps_n = 1e-15
) {
  arma::vec residual(max_iters, arma::fill::zeros);
  arma::vec normalizedResid(max_iters, arma::fill::zeros);
  arma::vec taus(max_iters, arma::fill::zeros);
  arma::vec fVals(max_iters, arma::fill::zeros);
  arma::vec objective(max_iters + 1, arma::fill::zeros);
  
  int totalBacktracks = 0, backtrackCount = 0;
  arma::vec x1 = x0;
  arma::vec d1 = x1;
  
  double f1 = as<double>(f(d1));
  fVals(0) = f1;
  arma::vec gradf1 = as<arma::vec>(gradf(d1));
  
  arma::mat iterates;
  if (recordIterates) {
    iterates = arma::mat(x0.n_elem, max_iters + 1, arma::fill::zeros);
    iterates.col(0) = x1;
  }
  
  double maxResidual = -arma::datum::inf;
  double minObjectiveValue = arma::datum::inf;
  objective(0) = f1 + as<double>(g(x0));
  
  arma::vec bestObjectiveIterate = x1;
  
  for (int i = 0; i < max_iters; ++i) {
    x0 = x1;
    arma::vec gradf0 = gradf1;
    double tau0 = tau1;
    
    arma::vec x1hat = x0 - tau0 * gradf0;
    x1 = as<arma::vec>(proxg(x1hat, tau0));
    arma::vec Dx = x1 - x0;
    
    d1 = x1;
    f1 = as<double>(f(d1));
    
    if (backtrack) {
      double M = fVals.subvec(std::max(i - w, 0), i - 1).max();
      backtrackCount = 0;
      
      while ((f1 - 1e-12 > M + arma::dot(Dx, gradf0) + 0.5 * arma::norm(Dx, 2) * arma::norm(Dx, 2) / tau0) &&
             (backtrackCount < 20)) {
        tau0 *= stepsizeShrink;
        x1hat = x0 - tau0 * gradf0;
        x1 = as<arma::vec>(proxg(x1hat, tau0));
        d1 = x1;
        f1 = as<double>(f(d1));
        Dx = x1 - x0;
        backtrackCount++;
      }
      totalBacktracks += backtrackCount;
    }
    
    taus(i) = tau0;
    residual(i) = arma::norm(Dx, 2) / tau0;
    maxResidual = std::max(maxResidual, residual(i));
    
    double normalizer = std::max(arma::norm(gradf0, 2), arma::norm(x1 - x1hat, 2) / tau0) + eps_n;
    normalizedResid(i) = residual(i) / normalizer;
    
    fVals(i) = f1;
    objective(i + 1) = f1 + as<double>(g(x1));
    
    if (objective(i + 1) < minObjectiveValue) {
      bestObjectiveIterate = x1;
      minObjectiveValue = std::min(minObjectiveValue, objective(i + 1));
    }
    
    gradf1 = as<arma::vec>(gradf(d1));
    arma::vec Dg = gradf1 + (x1hat - x0) / tau0;
    double dotprod = arma::dot(Dx, Dg);
    
    if (std::abs(dotprod) < 1e-15) {
      break;
    }
    
    double tau_s = arma::dot(Dx, Dx) / dotprod;
    double tau_m = dotprod / arma::dot(Dg, Dg);
    tau_m = std::max(tau_m, 0.0);
    
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


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
