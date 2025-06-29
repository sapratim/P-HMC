// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Efficient phi function
inline arma::vec phi(double del, double x_i, const arma::vec &t) {
  arma::vec val = abs((x_i - t) / del);
  return (1.0 - val) % (val <= 1.0);
}

// Efficient target function
double target(const arma::vec &chi,
              const List &x,
              const arma::vec &c,
              const arma::vec &t,
              const arma::mat &inv_cov,
              const arma::vec &ns,
              double delta_m) {

  if(!all(chi> 0)) return R_NegInf ;
  int N0 = ns.n_elem;
  double ret = -0.5 * dot(chi, inv_cov * chi) - N0 * dot(c, chi);

  for (int j = 0; j < N0; ++j) {
    NumericVector temp = x[j];
    int nj = ns[j];
    for (int i = 0; i < nj; ++i) {
      arma::vec phi_vec = phi(delta_m, temp[i], t);
      ret += std::log(dot(phi_vec, chi));
    }
  }
  return ret;
}

// Efficient gradient of log-posterior
arma::vec grad_logpi_dur(const arma::vec &chi,
                         const List &x,
                         const arma::vec &c,
                         const arma::vec &t,
                         const arma::mat &inv_cov,
                         const arma::vec &ns,
                         double lambda,
                         double delta_m) {

  int p = chi.n_elem;
  int N0 = ns.n_elem;
  arma::vec track(p, fill::zeros);

  for (int j = 0; j < N0; ++j) {
    NumericVector temp = x[j];
    int nj = ns[j];

    for (int i = 0; i < nj; ++i) {
      arma::vec phi_vec = phi(delta_m, temp[i], t);
      double denom = dot(phi_vec, chi);
      track += phi_vec / denom;
    }
  }

  arma::vec prox = chi % (chi > 0);
  return -N0 * c + track - ((chi - prox) / lambda + inv_cov * chi);
}

// [[Rcpp::export]]
List cox_hmc_cpp(int N,
                 const arma::vec &init,
                 const arma::vec &ns,
                 const List &x,
                 const arma::vec &c,
                 const arma::vec &t,
                 const arma::mat &sqrt_cov,
                 const arma::mat &inv_cov,
                 double lambda,
                 double eps_hmc,
                 int L,
                 double delta_m) {

  int nvar = init.n_elem;
  arma::mat samp_hmc(N, nvar);
  arma::vec samp = init;
  samp_hmc.row(0) = samp.t();

  int accept = 0;
  arma::mat mom_mat = randn<arma::mat>(N, nvar);
  double  count = 0;
  double count2 = 0;
  for (int i = 1; i < N; ++i) {

    arma::vec p_prop = sqrt_cov * mom_mat.row(i).t();
    arma::vec U_samp = -grad_logpi_dur(samp, x, c, t, inv_cov, ns, lambda, delta_m);
    arma::vec p_current = p_prop - 0.5 * eps_hmc * U_samp;
    arma::vec q_current = samp;

    for (int j = 0; j < L; ++j) {
      samp += eps_hmc * inv_cov * p_current;
      if(!all(samp> 0)) count = count + 1; 
      U_samp = -grad_logpi_dur(samp, x, c, t, inv_cov, ns, lambda, delta_m);
      if (j != L - 1) {
        p_current -= eps_hmc * U_samp;
      }
    }

    p_current -= 0.5 * eps_hmc * U_samp;
    p_current = -p_current;

    if (all(samp >= 0)) {
      double U_curr = -target(q_current, x, c, t, inv_cov, ns, delta_m);
      double U_prop = -target(samp, x, c, t, inv_cov, ns, delta_m);
      double K_curr = 0.5 * dot(p_prop, inv_cov * p_prop);
      double K_prop = 0.5 * dot(p_current, inv_cov * p_current);
      double log_acc_prob = U_curr - U_prop + K_curr - K_prop;

      if (std::log(R::runif(0, 1)) <= log_acc_prob) {
        samp_hmc.row(i) = samp.t();
        accept++;
      } else {
        samp_hmc.row(i) = q_current.t();
        samp = q_current;
      }
    } else {
      count2 = count2 + 1;
      samp_hmc.row(i) = q_current.t();
      samp = q_current;
    }
    if (i % (N / 10) == 0) {
      Rcout << "Iter: " << i << "  Accept Rate: " << static_cast<double>(accept) / i << std::endl;
    }
  }
  Rcout << "outside:" << count/(L*N) << std::endl;
  Rcout << "outside final:" << count2/(N) << std::endl;
  double acc_rate = static_cast<double>(accept) / N;
  Rcout << "Final Acceptance Rate: " << acc_rate << std::endl;

  return List::create(
    Named("samples") = samp_hmc,
    Named("acceptance_rate") = acc_rate
  );
}



// Helper to check if all elements in vector are positive
inline bool is_positive(const arma::vec& v) {
  return all(v >= 0);
}

// [[Rcpp::export]]
List cox_bf_cpp(int N,
                const arma::vec& init,
                const arma::vec& ns,
                const List& x,
                const arma::vec& c,
                const arma::vec& t,
                const arma::mat& cov,
                double eta,
                const arma::mat& sqrt_cov,
                int m,
                double delta_m) {
  
  int p = init.n_elem;
  arma::mat chi(N, p, fill::zeros);
  arma::vec log_post(N, fill::zeros);
  arma::vec bernoulli_loops(N, fill::zeros);
  
  chi.row(0) = init.t();
  log_post(0) = target(init, x, c, t, cov, ns, delta_m);
  int accept_rate = 0;
  
  RNGScope scope; // ensure RNG works in Rcpp
  
  for(int i = 1; i < N; i++) {
    int accept = 0;
    
    if(i % (N/10) == 0) Rcout << "Iteration: " << i << std::endl;
    
    // Accept-Reject step
    arma::vec z(p);
    while(!accept) {
      z = chi.row(i-1).t() + sqrt(eta) * (sqrt_cov * arma::randn<arma::vec>(p));
      if(is_positive(z)) {
        accept = 1;
      }
    }
    arma::vec y = z;
    
    double c1 = target(y, x, c, t, cov, ns, delta_m);
    double c2 = target(chi.row(i-1).t(), x, c, t, cov, ns, delta_m);
    
    double exp_diff = std::exp(c1 - c2);
    double C = exp_diff / (1 + exp_diff);
    
    int bern_loops = 0;
    bool accepted = false;
    
    while(!accepted) {
      bern_loops++;
      
      int C1 = R::rbinom(1, C);
      
      if(C1 == 1) {
        arma::vec m1 = chi.row(i-1).t() + sqrt(eta) * (sqrt_cov * arma::randn<arma::vec>(m));
        int p_x = is_positive(m1) ? 1 : 0;
        int C2 = R::rbinom(1, p_x);
        if(C2 == 1) {
          chi.row(i) = y.t();
          log_post(i) = c1;
          accepted = true;
        }
      } else {
        arma::vec m2 = y + sqrt(eta) * (sqrt_cov * arma::randn<arma::vec>(m));
        int p_y = is_positive(m2) ? 1 : 0;
        int C2 = R::rbinom(1, p_y);
        if(C2 == 1) {
          chi.row(i) = chi.row(i-1);
          log_post(i) = c2;
          accepted = true;
        }
      }
    }
    bernoulli_loops(i) = bern_loops;
    
    if(chi(i, 0) == y(0)) accept_rate++;
  }
  
  double bern_loops_avg = mean(bernoulli_loops);
  double acc_rate = (double)accept_rate / N;
  
  return List::create(
    Named("chi") = chi,
    Named("bernoulli_loops") = bernoulli_loops,
    Named("accept_rate") = acc_rate,
    Named("log_post") = log_post
  );
}
