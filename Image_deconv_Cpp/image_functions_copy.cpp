// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// One-dimensional DWT function
void dwt(const vec& Vin, int M, int L, const vec& h, const vec& g, vec& Wout, vec& Vout) {
  int half_M = M / 2;
  Wout.zeros(half_M);
  Vout.zeros(half_M);
  
  for (int t = 0; t < half_M; ++t) {
    int u = 2 * t + 1;  // Start at odd indices
    Wout[t] = h[0] * Vin[u];
    Vout[t] = g[0] * Vin[u];
    for (int n = 1; n < L; ++n) {
      u -= 1;
      if (u < 0) u = M - 1;  // Periodic boundary
      Wout[t] += h[n] * Vin[u];
      Vout[t] += g[n] * Vin[u];
    }
  }
}

// Single-level 2D DWT function
Rcpp::List two_D_dwt_single(const mat& X, int M, int N, int L, const vec& h, const vec& g) {
  mat Low(M / 2, N, fill::zeros);
  mat High(M / 2, N, fill::zeros);
  
  // Step 1: DWT on columns
  vec Wout, Vout;
  for (int i = 0; i < N; ++i) {
    vec col = X.col(i);
    dwt(col, M, L, h, g, Wout, Vout);
    Low.col(i) = Vout;
    High.col(i) = Wout;
  }
  
  // Step 2: DWT on rows
  mat LL(M / 2, N / 2, fill::zeros);
  mat LH(M / 2, N / 2, fill::zeros);
  mat HL(M / 2, N / 2, fill::zeros);
  mat HH(M / 2, N / 2, fill::zeros);
  
  for (int i = 0; i < M / 2; ++i) {
    dwt(Low.row(i).t(), N, L, h, g, Wout, Vout);  // Transpose to row vector
    LL.row(i) = Vout.t();
    HL.row(i) = Wout.t();
    
    dwt(High.row(i).t(), N, L, h, g, Wout, Vout);
    LH.row(i) = Vout.t();
    HH.row(i) = Wout.t();
  }
  
  return Rcpp::List::create(
    Rcpp::Named("LL") = LL,
    Rcpp::Named("LH") = LH,
    Rcpp::Named("HL") = HL,
    Rcpp::Named("HH") = HH
  );
}

// [[Rcpp::export]]
Rcpp::List two_D_dwt_multi(const mat& X, int M, int N, int L, const vec& h, const vec& g, int J) {
  if (J < 1) Rcpp::stop("Number of levels must be at least 1.");
  int min_dim = std::min(M, N);
  int max_levels = 0;
  while (min_dim >= 2) { min_dim /= 2; max_levels++; }
  if (J > max_levels) Rcpp::stop("Requested levels exceed maximum possible for input dimensions.");
  
  Rcpp::List result;
  mat current_X = X;
  int current_M = M;
  int current_N = N;
  
  for (int j = 0; j < J; ++j) {
    Rcpp::List subbands = two_D_dwt_single(current_X, current_M, current_N, L, h, g);
    result["LH_" + std::to_string(j + 1)] = subbands["LH"];
    result["HL_" + std::to_string(j + 1)] = subbands["HL"];
    result["HH_" + std::to_string(j + 1)] = subbands["HH"];
    if (j == J - 1) result["LL"] = subbands["LL"];
    current_X = as<mat>(subbands["LL"]);
    current_M /= 2;
    current_N /= 2;
  }
  
  return result;
}

// One-dimensional IDWT function
void idwt(const vec& Win, const vec& Vin, int M, int L, const vec& h, const vec& g, vec& Xout) {
  Xout.zeros(2 * M);
  int m = -2, n = -1;
  
  for (int t = 0; t < M; ++t) {
    m += 2;
    n += 2;
    int u = t;
    int i = 1, j = 0;
    Xout[m] = h[i] * Win[u] + g[i] * Vin[u];
    Xout[n] = h[j] * Win[u] + g[j] * Vin[u];
    if (L > 2) {
      for (int l = 1; l < L / 2; ++l) {
        u = (u + 1) % M;  // Periodic boundary
        i += 2;
        j += 2;
        Xout[m] += h[i] * Win[u] + g[i] * Vin[u];
        Xout[n] += h[j] * Win[u] + g[j] * Vin[u];
      }
    }
  }
}

// [[Rcpp::export]]
mat two_D_idwt_single(const mat& LL, const mat& LH, const mat& HL, const mat& HH, 
                      int M, int N, int L, const vec& h, const vec& g) {
  // Step 1: Reconstruct rows (M x N -> M x 2N)
  mat Low(M, 2 * N, fill::zeros);    // M rows, 2N cols
  mat High(M, 2 * N, fill::zeros);   // M rows, 2N cols
  vec Win(N), Vin(N), Xout(2 * N);   // Xout has 2N elements
  
  for (int i = 0; i < M; ++i) {
    Win = HL.row(i).t();             // N elements
    Vin = LL.row(i).t();             // N elements
    idwt(Win, Vin, N, L, h, g, Xout); // Xout: 2N elements
    Low.row(i) = Xout.t();           // Assign to row i (2N cols)
    
    Win = HH.row(i).t();             // N elements
    Vin = LH.row(i).t();             // N elements
    idwt(Win, Vin, N, L, h, g, Xout); // Xout: 2N elements
    High.row(i) = Xout.t();          // Assign to row i (2N cols)
  }
  
  // Step 2: Reconstruct columns (M x 2N -> 2M x 2N)
  mat image(2 * M, 2 * N, fill::zeros);  // 2M rows, 2N cols
  vec Win2(M), Vin2(M), Xout2(2 * M);    // Xout2 has 2M elements
  
  for (int i = 0; i < 2 * N; ++i) {      // Iterate over all 2N columns
    Win2 = High.col(i);                  // M elements
    Vin2 = Low.col(i);                   // M elements
    idwt(Win2, Vin2, M, L, h, g, Xout2); // Xout2: 2M elements
    image.col(i) = Xout2;                // Assign 2M elements to column i
  }
  
  return image;
}

// [[Rcpp::export]]
mat two_D_idwt_multi(Rcpp::List y, int J, int L, const vec& h, const vec& g, int digits = 7) {
  if (J < 1) Rcpp::stop("Number of levels must be at least 1.");
  if (!y.containsElementNamed("LL")) Rcpp::stop("Input list 'y' does not contain 'LL'.");
  
  mat y_in = as<mat>(y["LL"]);
  int m = y_in.n_rows;
  int n = y_in.n_cols;
  
  for (int j = J; j >= 1; --j) {
    std::string LH_key = "LH_" + std::to_string(j);
    std::string HL_key = "HL_" + std::to_string(j);
    std::string HH_key = "HH_" + std::to_string(j);
    
    if (!y.containsElementNamed(LH_key.c_str()) || 
        !y.containsElementNamed(HL_key.c_str()) || 
        !y.containsElementNamed(HH_key.c_str())) {
        Rcpp::stop("Missing subband for level " + std::to_string(j));
    }
    
    mat LH = as<mat>(y[LH_key]);
    mat HL = as<mat>(y[HL_key]);
    mat HH = as<mat>(y[HH_key]);
    
    y_in = two_D_idwt_single(y_in, LH, HL, HH, m, n, L, h, g);
    m = y_in.n_rows;
    n = y_in.n_cols;
  }
  
  // Mimic R's zapsmall fully
  double mx = max(abs(y_in(find_finite(y_in))));
  if (mx == 0.0) return y_in;
  double tol = mx * std::pow(10.0, -digits);
  // Step 1: Zero out small values
  y_in.elem(find(abs(y_in) <= tol)).zeros();
  // Step 2: Round all finite values to 'digits' decimal places
  for (uword i = 0; i < y_in.n_elem; ++i) {
    if (std::isfinite(y_in[i])) {
      y_in[i] = std::round(y_in[i] * std::pow(10.0, digits)) / std::pow(10.0, digits);
    }
  }
  
  return y_in;
}

// [[Rcpp::export]]
double wavelet_l1_cpp(const vec& image_vec, int dimen, const vec& h, const vec& g, int nlev = 3) {
  mat image_mat = reshape(image_vec, dimen, dimen);  // Armadillo reshape
  Rcpp::List trans = two_D_dwt_multi(image_mat, dimen, dimen, h.n_elem, h, g, nlev);
  
  double wave_sum = 0.0;
  Rcpp::CharacterVector names = trans.names();
  for (int i = 0; i < trans.size(); ++i) {
    wave_sum += sum(sum(abs(as<mat>(trans[i]))));
  }
  
  return wave_sum;
}

// [[Rcpp::export]]
vec softthreshold(const vec& u, double pen) {
  vec result = sign(u) % max(abs(u) - pen, zeros<vec>(u.n_elem));
  return result;
}

// [[Rcpp::export]]
vec convolve_image(const vec& image, int img_rows, int img_cols, const mat& filter) {
  int filter_rows = filter.n_rows;
  int filter_cols = filter.n_cols;
  int pad_rows = filter_rows / 2;
  int pad_cols = filter_cols / 2;
  
  vec output(image.n_elem, fill::zeros);
  mat img_mat = reshape(image, img_cols, img_rows).t();
  
  for (int i = 0; i < img_rows; ++i) {
    for (int j = 0; j < img_cols; ++j) {
      double sum = 0.0;
      for (int fi = 0; fi < filter_rows; ++fi) {
        for (int fj = 0; fj < filter_cols; ++fj) {
          int ni = std::max(0, std::min(i + fi - pad_rows, img_rows - 1));
          int nj = std::max(0, std::min(j + fj - pad_cols, img_cols - 1));
          sum += img_mat(ni, nj) * filter(fi, fj);
        }
      }
      output[i * img_cols + j] = sum;
    }
  }
  
  return output;
}

// [[Rcpp::export]]
double f(const vec& z, int dimen, const mat& H, const vec& y, double sigma2, 
         const vec& x_true, double lamb) {
  vec H_z = convolve_image(z, dimen, dimen, H);
  return accu(square(H_z - y)) / (2 * sigma2) + accu(square(x_true - z)) / (2 * lamb);
}

// [[Rcpp::export]]
vec gradf(const vec& z, int dimen, const mat& H, const vec& y, double sigma2, 
          const vec& x_true, double lamb) {
  vec H_z = convolve_image(z, dimen, dimen, H);
  mat H_t = H.t();
  vec t1 = convolve_image(H_z - y, dimen, dimen, H_t) / sigma2;
  vec t2 = (x_true - z) / lamb;
  return t1 + t2;
}

// [[Rcpp::export]]
double g(const vec& z, int dimen, const vec& h, const vec& g_wave, double beta_pen) {
  return beta_pen * wavelet_l1_cpp(z, dimen, h, g_wave, 3);
}

// [[Rcpp::export]]
vec proxg(const vec& z_vec, int dimen, const vec& h, const vec& g_wave, double tau_fasta, 
          double beta_pen) {
  mat z = reshape(z_vec, dimen, dimen);
  Rcpp::List wave_trans = two_D_dwt_multi(z, dimen, dimen, h.n_elem, h, g_wave, 3);
  
  Rcpp::CharacterVector names = wave_trans.names();
  for (int i = 0; i < wave_trans.size(); ++i) {
    std::string name = as<std::string>(names[i]);
    if (name != "LL") {
      mat subband = as<mat>(wave_trans[i]);
      vec subband_vec = vectorise(subband);
      vec thresh_vec = softthreshold(subband_vec, beta_pen * tau_fasta);
      wave_trans[i] = reshape(thresh_vec, subband.n_rows, subband.n_cols);
    }
  }
  
  return vectorise(two_D_idwt_multi(wave_trans, 3, h.n_elem, h, g_wave));
}

// [[Rcpp::export]]
Rcpp::List fasta_cpp(const vec& x0, double tau1, const mat& H, const vec& y,
                     const vec& h, const vec& g_wave, double sigma2, const vec& x_true,
                     double lamb, double beta_pen, int dimen, int nlev = 3, int max_iters = 100,
                     int w = 10, bool backtrack = true, bool recordIterates = false,
                     double stepsizeShrink = 0.5, double eps_n = 1e-15) {
  
  const int n = x0.n_elem;
  vec residual(max_iters, fill::zeros);
  vec normalizedResid(max_iters, fill::zeros);
  vec taus(max_iters, fill::zeros);
  vec fVals(max_iters, fill::zeros);
  vec objective(max_iters + 1, fill::zeros);
  mat iterates(recordIterates ? n : 0, recordIterates ? max_iters + 1 : 0, fill::zeros);
  
  if (recordIterates) iterates.col(0) = x0;
  
  int totalBacktracks = 0;
  double maxResidual = -datum::inf;
  double minObjectiveValue = datum::inf;
  
  vec x1 = x0;
  vec d1(n);
  vec x0_old(n);
  vec gradf0(n);
  vec x1hat(n);
  vec Dx(n);
  vec Dg(n);
  vec gradf1(n);
  vec bestObjectiveIterate = x1;
  
  d1 = x1;
  double f1 = f(d1, dimen, H, y, sigma2, x_true, lamb);
  fVals[0] = f1;
  gradf1 = gradf(d1, dimen, H, y, sigma2, x_true, lamb);
  objective[0] = f1 + g(x1, dimen, h, g_wave, beta_pen);
  minObjectiveValue = objective[0];
  
  int i;
  for (i = 0; i < max_iters; ++i) {
    x0_old = x1;
    gradf0 = gradf1;
    double tau0 = tau1;
    
    x1hat = x0_old - tau0 * gradf0;
    x1 = proxg(x1hat, dimen, h, g_wave, tau0, beta_pen);
    Dx = x1 - x0_old;
    d1 = x1;
    f1 = f(d1, dimen, H, y, sigma2, x_true, lamb);
    
    double Dx_Dx = dot(Dx, Dx);
    
    if (backtrack) {
      int startIdx = std::max(0, i - w);
      int endIdx = i;  // Include current fVals[i]
      double M = (i == 0) ? fVals[0] : fVals.subvec(startIdx, endIdx).max();
      int backtrackCount = 0;
      double Dx_gradf0 = dot(Dx, gradf0);
      while ((f1 - 1e-12 > M + Dx_gradf0 + 0.5 * Dx_Dx / tau0) && (backtrackCount < 20)) {
        tau0 *= stepsizeShrink;
        x1hat = x0_old - tau0 * gradf0;
        x1 = proxg(x1hat, dimen, h, g_wave, tau0, beta_pen);
        d1 = x1;
        f1 = f(d1, dimen, H, y, sigma2, x_true, lamb);
        Dx = x1 - x0_old;
        Dx_gradf0 = dot(Dx, gradf0);
        Dx_Dx = dot(Dx, Dx);
        ++backtrackCount;
      }
      
      totalBacktracks += backtrackCount;
    }
    
    taus[i] = tau0;
    residual[i] = sqrt(Dx_Dx) / tau0;
    maxResidual = std::max(maxResidual, residual[i]);
    
    double normalizer = std::max(norm(gradf0), norm(x1 - x1hat) / tau0) + eps_n;
    normalizedResid[i] = residual[i] / normalizer;
    fVals[i] = f1;
    objective[i + 1] = f1 + g(x1, dimen, h, g_wave, beta_pen);
    
    double newObjectiveValue = objective[i + 1];
    if (newObjectiveValue < minObjectiveValue) {
      bestObjectiveIterate = x1;
      minObjectiveValue = newObjectiveValue;
    }
    
    gradf1 = gradf(d1, dimen, H, y, sigma2, x_true, lamb);
    Dg = gradf1 + (x1hat - x0_old) / tau0;
    double dotprod = dot(Dx, Dg);
    //Rcout << "The value of Dx_Dx : " << std::abs(dotprod) << "\n";
    if (std::abs(dotprod) < eps_n) break;
    
    double tau_s = Dx_Dx / dotprod;
    double tau_m = std::max(dotprod / dot(Dg, Dg), 0.0);
    tau1 = (2 * tau_m > tau_s) ? tau_m : tau_s - 0.5 * tau_m;
    if (tau1 <= 0 || std::isinf(tau1) || std::isnan(tau1)) tau1 = tau0 * 1.5;
    
    if (recordIterates) iterates.col(i + 1) = x1;
  }
  
  int final_iter = (i < max_iters) ? i : (max_iters - 1); // Cap final_iter
  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("x") = bestObjectiveIterate,
    Rcpp::Named("objective") = objective.subvec(0, final_iter + 1), // Include initial + iterations
    Rcpp::Named("fVals") = fVals.subvec(0, final_iter),
    Rcpp::Named("totalBacktracks") = totalBacktracks,
    Rcpp::Named("residual") = residual.subvec(0, final_iter),
    Rcpp::Named("taus") = taus.subvec(0, final_iter),
    Rcpp::Named("iternumber") = final_iter + 1
  );
  
  if (recordIterates) result["iterates"] = iterates.cols(0, final_iter + 1);
  
  return result;
}