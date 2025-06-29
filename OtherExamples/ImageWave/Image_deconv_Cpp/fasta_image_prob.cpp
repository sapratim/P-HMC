// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Convolution function with mirror padding
// [[Rcpp::export]]
NumericVector convolve_image(NumericVector image, int img_rows, 
                             int img_cols, NumericMatrix filter) {
  int filter_rows = filter.nrow();
  int filter_cols = filter.ncol();
  int pad_rows = filter_rows / 2;
  int pad_cols = filter_cols / 2;
  
  NumericVector output(image.size(), 0.0);
  
  for (int i = 0; i < img_rows; i++) {
    for (int j = 0; j < img_cols; j++) {
      double sum = 0.0;
      
      for (int fi = 0; fi < filter_rows; fi++) {
        for (int fj = 0; fj < filter_cols; fj++) {
          int ni = i + fi - pad_rows;
          int nj = j + fj - pad_cols;
          
          // Mirror padding instead of zero-padding
          ni = std::max(0, std::min(ni, img_rows - 1));
          nj = std::max(0, std::min(nj, img_cols - 1));
          
          sum += image[ni * img_cols + nj] * filter(fi, fj);
        }
      }
      
      output[i * img_cols + j] = sum;
    }
  }
  
  return output;
}
