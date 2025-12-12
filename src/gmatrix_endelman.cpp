// Thiago de Paula Oliveira
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "gmatrix_utils.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat Gmatrix_Endelman(const arma::mat& M, int ploidy) {
  if (ploidy != 4) stop("'Endelman' method is just implemented for ploidy=4");
  const uword n = M.n_rows, m = M.n_cols;
  
  vec freq = mean(M, 0).t() / 4.0; // colMeans / ploidy
  vec freq_comp = 1.0 - freq;
  const double SixPQ = 6.0 * dot(square(freq), square(freq_comp));
  
  mat SNP = M;
  rowvec f  = freq.t();
  rowvec f2 = square(f);
  
  mat rep1 = repmat(f2, n, 1);      // (Frequency[,1]^2) repeated by rows
  mat rep2 = repmat(f,  n, 1);      // Frequency[,1]
  
  SNP = 6.0 * rep1 - 3.0 * rep2 % M + 0.5 * M % (M - 1.0);
  
  return (SNP * SNP.t()) / SixPQ;
}
