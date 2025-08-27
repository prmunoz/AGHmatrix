// Thiago de Paula Oliveira
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "gmatrix_utils.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat Gmatrix_Vitezica(arma::mat M, const arma::mat& FreqP) {
  const uword n = M.n_rows, m = M.n_cols;
  if (FreqP.n_rows != n || FreqP.n_cols != m)
    stop("FreqP must be n x m matrix (replicated column-wise frequencies).");
  
  rowvec F = FreqP.row(0); // vector of FreqP per marker
  vec TwoPQ = 2.0 * F.t() % (1.0 - F.t());
  const double denom = accu(square(TwoPQ));
  
  // replace NA by 0 BEFORE recoding
  na_to_zero(M);
  
  mat X(n, m, fill::zeros);
  for (uword j = 0; j < m; ++j) {
    const double p = F[j];
    const double term0 = -2.0 * p * p;
    const double term1 =  2.0 * p * (1.0 - p);
    const double term2 = -2.0 * (1.0 - p) * (1.0 - p);
    for (uword i = 0; i < n; ++i) {
      const double v = M(i,j);
      if (v == 0.0)      X(i,j) = term0;
      else if (v == 1.0) X(i,j) = term1;
      else if (v == 2.0) X(i,j) = term2;
      else               X(i,j) = 0.0; // safety for unusual codes
    }
  }
  
  return (X * X.t()) / denom;
}
