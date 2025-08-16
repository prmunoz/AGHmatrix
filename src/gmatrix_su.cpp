// Thiago de Paula Oliveira
#include <RcppArmadillo.h>
#include "gmatrix_utils.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat Gmatrix_Su(const arma::mat& M) {
  const uword n = M.n_rows, m = M.n_cols;
  
  vec p(m, fill::zeros);
  for (uword j = 0; j < m; ++j) p[j] = p2_allele_freq(M.col(j));
  
  vec TwoPQ = 2.0 * p % (1.0 - p);
  const double denom = accu(TwoPQ % (1.0 - TwoPQ));
  
  // build H where 1=heterozygote, 0 otherwise, then H - TwoPQ
  mat Z(n, m, fill::zeros);
  for (uword j = 0; j < m; ++j) {
    for (uword i = 0; i < n; ++i) {
      const double v = M(i,j);
      if (!std::isnan(v)) {
        const double h = (v == 1.0) ? 1.0 : 0.0;
        Z(i,j) = h - TwoPQ[j];
      }
    }
  }
  // only necessary if there were all-NA columns
  na_to_zero(Z); 
  
  return (Z * Z.t()) / denom;
}
