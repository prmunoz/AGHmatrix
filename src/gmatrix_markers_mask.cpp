// Thiago de Paula Oliveira
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp; using namespace arma;

// [[Rcpp::export]]
arma::mat Gmatrix_MarkersMask(const arma::mat& M) {
  const uword n = M.n_rows, m = M.n_cols;
  // Build a binary mask B where B(i,j)=1 if !is.na(M(i,j)), else 0
  mat B(n, m, fill::zeros);
  for (uword j = 0; j < m; ++j)
    for (uword i = 0; i < n; ++i)
      if (!std::isnan(M(i,j))) B(i,j) = 1.0;
  
  return B * B.t();  // tcrossprod(mask, mask)
}
