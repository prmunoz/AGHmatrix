// Thiago de Paula Oliveira
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "gmatrix_utils.h"
using namespace Rcpp; using namespace arma;

// [[Rcpp::export]]
arma::mat Gmatrix_unweighted(arma::mat Z, const double K) {
  if (!(K > 0.0)) stop("K must be > 0");
  
  // NA -> 0 (match R path)
  for (uword j = 0; j < Z.n_cols; ++j)
    for (uword i = 0; i < Z.n_rows; ++i)
      if (std::isnan(Z(i,j))) Z(i,j) = 0.0;
  
  return (Z * Z.t()) / K;
}
