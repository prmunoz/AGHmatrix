// Thiago de Paula Oliveira
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp; using namespace arma;

// [[Rcpp::export]]
arma::mat Gmatrix_general_weighted(arma::mat Z,
                                  const arma::vec& weights,
                                  const double K) {
  if (weights.n_elem != Z.n_cols) stop("weights length mismatch");
  if (!(K > 0.0)) stop("K must be > 0");
  
  for (uword j = 0; j < Z.n_cols; ++j)
    for (uword i = 0; i < Z.n_rows; ++i)
      if (std::isnan(Z(i,j))) Z(i,j) = 0.0;
      
      mat W = diagmat(weights);
      return (Z * W * Z.t()) / K;
}