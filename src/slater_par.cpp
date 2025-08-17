// Thiago de Paula Oliveira
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat slater_par_cpp(const arma::mat& X, const unsigned int ploidy){
  const uword n = X.n_rows, m = X.n_cols;
  const uword glev = static_cast<uword>(ploidy) + 1u;
  mat out(n, m * glev, fill::zeros);
  
  for (uword j = 0; j < m; ++j){
    const double* col = X.colptr(j);
    for (uword i = 0; i < n; ++i){
      const double v = col[i];
      if (!std::isnan(v)) {
        // Accept only valid discrete dosages 0 ... ploidy
        const long dv = static_cast<long>(v);
        if (dv >= 0 && dv <= static_cast<long>(ploidy)){
          // column index for genotype level dv
          const uword cj = j * glev + static_cast<uword>(dv);
          out(i, cj) = 1.0;
        }
      }
      // else: NA stays 0 (effectively treats NA as absence)
    }
  }
  return out;
}
