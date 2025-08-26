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
        // accept only exact integers 0..ploidy (legacy semantics)
        const double vr = std::round(v);
        if (std::fabs(v - vr) < 1e-12) {
          const long dv = static_cast<long>(vr);
          if (dv >= 0 && dv <= static_cast<long>(ploidy)) {
            const uword cj = j * glev + static_cast<uword>(dv);
            out(i, cj) = 1.0;
          }
        }
      }
      // else: NA stays 0 (effectively treats NA as absence)
    }
  }
  return out;
}
