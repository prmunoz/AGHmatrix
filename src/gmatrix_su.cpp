// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// allele frequency p for a diploid-coded column (0/1/2)
static inline double diploid_alt_freq(const vec& col) {
  double c2 = 0.0, c1 = 0.0; uword nobs = 0;
  const uword n = col.n_elem;
  for (uword i = 0; i < n; ++i) {
    const double v = col[i];
    if (!std::isnan(v)) {
      if (v == 2.0)      c2 += 1.0;
      else if (v == 1.0) c1 += 1.0;
      // v==0 contributes only to denominator
      nobs += 1;
    }
  }
  if (nobs == 0) return NA_REAL; // no data at this marker
  return (2.0 * c2 + c1) / (2.0 * static_cast<double>(nobs));
}

// [[Rcpp::export]]
arma::mat Gmatrix_Su(const arma::mat& M) {
  const uword n = M.n_rows, m = M.n_cols;
  
  // 1) p_j and TwoPQ_j per marker
  vec p(m, fill::zeros);
  for (uword j = 0; j < m; ++j) {
    p[j] = diploid_alt_freq(M.col(j));
  }
  vec TwoPQ = 2.0 * p % (1.0 - p);
  
  // 2) Denominator sum_j TwoPQ_j * (1 - TwoPQ_j)
  const double denom = accu(TwoPQ % (1.0 - TwoPQ));
  
  // 3) Build Z = H - TwoPQ where H(i,j) = 1 if genotype==1, else 0; NA -> 0
  mat Z(n, m, fill::zeros);
  for (uword j = 0; j < m; ++j) {
    const double tpq = TwoPQ[j];
    for (uword i = 0; i < n; ++i) {
      const double v = M(i, j);
      if (!std::isnan(v)) {
        const double h = (v == 1.0) ? 1.0 : 0.0;
        Z(i, j) = h - tpq;
      } // else leave 0.0
    }
  }
  
  // 4) G = Z %*% t(Z) / denom
  return (Z * Z.t()) / denom;
}
