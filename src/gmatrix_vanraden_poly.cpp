// Thiago de Paula Oliveira
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Internal: build centred Z (NA -> 0), and compute K
template <bool PL_CORR, bool RATIO>
  static void build_Z_and_K(const mat& X, const uword ploidy,
                            mat& Z, double& K) {
    const uword n = X.n_rows, m = X.n_cols;
    Z.set_size(n, m);
    K = 0.0;
    
    for (uword j = 0; j < m; ++j) {
      // mean, count over non-NA
      double sum = 0.0; uword cnt = 0;
      for (uword i = 0; i < n; ++i) {
        double v = X(i, j);
        if (!std::isnan(v)) { sum += v; ++cnt; }
      }
      const double mean = (cnt > 0) ? (sum / static_cast<double>(cnt)) : std::numeric_limits<double>::quiet_NaN();
      
      // centre, accumulate sumsq for variance when needed
      double sumsq = 0.0; uword cnt_var = 0;
      for (uword i = 0; i < n; ++i) {
        double v = X(i, j);
        if (!std::isnan(v)) {
          double z = v - mean;
          Z(i, j) = z;
          sumsq += z * z;
          ++cnt_var;
        } else {
          Z(i, j) = 0.0; // NA -> 0 for the cross-product
        }
      }
      
      if constexpr (PL_CORR) {
        if constexpr (RATIO) {
          // K += (1/ploidy) * mean * (1 - mean), skip if no data
          if (cnt > 0) K += (1.0 / static_cast<double>(ploidy)) * mean * (1.0 - mean);
        } else {
          // f = mean / ploidy; K += ploidy * f * (1 - f)
          if (cnt > 0) {
            double f = mean / static_cast<double>(ploidy);
            K += static_cast<double>(ploidy) * f * (1.0 - f);
          }
        }
      } else {
        // sample variance across rows
        if (cnt_var > 1) {
          double var_j = sumsq / static_cast<double>(cnt_var - 1);
          K += var_j;
        }
      }
    }
  }

// [[Rcpp::export]]
arma::mat Gmatrix_vanraden_poly_unweighted(const arma::mat& X,
                                           const unsigned int ploidy,
                                           const bool ratio,
                                           const bool ploidy_correction) {
  mat Z; double K = 0.0;
  if (ploidy_correction) {
    if (ratio) build_Z_and_K<true, true>(X, ploidy, Z, K);
    else       build_Z_and_K<true, false>(X, ploidy, Z, K);
  } else {
    if (ratio) build_Z_and_K<false, true>(X, ploidy, Z, K);
    else       build_Z_and_K<false, false>(X, ploidy, Z, K);
  }
  // G = (Z %*% t(Z)) / K
  mat G = (K > 0.0) ? (Z * Z.t()) / K : mat(Z.n_rows, Z.n_rows, fill::zeros);
  return G;
}

// [[Rcpp::export]]
arma::mat Gmatrix_vanraden_poly_weighted(const arma::mat& X,
                                        const arma::vec& w,
                                        const unsigned int ploidy,
                                        const bool ratio,
                                        const bool ploidy_correction) {
  if (X.n_cols != w.n_elem) stop("weights length does not match #markers");
  mat Z; double K = 0.0;
  if (ploidy_correction) {
    if (ratio) build_Z_and_K<true, true>(X, ploidy, Z, K);
    else       build_Z_and_K<true, false>(X, ploidy, Z, K);
  } else {
    if (ratio) build_Z_and_K<false, true>(X, ploidy, Z, K);
    else       build_Z_and_K<false, false>(X, ploidy, Z, K);
  }
  // diag(w): Zw = Z; Zw.each_row() %= w.t();
  mat Zw = Z;
  Zw.each_row() %= w.t();
  mat G = (K > 0.0) ? (Zw * Z.t()) / K : mat(Z.n_rows, Z.n_rows, fill::zeros);
  return G;
}