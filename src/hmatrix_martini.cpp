// Thiago de Paula Oliveira
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// SPD inverse (Cholesky when possible; fallback to pinv)
static mat inv_sympd_safe(const mat& X) {
  mat Y;
  // SPD path
  if (inv_sympd(Y, X)) return Y;   
  // Pseudo-inverse (numerically safer than inv for near-singular)
  return pinv(X);
}

//' Compute Martini block corrections (H11, H12, H21, H22 - A22)
//'
//' @param A12 n1 x n2 block from A
//' @param A22 n2 x n2 block from A (genotyped)
//' @param G22 n2 x n2 block from G (genotyped)
//' @param tau scalar
//' @param omega scalar
//' @return list(H11, H12, H21, H22corr) where H22corr = H22 - A22
// [[Rcpp::export]]
Rcpp::List H_martini_blocks(const arma::mat& A12,
                            const arma::mat& A22,
                            const arma::mat& G22,
                            const double tau,
                            const double omega) {
  const uword n1 = A12.n_rows;
  const uword n2 = A12.n_cols;
  
  if (A22.n_rows != n2 || A22.n_cols != n2)
    stop("A22 must be n2 x n2 where n2 = ncol(A12)");
  if (G22.n_rows != n2 || G22.n_cols != n2)
    stop("G22 must be n2 x n2 and align with A22");
  
  // Inverses via SPD path
  const mat A22inv = inv_sympd_safe(A22);
  const mat G22inv = inv_sympd_safe(G22);
  
  // S = tau * G22inv + (1 - omega) * A22inv ; H22 = inv(S)
  const mat S = tau * G22inv + (1.0 - omega) * A22inv;
  const mat H22 = inv_sympd_safe(S);
  
  // delta = H22 - A22
  const mat Delta = H22 - A22;
  
  // Precompute A12 * A22inv and A22inv * A12.t() to reuse
  const mat A12_A22inv   = A12 * A22inv;         // n1 x n2
  const mat A22inv_A21   = A22inv * A12.t();     // n2 x n1
  
  // Blocks
  // H12 = A12 * A22inv * delta
  const mat H12 = A12_A22inv * Delta;            // n1 x n2
  // H21 = delta * A22inv * A21
  const mat H21 = Delta * A22inv_A21;            // n2 x n1
  // H11 = A12 * A22inv * delta * A22inv * A21
  const mat H11 = A12_A22inv * Delta * A22inv_A21; // n1 x n1
  
  // H22corr = H22 - A22
  const mat H22corr = Delta;                     // n2 x n2
  
  return Rcpp::List::create(
    _["H11"]     = H11,
    _["H12"]     = H12,
    _["H21"]     = H21,
    _["H22corr"] = H22corr
  );
}
