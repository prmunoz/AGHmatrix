// Thiago de Paula Oliveira
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "gmatrix_utils.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat Gmatrix_Slater(const arma::mat& M,          // n x p (after slater_par in R)
                         const arma::vec& Frequency,  // length p (colMeans of M)
                         const arma::mat& FreqP,      // n x p (can be ignored/rebuilt)
                         const int NumberMarkers)     // legacy denominator (pre-drop)
{
  const uword n = M.n_rows, p = M.n_cols;
  if (Frequency.n_elem != p) stop("Frequency length mismatch");
  if (NumberMarkers <= 0)    stop("NumberMarkers must be > 0");
  
  // --- drop only Frequency == 0 (legacy) ---
  uvec keep = find(Frequency != 0.0);
  mat   M_use    = M.cols(keep);
  vec   F_use    = Frequency.elem(keep);
  rowvec F_row   = F_use.t();
  rowvec Fpq_row = F_row % (1.0 - F_row);       // may contain zeros (for F==1)
  
  // Build FreqP/FreqPQ for kept columns
  mat FreqP_use  = repmat(F_row, n, 1);
  mat FreqPQ_use = repmat(Fpq_row, n, 1);
  
  // --- Z standardisation; set NaN/Inf -> 0 (legacy) ---
  mat Z = (M_use - FreqP_use) / sqrt(FreqPQ_use);
  na_to_zero(Z);
  
  // --- cross-product divided by NumberMarkers (legacy denominator) ---
  mat G = (Z * Z.t()) / static_cast<double>(NumberMarkers);
  
  // --- diagonal terms; zero-out elementwise NaN before summing (legacy) ---
  mat Gii_mat = (square(M_use) - 2.0 * (FreqP_use % M_use) + square(FreqP_use)) / FreqPQ_use;
  na_to_zero(Gii_mat);  // IMPORTANT: elementwise, before sum
  vec Gii = 1.0 + sum(Gii_mat, 1) / static_cast<double>(NumberMarkers);
  
  G.diag() = Gii;
  return G;
}