// Thiago de Paula Oliveira
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void check_matrix(SEXP obj) {
  if (!Rf_isMatrix(obj)) {
    stop("Input must be a matrix.");
  }
}
