// Thiago de Paula Oliveira
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List preprocessSNPmatrix(NumericMatrix mat, double missingValue) {
  int nrow = mat.nrow(), ncol = mat.ncol();
  NumericMatrix freq(ncol, 2);
  
  for (int j = 0; j < ncol; ++j) {
    double ref = 0, alt = 0, valid = 0;
    for (int i = 0; i < nrow; ++i) {
      if (mat(i, j) == missingValue) {
        mat(i, j) = NA_REAL;
      } else {
        double val = mat(i, j);
        if (!R_IsNA(val)) {
          if (val == 0) ref += 2;
          else if (val == 1) { ref += 1; alt += 1; }
          else if (val == 2) alt += 2;
          valid += 1;
        }
      }
    }
    freq(j, 0) = ref / (2 * valid);
    freq(j, 1) = alt / (2 * valid);
  }
  
  return List::create(Named("SNPmatrix") = mat,
                      Named("Frequency") = freq);
}
