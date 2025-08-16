// Thiago de Paula Oliveira
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix buildA_ploidy2_cpp(IntegerVector sire, IntegerVector dire, int n) {
  NumericMatrix A(n, n);
  
  A(0, 0) = 1.0;
  
  for (int i = 1; i < n; ++i) {
    int si = sire[i] - 1; // 1-based to 0-based
    int di = dire[i] - 1;
    
    if (sire[i] == 0 && dire[i] == 0) {
      A(i, i) = 1.0;
    } else if (sire[i] == 0) {
      A(i, i) = 1.0;
      for (int j = 0; j < i; ++j) {
        A(j, i) = A(i, j) = 0.5 * A(j, di);
      }
    } else if (dire[i] == 0) {
      A(i, i) = 1.0;
      for (int j = 0; j < i; ++j) {
        A(j, i) = A(i, j) = 0.5 * A(j, si);
      }
    } else {
      A(i, i) = 1.0 + 0.5 * A(di, si);
      for (int j = 0; j < i; ++j) {
        A(j, i) = A(i, j) = 0.5 * (A(j, si) + A(j, di));
      }
    }
  }
  
  return A;
}
