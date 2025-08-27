// Thiago de Paula Oliveira
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix buildDominanceMatrix_cpp(NumericMatrix A, IntegerVector s, IntegerVector d) {
  int n = A.nrow();
  NumericMatrix D(n, n);
  
  for (int i = 0; i < n; ++i) {
    int si = s[i] - 1;
    int di = d[i] - 1;
    
    for (int j = i; j < n; ++j) {
      int sj = s[j] - 1;
      int dj = d[j] - 1;
      
      double u1 = (si >= 0 && sj >= 0) ? A(si, sj) : 0.0;
      double u2 = (di >= 0 && dj >= 0) ? A(di, dj) : 0.0;
      double u3 = (si >= 0 && dj >= 0) ? A(si, dj) : 0.0;
      double u4 = (sj >= 0 && di >= 0) ? A(sj, di) : 0.0;
      
      D(i, j) = D(j, i) = 0.25 * (u1 * u2 + u3 * u4);
    }
  }
  
  for (int i = 0; i < n; ++i) {
    D(i, i) = 1.0;
  }
  
  return D;
}