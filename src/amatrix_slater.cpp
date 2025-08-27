// Thiago de Paula Oliveira
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix buildA_slater_cpp(IntegerVector s, IntegerVector d, double w) {
  int n = s.size();
  NumericMatrix A(n, n);
  A(0, 0) = (1.0 + w) / 4.0;
  
  for (int i = 1; i < n; ++i) {
    int si = s[i] - 1;
    int di = d[i] - 1;
    
    if (s[i] == 0 && d[i] == 0) {
      A(i, i) = (1.0 + w) / 4.0;
    } else if (s[i] == 0) {
      A(i, i) = (5 + 7*w + 4*A(di, di)*(1-w)) / 24.0;
      for (int j = 0; j < i; ++j) {
        A(j, i) = A(i, j) = 0.5 * A(j, di);
      }
    } else if (d[i] == 0) {
      A(i, i) = (5 + 8*w + 4*A(si, si)*(1-w)) / 24.0;
      for (int j = 0; j < i; ++j) {
        A(j, i) = A(i, j) = 0.5 * A(j, si);
      }
    } else {
      A(i, i) = (1 + 2*w + (1-w)*(A(si, si) + A(di, di)) + 3*A(si, di)) / 6.0;
      for (int j = 0; j < i; ++j) {
        A(j, i) = A(i, j) = 0.5 * (A(j, si) + A(j, di));
      }
    }
  }
  
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      A(i, j) *= 4.0;
  
  return A;
}