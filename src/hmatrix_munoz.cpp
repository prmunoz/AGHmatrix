// Thiago de Paula Oliveira
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <unordered_map>
#include <vector>
#include <cmath>
using namespace Rcpp;
using namespace arma;

struct Agg {
  double sum  = 0.0;
  double sum2 = 0.0;
  std::size_t n = 0;
};

//' Compute var/mean of G by classes of rounded A, and fill varA/meanG
//'
//' @param G square matrix
//' @param A square matrix (same dim as G)
//' @param round_digits integer number of decimal places to round A
//' @return list(varA, meanG, classes, var_by_class, mean_by_class)
// [[Rcpp::export]]
Rcpp::List munoz_var_mean_by_Aclass(const arma::mat& G,
                                    const arma::mat& A,
                                    const int round_digits) {
  const uword n = G.n_rows;
  if (G.n_cols != n || A.n_rows != n || A.n_cols != n)
    stop("G and A must be square and of the same dimensions");
  
  const double pow10 = std::pow(10.0, static_cast<double>(std::max(0, round_digits)));
  
  // Map rounded A value -> class id
  std::unordered_map<double, int> cls_map;
  std::vector<double> cls_vals; 
  
  // Build class ids and aggregates
  std::vector<Agg> agg; agg.reserve(64);
  
  auto get_class_id = [&](double aij) -> int {
    // Round to round_digits
    const double rounded = std::round(aij * pow10) / pow10;
    auto it = cls_map.find(rounded);
    if (it == cls_map.end()) {
      int id = static_cast<int>(cls_vals.size());
      cls_map.emplace(rounded, id);
      cls_vals.push_back(rounded);
      agg.push_back(Agg{});
      return id;
    } else {
      return it->second;
    }
  };
  
  // We also keep the class id per cell to avoid recomputation in the 
  // second fill step
  arma::Mat<int> class_id(n, n);
  class_id.fill(-1);
  
  for (uword i = 0; i < n; ++i) {
    for (uword j = 0; j < n; ++j) {
      const double aij = A(i, j);
      // skip NA in classing
      if (std::isnan(aij)) continue;  
      const int cid = get_class_id(aij);
      class_id(i, j) = cid;
      
      const double gij = G(i, j);
      // skip NA in G for stats
      if (std::isnan(gij)) continue;
      
      Agg &ag = agg[cid];
      ag.sum  += gij;
      ag.sum2 += gij * gij;
      ag.n    += 1u;
    }
  }
  
  // Compute per-class mean and unbiased variance
  const int C = static_cast<int>(agg.size());
  arma::vec mean_by_class(C, fill::value(NAN));
  arma::vec var_by_class (C, fill::value(NAN));
  for (int c = 0; c < C; ++c) {
    const Agg &ag = agg[c];
    if (ag.n >= 1u) {
      const double mu = ag.sum / static_cast<double>(ag.n);
      mean_by_class[c] = mu;
      if (ag.n >= 2u) {
        const double ss = ag.sum2 - static_cast<double>(ag.n) * mu * mu;
        // sample variance
        var_by_class[c]  = ss / static_cast<double>(ag.n - 1u);  
      }
    }
  }
  
  // Fill matrices varA and meanG using class_id
  arma::mat varA(n, n, fill::value(NAN));
  arma::mat meanG(n, n, fill::value(NAN));
  for (uword i = 0; i < n; ++i) {
    for (uword j = 0; j < n; ++j) {
      const int cid = class_id(i, j);
      if (cid >= 0) {
        varA (i, j) = var_by_class[cid];
        meanG(i, j) = mean_by_class[cid];
      }
    }
  }
  
  // Return classes in ascending numeric order
  // We will sort and reorder the per-class vectors.
  std::vector<int> order(C);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(),
            [&](int a, int b){ return cls_vals[a] < cls_vals[b]; });

  arma::vec classes(C);
  arma::vec mean_sorted(C);
  arma::vec var_sorted(C);
  for (int k = 0; k < C; ++k) {
    const int c = order[k];
    classes[k]     = cls_vals[c];
    mean_sorted[k] = mean_by_class[c];
    var_sorted[k]  = var_by_class[c];
  }

  return Rcpp::List::create(
    _["varA"]          = varA,
    _["meanG"]         = meanG,
    _["classes"]       = classes,
    _["mean_by_class"] = mean_sorted,
    _["var_by_class"]  = var_sorted
  );
}
