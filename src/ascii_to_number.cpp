// Thiago de Paula Oliveira
// ascii_to_number.cpp
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List ascii_to_number(CharacterMatrix pedigree_data, std::string unk = "0") {
  int n = pedigree_data.nrow();
  if (pedigree_data.ncol() != 3) {
    Rcpp::stop("Data must have exactly 3 columns");
  }
  
  // Check for duplicate IDs
  std::set<std::string> seen;
  for (int i = 0; i < n; ++i) {
    std::string id = Rcpp::as<std::string>(pedigree_data(i, 0));
    if (!seen.insert(id).second) {
      Rcpp::stop("Duplicate individual ID found: " + id);
    }
  }
  
  // Create full list of names with unk first (like legacy R)
  CharacterVector ind_data(n + 1);
  ind_data[0] = unk;
  for (int i = 0; i < n; ++i) {
    ind_data[i + 1] = pedigree_data(i, 0);
  }
  
  CharacterVector sire_data = pedigree_data(_, 1);
  CharacterVector dire_data = pedigree_data(_, 2);
  IntegerVector sire(n), dire(n);
  
  for (int i = 0; i < n; ++i) {
    auto it_sire = std::find(ind_data.begin(), ind_data.end(), sire_data[i]);
    auto it_dire = std::find(ind_data.begin(), ind_data.end(), dire_data[i]);
    
    sire[i] = (it_sire != ind_data.end()) ? std::distance(ind_data.begin(), it_sire) : NA_INTEGER;
    dire[i] = (it_dire != ind_data.end()) ? std::distance(ind_data.begin(), it_dire) : NA_INTEGER;
  }
  
  // Return names of individuals (no 'unk')
  CharacterVector ind_names = pedigree_data(_, 0);
  
  return List::create(
    _["sire"] = sire,
    _["dire"] = dire,
    _["ind_data"] = ind_names
  );
}
