// Thiago de Paula Oliveira
// ascii_to_number.cpp
#include <Rcpp.h>
#include <unordered_map>
using namespace Rcpp;

// [[Rcpp::export]]
List ascii_to_number(CharacterMatrix pedigree_data, std::string unk = "0") {
  int n = pedigree_data.nrow();
  if (pedigree_data.ncol() != 3) {
    Rcpp::stop("Data must have exactly 3 columns");
  }
  
  // Map to check duplicates and index individuals
  std::unordered_map<std::string, int> id_map;
  id_map[unk] = 0;
  
  CharacterVector ind_names(n);
  for (int i = 0; i < n; ++i) {
    std::string id = Rcpp::as<std::string>(pedigree_data(i, 0));
    if (id_map.find(id) != id_map.end()) {
      Rcpp::stop("Duplicate individual ID found: " + id);
    }
    id_map[id] = i + 1;
    ind_names[i] = id;
  }
  
  IntegerVector sire(n), dire(n);
  
  for (int i = 0; i < n; ++i) {
    std::string sire_id = Rcpp::as<std::string>(pedigree_data(i, 1));
    std::string dire_id = Rcpp::as<std::string>(pedigree_data(i, 2));
    
    auto it_sire = id_map.find(sire_id);
    auto it_dire = id_map.find(dire_id);
    
    sire[i] = (it_sire != id_map.end()) ? it_sire->second : NA_INTEGER;
    dire[i] = (it_dire != id_map.end()) ? it_dire->second : NA_INTEGER;
  }
  
  return List::create(
    _["sire"] = sire,
    _["dire"] = dire,
    _["ind_data"] = ind_names
  );
}
