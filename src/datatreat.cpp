// Thiago de Paula Oliveira

#include <Rcpp.h>
using namespace Rcpp;
List ascii_to_number(CharacterMatrix pedigree_data, std::string unk = "0");

// [[Rcpp::export]]
List datatreat_cpp(CharacterMatrix data, int n_max = 50, std::string unk = "0", bool save = false) {
  int indicator = 0;
  IntegerVector k(2, NA_INTEGER);
  CharacterMatrix new_data = clone(data);
  int n = data.nrow();
  
  for (int i = 1; i <= n_max; ++i) {
    if (i > 1) {
      data = new_data;
    }
    
    List pedigree = ascii_to_number(data, unk);
    IntegerVector sire = pedigree["sire"];
    IntegerVector dire = pedigree["dire"];
    IntegerVector ind = seq(0, n - 1);
    IntegerVector right_pos(n, NA_INTEGER);
    
    IntegerVector parent = sire;
    std::string parent_ind = "sire";
    
    if (indicator % 2 == 1) {
      parent = dire;
      parent_ind = "dire";
    }
    
    for (int j = 0; j < n; ++j) {
      if (std::find(ind.begin(), ind.begin() + j + 1, parent[j]) == ind.begin() + j + 1 &&
          !IntegerVector::is_na(parent[j])) {
          auto it = std::find(ind.begin(), ind.end(), parent[j]);
        if (it != ind.end()) {
          right_pos[j] = std::distance(ind.begin(), it);
        }
      }
    }
    
    std::vector<int> error;
    for (int j = 0; j < n; ++j) {
      if (!IntegerVector::is_na(parent[j]) && parent[j] > j) {
        error.push_back(j);
      }
    }
    
    if (save) {
      Rcout << "iteration #" << i << " (" << parent_ind << ")\n";
      for (int e : error) Rcout << e << " ";
      Rcout << "\n";
    }
    
    if (!error.empty()) {
      std::vector<int> after, before;
      for (int j = 0; j < n; ++j) {
        if (!IntegerVector::is_na(right_pos[j])) {
          after.push_back(j);
          before.push_back(right_pos[j]);
        }
      }
      
      for (size_t j = 0; j < before.size(); ++j) {
        if (j == 0 || before[j] != before[j - 1]) {
          std::swap(ind[after[j]], ind[before[j]]);
        }
      }
    }
    
    CharacterMatrix temp_data(n, 3);
    for (int j = 0; j < n; ++j) {
      temp_data(j, 0) = data(ind[j], 0);
      temp_data(j, 1) = data(ind[j], 1);
      temp_data(j, 2) = data(ind[j], 2);
    }
    
    new_data = temp_data;
    int last_k = k[0];
    k[0] = error.size();
    
    if (k[0] == 0) indicator++;
    
    if (i != 1 && k[0] == 0 && last_k == 0) {
      Rcout << "Your data was chronologically organized with success.\n";
      if (save) {
        Rcpp::Function write_table("write.table");
        write_table(new_data, Named("file") = "orgped.txt", Named("quote") = false,
                    Named("row.names") = false, Named("col.names") = false);
      }
      return pedigree;
    }
    
    if (i == n_max) {
      Rcout << "Your data was not chronologically organized with success.\n";
      if (save) {
        Rcpp::Function write_table("write.table");
        write_table(new_data, Named("file") = "orgped.txt", Named("quote") = false,
                    Named("row.names") = false, Named("col.names") = false);
      }
      return pedigree;
    }
  }
  
  // Fallback in case loop exits unexpectedly
  stop("Unexpected termination of sorting loop.");
}

