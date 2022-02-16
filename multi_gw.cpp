#include <Rcpp.h> 
#include <list>
#include <iostream>
#include <iterator>

using namespace Rcpp;

// replicative age?
struct cell {
  int type;
  int gen;
};

//[[Rcpp::export]]
NumericVector age_prop(double p, int n, int k, double q, int hours, int st_type){
  // some starting population
  NumericVector ages (hours+1);
  
  std::list<cell> Z;
  Z.push_front({ st_type, 0 });
  
  for(int j = 1; j <= hours; j++) {
    for(auto it = Z.begin(); it != Z.end(); ) {
      double U = R::runif(0,1);
      int def_m = it->type + n;
      // the cell divides
      if(U < q) {
        int def = R::rbinom(def_m, p);
        def_m -= def;
        if(def <= k) {
          Z.push_front({ def, j });
        }
      }
      if(def_m <= k){
        it->type = def_m;
        ++it;
      }
      else {
        ages[j-(it->gen)]++;
        it = Z.erase(it);
      }
    }
  }
  
  return ages;
}

/**
// [[Rcpp::export]]
NumericVector a_row(int i, double p, int n, int k) {
  NumericVector row (k);
  for (int j = 0; j <= k; j++)
      row[j] = R::dbinom(i+n-j, n+i, p,0) + R::dbinom(j, n+i, p,0); 
  return row;
}
// [[Rcpp::export]]
NumericVector b_row(int i, int n, int k, double q) {
  NumericVector b (k);
  if(i+n <= k) {
    b[i+n] = 1-q;
  }
  return b;
}

//[[Rcpp::export]]
NumericMatrix mat(double p, int n, int k, double q, int gen) {
  NumericMatrix M (k+1);
  for(int i = 0; i <= k; i++){
    M(i,_) = q*a_row(i, p, n, k) + b_row(i, n, k, q);
  }
  return(M);
}
 **/