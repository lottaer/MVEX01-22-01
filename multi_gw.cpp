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

// New simulation - we remember the ages of cells

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
        // The cell dies immediately
        else {
          ages[0]++;
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

// This is the old simulation - does not remember age

void Z_def(int i, double p, int n, int k, double q, NumericVector Z) {
  double U = R::runif(0,1);
  int def_m = i + n;
  if(U < q) {
    int def = R::rbinom(def_m, p);
    def_m -= def;
    // daughter cell
    if(def <= k){
      Z[def]++;
    }
  }
  if(def_m <= k) {
    Z[def_m]++;
  }
}

NumericVector Z_t(double p, int n, int k, double q, NumericVector curr, NumericVector nxt) {
  for(int i = 0; i <= k; i++) {
    for(int j = 0; j < curr[i]; j++ ) {
      Z_def(i, p, n, k, q, nxt);
    }
  }
  return(nxt);
}

//[[Rcpp::export]]
NumericMatrix Z_mat(double p, int n, int k, double q, int hours, NumericVector start) {
  // gör nåt 
  NumericMatrix Z_m (hours+1, k+1);
  Z_m(0, _) = start;
  for(int i = 1; i <= hours; i++) {
    Z_m(i,_) = Z_t(p, n, k, q, Z_m(i-1,_), Z_m(i,_));
  }
  return(Z_m);
}

// This function constructs the rows of the matrix
NumericVector M_row_cpp(int i, double p, int n, int k, double q, NumericVector Z_r) {
  Z_def(i, p, n, k ,q, Z_r);
  return(Z_r);
}


// This functions builds the M matrix from independent rows
//[[Rcpp::export]]
NumericMatrix M_est_cpp(double p, int n, int k, double q) {
  NumericMatrix M (k+1, k+1);
  for(int i = 0; i <= k; i++) {
    M(i,_) = M_row_cpp(i, p, n, k, q, M(i, _));
  }
  return(M);
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
