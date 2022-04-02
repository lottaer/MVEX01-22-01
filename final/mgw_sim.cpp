#include <Rcpp.h> 
#include <list>
#include <iostream>
#include <iterator>

using namespace Rcpp;

// For age distribution
struct cell {
  int type;   // keep track of biological age over time
  int age;    // to calculate age at death
};

// ------------------- POPULATIONS --------------------------

// applied to each cell in the population
void Z_def(int i, double p, int n, int k, double q, NumericMatrix::Row Z) {
  double U = R::runif(0,1);
  int def_m = i + n;
  if(U < q) {
    int def = R::rbinom(def_m, p); // inherited by daughter
    def_m -= def;
    // daughter cell
    if(def <= k){
      Z[def]++;
    }
  }
  // mother cell 
  if(def_m <= k) {
    Z[def_m]++;
  }
}

// finds Z_bar over n generations, returns matrix
//[[Rcpp::export]] 
NumericMatrix Z_mat(double p, int n, int k, double q, int hours, NumericVector start) {
  NumericMatrix Z_m (hours+1, k+1);
  Z_m(0, _) = start;
  for(int gen = 1; gen <= hours; gen++) {
    for(int i = 0; i <= k; i++) {
      for(int j = 0; j < Z_m(gen-1, i); j++) {
        Z_def(i, p, n, k, q, Z_m(gen,_));
      }
    }
  }
  return(Z_m);
}

// Simulation of extinction probability 
//[[Rcpp::export]]
int Z_ext(double p, int n, int k, double q, int hours, NumericVector start) {
  NumericMatrix Z_m (hours+1, k+1);
  Z_m(0, _) = start;
  for(int gen = 1; gen <= hours; gen++) {
    for(int i = 0; i <= k; i++) {
      for(int j = 0; j < Z_m(gen-1, i); j++) {
        Z_def(i, p, n, k, q, Z_m(gen,_));
        // if sz > 5000 we are assuming that the popualtion survives
        if(sum(Z_m(gen,_) >= 5000)) {
          return(sum(Z_m(gen,_)));
        }
      }
    }
  }
  return(sum(Z_m(hours,_)));
}

// simulation of mean matrix 
//[[Rcpp::export]]
NumericMatrix M_est_cpp(double p, int n, int k, double q) {
  NumericMatrix M (k+1, k+1);
  for(int i = 0; i <= k; i++) {
    Z_def(i, p, n, k, q, M(i,_));
  }
  return(M);
}

// ----------------- INDIVIDUALS ----------------------

// Follows a mothercell until it dies
//[[Rcpp::export]]
int cell_age(int i, double p, int n, int k, double q){
  cell mth = {i, 0};
  while(mth.type <= k) {
    int def_m = mth.type+n;
    double U = R::runif(0,1);
    if(U < q) {
      int def = R::rbinom(def_m, p);
      def_m -= def;
    }
    if(def_m <= k) {
      mth.type = def_m;
      mth.age++;
      //bio_age.push_back(def_m);
    }
    else{
      return mth.age;
      //return bio_age;
    }
  }
}

// Returns numeric vector with all age counts
//[[Rcpp::export]]
NumericVector cell_ages(int i, double p, int n, int k, double q, int trials) {
  int max_age = 0;
  std::unordered_map<int, int> ages; // age and count
  for(int ii = 1; ii <= trials; ii++) {
    int age = cell_age(i, p, n, k,q);
    // update longest lifelength observed
    max_age = (age > max_age) ? age : max_age; 
    ages[age]++;
  }
  NumericVector age_vec (max_age+1);
  for(auto const &it: ages) {
    age_vec[it.first] = it.second;
  }
  return(age_vec);
}