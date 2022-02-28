#include <Rcpp.h> 
#include <list>
#include <iostream>
#include <iterator>
#include <unordered_map>
using namespace Rcpp;

struct cell {
  int type; 
  int age; // replicative age
};

// Checks lifelength of a cell starting from type i
int cell_age(int i, double p, int n, int k){
  cell mth = {i, 0};
  // change to avoid inf loop?
  int iter = 0;
  while(mth.type <= k) {
    iter++;
    int def_m = mth.type+n;
    int def = R::rbinom(def_m, p);
    def_m -= def;
    if(def_m <= k) {
      mth.type = def_m;
      mth.age++;
    }
    else{
      return mth.age;
    }
  }
}

// (TEST) Lifelengts of mother cells 
//[[Rcpp::export]]
NumericVector cell_ages(int i, double p, int n, int k, int iter) {
  int max_age = 0;
  std::unordered_map<int, int> ages; // age and count
  for(int ii = 1; ii <= iter; ii++) {
    int age = cell_age(i, p, n, k);
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


// Simulation core
void Z_def(int i, double p, int n, int k, NumericMatrix::Row Z){
  // do something
  int def_m = i + n;
  int def = R::rbinom(def_m, p);
  def_m -= def;
  // daughter cell
  if(def <= k){
    Z[def]++;
  }
  // mother cell
  if(def_m <= k) {
    Z[def_m]++;
  }
}

// Save simulation in NumericMatrix to then visualize
//[[Rcpp::export]]
NumericMatrix Z_mat(double p, int n, int k, int hours, NumericVector start) {
  NumericMatrix Z_m (hours+1, k+1);
  Z_m(0, _) = start;
  for(int gen = 1; gen <= hours; gen++) {
    for(int i = 0; i <= k; i++) {
      for(int j = 0; j < Z_m(gen-1, i); j++) {
        Z_def(i, p, n, k, Z_m(gen,_));
      }
    }
  }
  return(Z_m);
}

// Simulation of reproduction matrix
//[[Rcpp::export]]
NumericMatrix M_est(double p, int n, int k) {
  NumericMatrix M (k+1, k+1);
  for(int i = 0; i <= k; i++) {
    Z_def(i, p, n, k, M(i,_));
  }
  return(M);
}

// Simulation of extinction probability (only started)
//[[Rcpp::export]]
int Z_ext(double p, int n, int k, int hours, NumericVector start) {
  NumericMatrix Z_m (hours+1, k+1);
  Z_m(0, _) = start;
  for(int gen = 1; gen <= hours; gen++) {
    for(int i = 0; i <= k; i++) {
      for(int j = 0; j < Z_m(gen-1, i); j++) {
        Z_def(i, p, n, k, Z_m(gen,_));
        if(sum(Z_m(gen,_) >= 10000)) {
          return(sum(Z_m(gen,_)));
        }
      }
    }
  }
  return(sum(Z_m(hours,_)));
}