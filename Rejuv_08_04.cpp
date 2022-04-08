#include <Rcpp.h> 
#include <list>
#include <iostream>
#include <iterator>
#include <unordered_map>
#include <array>
#include <vector>
#include <deque>
using namespace Rcpp;

// defines a cell, save type and age for every object
struct cell {
  int type; 
  int age; // replicative age
};

// TODO
// Follows a mothercell until it dies
//[[Rcpp::export]]
int cell_age(int i, double p, int n, int k){
  cell mth = {i, 0};
  //std::list<int> bio_age;
  //bio_age.push_back(i);
  // change to avoid inf loop?
  while(mth.type <= k) {
    int def_m = mth.type+n;
    int def = R::rbinom(def_m, p);
    def_m -= def;
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

//----------------------------Rejuvination-Index using  deque-----------------
// Functions for calculations of rejuvination index
//Saves daughter types through a mothers life, simulates their age
// and returns the difference in ages between mother and each daugther
// Just nu kan man ej ändra q
//[[Rcpp::export]]
NumericVector deqrls(int i, double p, int n, int k){
  std::deque<int> dage;
  NumericVector drls;
  cell mth = {i, 0};
  while(mth.type <= k) {
    int def_m = mth.type+n;
    int def = R::rbinom(def_m, p);
    def_m -= def;
    dage.push_back(def);
    if(def_m <= k) {
      mth.type = def_m;
      mth.age++;
    }
    else{
      //Simulate all daughter ages and return difference in age to mother
      int mage = dage.size();
      for (int i = 0; i < dage.size(); i++) {
        drls.push_back(cell_age(dage[i], p, n, k)-mage);
      }
      return drls;
    }  
  }
}


//---------------------------------------------------------------

// This function is applied to all the cells
// Let the cell divide and look at which types it becomes and verify that they survive
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

// Save all timesteps into a matrix
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
