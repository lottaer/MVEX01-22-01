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
  cell mth = {i, 1};
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
std::vector<double> deqrls(int i, double p, int n, int k){
  std::deque<int> dage;
  std::vector<double> drls;
  cell mth = {i, 0};
  while(mth.type <= k) {
    int def_m = mth.type+n;
    int def = R::rbinom(def_m, p);
    def_m -= def;
    // Add only if daughter is alive
    if (def <= k){
    dage.push_back(def);
    }
    if(def_m <= k) {
      mth.type = def_m;
      mth.age++;
    }
    else{
      //Simulate all daughter ages and return difference in age to mother
      int mage = dage.size();
      for (int i = 0; i < dage.size(); i++) {
        // Return normalized values
        //double x = (cell_age(dage[i], p, n, k));
        drls.push_back((double) (cell_age(dage[i], p, n, k)-mage)/mage);
      }
      return drls;
    }  
  }
}


//---------------------Repeat add all to same list-------------------------
// NEW version
//[[Rcpp::export]]
std::list<double> df_drls(int i, double p, int n, int k, int trials){
  std::list<double> drls;
  for (int ii = 0; ii < trials; ii++){
    // ONe trial
    std::deque<int> dage;
    cell mth = {i, 0};
    while(mth.type <= k) {
      int def_m = mth.type+n;
      int def = R::rbinom(def_m, p);
      def_m -= def;
      // Add only if daughter is alive
      if (def <= k){
        dage.push_back(def);
      }
      if(def_m <= k) {
        mth.type = def_m;
        mth.age++;
      }
      else{
        //Simulate all daughter ages and return difference in age to mother
        double mage = dage.size();
        for (int i = 0; i < dage.size(); i++) {
          double x = (cell_age(dage[i], p, n, k)-mage);
          drls.push_back((double) x/mage); //Dela med moder ålder
          //drls.push_back(cell_age(dage[i], p, n, k)-mage)
        }
        break;
        
      }  
    }
  }
  return drls;
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
