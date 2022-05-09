#include <Rcpp.h> 
#include <list>
#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <array>
#include <deque>
#include <unordered_set>

using namespace Rcpp;

// For age distribution
struct cell {
  int type;   // keep track of biological age over time
  int age;    // to calculate age at death
};

struct cellp {
  int type;   // biological age
  int age;    // chronological age
  bool alive; 
  cellp *mth; // remember mother cell
};

// ----------------------- SIMULATE INDIVIDUAL CELLS -------------------------

// Follows a mothercell until it dies and returns its chronological lifelength
//[[Rcpp::export]]
int cell_age(int i, double p, int n, int k, double q){
  cell mth = {i, 0};
  while(mth.type <= k) {
    mth.age++;
    int def_m = mth.type+n;
    double U = R::runif(0,1);
    if(U < q) {
      int def = R::rbinom(def_m, p);
      def_m -= def;
    }
    if(def_m <= k) {
      mth.type = def_m;
    }
    else{
      return mth.age;
    }
  }
  return(mth.age);
}

// Follows a mothercell until it dies and returns its replicative lifelength
//[[Rcpp::export]]
int repcell_age(int i, double p, int n, int k, double q){
  cell mth = {i, 0};
  while(mth.type <= k) {
    int def_m = mth.type+n;
    double U = R::runif(0,1);
    if(U < q) {
      mth.age++;
      int def = R::rbinom(def_m, p);
      def_m -= def;
    }
    if(def_m <= k) {
      mth.type = def_m;
    }
    else{
      return mth.age;
    }
  }
  return(mth.age);
}

//---------------------REJUVENATION INDEX-------------------------

// Simulate trials mothercells and simulate its daughter cells
// Return the difference in their lifelength
//[[Rcpp::export]]
std::list<double> df_drls(int i, double p, int n, int k, int trials){
  std::list<double> drls;
  for (int ii = 0; ii < trials; ii++){
    // ONe trial
    std::deque<int> dage;
    cell mth = {i, 0}; // replicative age
    while(mth.type <= k) {
      mth.age++;
      int def_m = mth.type+n;
      int def = R::rbinom(def_m, p);
      def_m -= def;
      // Add only if daughter is alive
      if (def <= k){
        dage.push_back(def);
      }
      if(def_m <= k) {
        mth.type = def_m;
      }
      else{
        //Simulate all daughter ages and return difference in age to mother
        double mage = mth.age;
        for (int i = 0; i < dage.size(); i++) {
          double x = (repcell_age(dage[i], p, n, k, 1)-mage);
          drls.push_back((double) x); // x/mage to divide by mothers age
          //drls.push_back(cell_age(dage[i], p, n, k)-mage)
        }
        break;
        
      }  
    }
  }
  return drls;
}

// ----------------- SIMULATE DISTRIBUTION ----------------------

// Returns vector of all age counts from simulation
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

// ------------------- SIMULATION OF MULTI GW  --------------------------

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

// Simulate the population over x number of hours 
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
        // if sz > 1000 we are assuming that the popualtion survives
        if(sum(Z_m(gen,_) >= 1000)) {
          return(sum(Z_m(gen,_)));
        }
      }
    }
  }
  return(sum(Z_m(hours,_)));
}

// simulation of mean value matrix
//[[Rcpp::export]]
NumericMatrix M_est_cpp(double p, int n, int k, double q) {
  NumericMatrix M (k+1, k+1);
  for(int i = 0; i <= k; i++) {
    Z_def(i, p, n, k, q, M(i,_));
  }
  return(M);
}

// Finding lifespan distribution of type zero cells 
//[[Rcpp::export]]
std::list<int> zeroage(int i, double p, int n, int k, double q, int hours) {
  std::list<cell> Z; // alive cells 
  Z.push_front({i, 0});
  for(int j = 1; j <= hours; j++) {
    for(auto it = Z.begin(); it != Z.end(); ) {
      double U = R::runif(0,1);
      int def_m = it->type + n;
      // the cell divides
      if(U < q) {
        int def = R::rbinom(def_m, p);
        def_m -= def;
        if(def <= k) {
          Z.push_front({def, 0});
        }
      }
      it->age++; // survived
      if(def_m <= k){
        it->type = def_m;
        ++it;
      }
      else {
        it = Z.erase(it); // remove from alive cells
      }
    }
  }
  std::list<int> zerotype;
  for(auto it = Z.begin(); it != Z.end(); ) {
    Rcout << "Im here";
    if (it->type == 0){
      zerotype.push_back(it->age);
    }
    ++it;
  }
  return zerotype;
}
