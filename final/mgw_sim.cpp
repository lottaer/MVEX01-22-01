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

struct cellp {
  int type;
  int age;
  bool alive;
  cellp *mth;
};

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
    }
    else{
      return mth.age;
    }
  }
  return(mth.age);
}

std::list<cell> get_daughters(cell &mth, double p, int n, int k, double q) {
  std::list<cell> dht; // daughter cells
  while(mth.type <= k) {
    int def;
    int def_m = mth.type+n;
    double U = R::runif(0,1);
    if(U < q) {
      def = R::rbinom(def_m, p);
      def_m -= def;
    }
    if(def_m <= k) {
      mth.type = def_m;
      mth.age++;
      dht.push_back({def, 0});
    }
    else{
      return(dht);
    }
  }
  return(dht);
}

// calculate the vector mean
//[[Rcpp::export]]
double rcpp_mean(std::vector<int> cells) {
  int sum = 0;
  int n = cells.size();
  for(auto const& value: cells) {
    sum += value;
  }
  return((double) sum/n);
}

// divide all elements in vector by scalar
//[[Rcpp::export]]
std::vector<double> divide(double mth_mean, std::vector<double> diff) {
  for(int i = 0; i < diff.size(); i++) {
    diff[i] = (double) diff[i]/mth_mean;
  }
  return(diff);
}

//[[Rcpp::export]]
std::vector<double> multi_rej(int i, double p, int n, int k, double q, int trials){
  
  int tot_age = 0;
  int non_empty = 0;
  std::vector<double> index;
  
  for(int tr = 1; tr <= trials; tr++) {
    
    cell mth = {i,0};    // init mother cell
    std::list<cell> dht = get_daughters(mth,p,n,k,q);
    std::vector<int> rej;
    
    // Lifelength of daughters
    for(cell daughter : dht) {
      // if no daughters skip (both needs to die)
      int rep_age = cell_age(daughter.type,p,n,k,q);
      rej.push_back(rep_age);
    }
    
    if(!rej.size()==0) {
      tot_age += mth.age;
      non_empty++;
      double diff = rcpp_mean(rej)-mth.age;
      index.push_back(diff);
    }
  }
  double mean = (double) tot_age/non_empty;
  return(divide(mean, index));
}

// calculates the difference between lifelength of mth and dth
// only makes sense in the special case
std::list<double> rej_index(std::list<cellp> dead) {
  std::list<double> rej;
  for(auto it = dead.begin(); it != dead.end(); ) {
    if(!(it->mth == NULL)) {
      if(!it->mth->alive) {
        double re = it->age - it->mth->age;
        it = dead.erase(it);
        rej.push_front(re);
      }
    }
    ++it;
  }
  return(rej);
}

//[[Rcpp::export]]
std::list<double> rej_cell(int i, double p, int n, int k, double q, int hours) {
  std::list<cellp> Z; // alive cells 
  Z.push_front({0, 0, true, NULL});
  
  std::list<cellp> dead; // add dead cells
  
  for(int j = 1; j <= hours; j++) {
    for(auto it = Z.begin(); it != Z.end(); ) {
      double U = R::runif(0,1);
      int def_m = it->type + n;
      // the cell divides
      if(U < q) {
        int def = R::rbinom(def_m, p);
        def_m -= def;
        if(def <= k) {
          Z.push_front({ def, j, true, &(*it) });
        }
        // The cell dies immediately
        else {
          dead.push_front({ def, 0, false, &(*it) });
        }
      }
      if(def_m <= k){
        it->type = def_m;
        it->age++; // survived
        ++it;
      }
      else {
        it->alive = false;
        dead.push_front(*it); // add to dead cells
        it = Z.erase(it); // remove from alive cells
      }
    }
  }
  return(rej_index(dead)); // return the dead cells
}

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
