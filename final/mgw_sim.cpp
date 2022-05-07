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

// Follows a mothercell until it dies
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


// Follows a mothercell until it dies replicative age (will be same as normal while q=1)
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
    mth.age++;
    // Add only if daughter is alive
    if (def <= k){
      dage.push_back(def);
    }
    if(def_m <= k) {
      mth.type = def_m;
    }
    else{
      //Simulate all daughter ages and return difference in age to mother
      int mage = mth.age;
      for (int i = 0; i < dage.size(); i++) {
        // Return normalized values
        //double x = (cell_age(dage[i], p, n, k));
        drls.push_back((double) (cell_age(dage[i], p, n, 1, k)-mage)/mage);
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
          drls.push_back((double) x); //Dela med moder ålder
          //drls.push_back(cell_age(dage[i], p, n, k)-mage)
        }
        break;
        
      }  
    }
  }
  return drls;
}

//---------------------------------------------------------------

std::list<cell> get_daughters(cell &mth, double p, int n, int k, double q) {
  std::list<cell> dht; // daughter cells
  while(mth.type <= k) {
    int def;
    int def_m = mth.type+n;
    double U = R::runif(0,1);
    if(U < q) {
      def = R::rbinom(def_m, p);
      def_m -= def;
      // The daughter cell survives
      if(def <= k){
        dht.push_back({def, 1});
      }
    }
    // check if the mother survives
    if(def_m <= k) {
      mth.type = def_m;
      mth.age++;
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

// Dela med moders livslängd eller mean av ALLA livslänger?
//[[Rcpp::export]]
std::vector<double> multi_rej(int i, double p, int n, int k, double q, int trials){
  int tot_age = 0;
  int non_empty = 0;
  std::vector<double> rej;
  for(int tr = 0; tr < trials; tr++) {
    cell mth = {i, 0}; // init mother cell
    std::list<cell> dht = get_daughters(mth,p,n,k,q);
    int m_age = mth.age;
    int d_age = 0;
    // Lifelength of daughters
    for(cell daughter : dht) {
      int rep_age = cell_age(daughter.type,p,n,k,q);
      rej.push_back(rep_age-m_age);
    }
  }
  return(rej);
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

// ----------------- INDIVIDUALS ----------------------

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
        // if sz > 1000 we are assuming that the popualtion survives
        if(sum(Z_m(gen,_) >= 1000)) {
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
