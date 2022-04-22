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
          Z.push_front({ def, 0});
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
    if (it->type == 0){
      zerotype.push_back(it->age);
    }
    ++it;
  }
  return zerotype;
}
