
# simple gw ---------------------------------------------------------------

vals <- c(0,1,2) # possible nr off offspring
p <- c(0.3, 0.3, 0.4) # offspring probability
n <- 10 # number of generations

# expected value
mu <- sum(vals*p)
expected <- mu^n

# simulation
simple_gw <- function(n) {
  # Z_0 = 1
  Z <- c(1, rep(0,n))
  for (i in 2:(n+1)) {
      if (Z[i-1] > 0) {
        Z[i] <- sum(sample(vals, size = Z[i-1], replace = T, prob = p))
      }
  }
  return(Z)
}

# number of trials
tr <- 10000
sim_mu <- mean(replicate(tr, simple_gw(n)[n+1]))

print("Expected value: ") 
print(expected)
print("Simulated value: ") 
print(sim_mu)

# multi type gw + reproduction matrix -------------------------------------

library(expm)

# estimate M trials time and take avarage
M_est <- function(p, n, k, q, hours, trials) {
  M <- Reduce('+', lapply(1:trials, function(i) M_mat(p, n, k, q)))/trials
  return(M %^% hours)
}

# build matrix from independent row
M_mat <- function(p, n, k, q) {
  M <- matrix(0, k+1, k+1)
  for(i in 0:k) {
    M[i+1,] <- M_row(i, p, n, k, q)
  }
  return(M)
}

M_row <- function(i, p, n, k, q) {
  M_i <- c(rep(0, k+1))

  U <- runif(1,0,1)
  def_m <- i+n # the mother cells ackumulates n new protein each hour
  
  # the cell divides
  if(U < q) {
    # daughter inherits defect proteins with pr = p
    def <- rmultinom(1, size = i+n, prob = c(1-p, p))
    def_m <- def[1] # defect proteins to mother cell
    if(def[2] <= k) {
      M_i[def[2]+1] <- 1 # survives
    }
  }
  
  # test mother cell 
  if(def_m <= k) {
    M_i[def_m+1] <- 1 # survives
  }
  
  return(M_i)
}

multi_gw <- function(p, n, k, q, hours, Z_0){
  Z <- Z_0
  for(i in 1:hours){
    # extinction Z_t is zero vector
    if(all(Z==0)) break
    # for all types that exist in this timestep run M_row
    Z <- rowSums(sapply(which(Z != 0), function(k_) { 
      rowSums(replicate(Z[k_], M_row(k_-1, p, n, k, q)))
    }))
  }
  return(Z)
}

# repeat multi type gw x nr of times and take mean
multi_sim <- function(p, n, k, q, hours, Z_0, trials) {
  # multi gw returns col vectors
  Z <- rowSums(replicate(trials, multi_gw(p, n, k, q, hours, Z_0)))/trials
  return(Z)
}




























