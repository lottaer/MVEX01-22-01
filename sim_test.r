
# simple gw ---------------------------------------------------------------
library(ggplot2)
library(tidyverse)

vals <- c(0,1,2) # possible nr off offspring
p <- c(0.3, 0.3, 0.4) # offspring probability
n <- 100 # number of generations

mu <- sum(vals*p)
expected <- mu^n

# simulation
simple_gw <- function(n, p) {
  # Z_0 = 1
  Z <- c(1, rep(0,n))
  for (i in 2:(n+1)) {
      if (Z[i-1] == 0) break
      Z[i] <- sum(sample(vals, size = Z[i-1], replace = T, prob = p))
  }
  return(Z[n+1])
}


# avarage size of population after n generations
sim_mu <- function(p, gen, trials) {
  mean(replicate(trials, simple_gw(n,p)))
}

# probability of extinction
ext_pr <- function(p, ext_time,  trials) {
  sum(replicate(trials, simple_gw(n, p)) == 0)/trials
}

# probability generating function
pgf <- function(s,p) {
  vals <- c(0:(length(p)-1))
  return(sum(s^vals*p))
}

# returns smallest positive (real) root to g(s)=s
pgf_root <- function(p)  {
  p[2] <- p[2]-1
  roots <- polyroot(p)
  real_root <- Re(roots)[abs(Im(roots)) < 1e-6]
  return(real_root[which.min(real_root > 0)])
}

### iterated function 
pgf_recurse <- function(x, p, n) {
  e <- pgf(x, p)
  i <- 1
  while(i < n) {
    i <- i+1
    e <- pgf(e, p)
  }
  return(e)
}
# TODO: convergence?

###########################################################################
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
                                            
# plot test ---------------------------------------------------------------

# test plot multi-type galton watson
# save output into txt file

sink(file = "multi-type_out.txt")
multi_gw(0.5, 2, 3, 0.2 ,15,c(1,0,0,0))
sink(file = NULL)

data <- read.table("multi-type_out.txt")
df <- data.frame(data[,2:ncol(data)])
df$x <- 0:(nrow(df)-1)

data_ggp <- data.frame(x = df$x,
                       y = c(df$V2, df$V3, df$V4, df$V5),
                       group = c(rep("type 0", nrow(data)),
                                 rep("type 1", nrow(data)),
                                 rep("type 2", nrow(data)),
                                 rep("type 3", nrow(data))))

ggp <- ggplot(data_ggp, aes(x, y, col = group)) +     
  geom_line(size = 1)
ggp  




























