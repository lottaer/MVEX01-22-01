library(expm)
library(ggplot2)
library(tidyverse)
library(pracma)
library(patchwork)
library(MASS)
library(easyGgplot2)
library(nleqslv)
library(tikzDevice)
library(comprehenr)

# Main file 

### -------------- SIMULATION ------------------

# simulation of multi gw, takes mean of trials simulations and returns data frame
multi_gw_mean <- function(p, n, k, q, hours, Z_0, trials) {
  if(length(Z_0) != k+1) {
    print("Enter starting vector of correct length")
    return(0)
  }
  Reduce('+', lapply(1:trials, function(i) Z_mat(p, n, k, q, hours, Z_0)))/trials
}

# simulation of the extinction probability
multi_gw_ext <- function(p,n,k,q,hours,Z_0,trials) {
  if(length(Z_0) != k+1) {
    print("Enter starting vector of correct length")
    return(0)
  }
  sum(replicate(trials, Z_ext(p,n,k,q,hours,Z_0))==0)/trials
}

# simulation of the reproduction matrix, core in cpp file
M_est <- function(p, n, k, q, hours, trials) {
  M <- Reduce('+', lapply(1:trials, function(i) M_est_cpp(p, n, k, q) %^% hours))/trials
  return(M)
}
                          
index <- function(i, p, n , k, q, trials) {
  mean(multi_rej(i, p, n, k, q, trials))
}

# simulation of age distribution
age_sim <- function(p, n, k, q, hours, Z_0, trials) {
  rowSums(replicate(trials, age_prop(p, n, k, q, hours, Z_0)))
}

# scale to distribution
age_df <- function(Z) {
  scaled <- Z/sum(Z);
  df <- as.data.frame(scaled)
  df$x <- 0:hours
  return(df)
}

### -------------- NUMERICAL -----------------------
# All numerical computations used in plots etc.

## ---------------- M-MATRIX -----------------------

# Computation of mean value matrix: M = qA + B
M_mat <- function(p,n,k,q,hours) {
  M <- matrix(0, k+1, k+1)
  for(i in 0:k) {
    M[i+1,] <- q*A_row(i, p, n, k) + B_row(i,n,k,q)
  }
  return(M %^% hours)
}

# expected size of population given start value
expected_sz <- function(p, n, k, q, hours, start) {
  M <- M_mat(p,n,k,q,hours)
  sum(start*M)
}

## ---------------- PERRON-FROBENIUS ---------------
# perron frobenius eigenvalue and eigenvector of M
pf_eigen <- function(M) {
  eigs = eigen(t(M))
  vals <- eigs$values
  vecs <- eigs$vectors
  ind <- which.max(Re(vals[abs(Im(vals)) < 1e-6]))
  return(list(Re(vals[ind]), c(normalize(Re(vecs[,ind])))))
}

# reproduction value
rho <- function(p,n,k,q) {
  M <- type_mat(p,n,k,q,1)
  pf_eigen(M)[[1]]
}

# stable type distribution
st_type <- function(p,n,k,q) {
  M <- type_mat(p,n,k,q,1)
  pf_eigen(M)[[2]]
}

## ------------- HELPER-FUNCTIONS ------------------

# rows in the A matrix
A_row <- function(i, p, n, k) {
  sapply(0:k, function(j) dbinom(i+n-j, n+i, p) + dbinom(j, n+i, p))
}
# if needed add the B matrix
B_row <- function(i, n, k, q) {
  b <- c(rep(0,k+1))
  if(i+n <= k) {
    b[i+n+1] <- 1-q
  }
  return(b)
}


# normalize st u*1 = 1
normalize <- function(vec) {
  return(vec/sum(vec))
}
