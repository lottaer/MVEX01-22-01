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

# proportion of daughters that lives longer than the mother
prop_rej <- function(data) {
  sum(data >= 0)/length(data)
}

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

### ---------------INDEX ---------------------------

scale_factor <- function(p,n,k,q,trials) {
  stable_dist <- st_type(p,n,k,q)
  EX <- rep(0, k+1)
  for(i in 0:k) {
    data <- cell_ages(i,p,n,k,q,trials)
    x <- 0:(length(data)-1)
    EX[i+1] <- mean(rep(x, times = data))
  }
  # weighted mean
  sum(stable_dist*EX)
}

### ---------------GEOM -----------------------------

# R-convention: The probability of the number Y = X - 1 of 
# failures before the first succes -> 0 included
MLE_geom <- function(data) {
  x <- 0:(length(data)-1)
  tmp <- rep(x, times = data)
  n <- length(tmp)
  n/(sum(tmp)+n) # maximum likelihood parameter
}

# The input data should be a vector from cell_ages
MLE_plot <- function(data) {
  plot(0:(length(data)-1), data/sum(data), type = 'l', ylab = '', xlab = 'Ålder')
  start <- readline(prompt = "Enter starting value: ")
  tail <- data[-(1:start)]
  mle <- MLE_geom(tail)
  print(mle)
  df <- data.frame(prop = tail/sum(tail), x = 0:(length(tail)-1))
  ggplot(df, aes(x, prop)) + geom_line(aes(x, prop), size = 1) +
    geom_line(aes(x, dgeom(x,mle)), color = "red") 
}

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
  M <- M_mat(p,n,k,q,1)
  pf_eigen(M)[[1]]
}

# stable type distribution
st_type <- function(p,n,k,q) {
  M <- M_mat(p,n,k,q,1)
  round(pf_eigen(M)[[2]],3)
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
