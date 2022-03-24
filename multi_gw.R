# multi type gw + reproduction matrix -------------------------------------

library(expm)
library(ggplot2)
library(tidyverse)
library(pracma)
library(patchwork)
library(MASS)
library(easyGgplot2)

### Section: Simulation

# simulation of the reproduction matrix, core in cpp file
M_est <- function(p, n, k, q, hours, trials) {
  M <- Reduce('+', lapply(1:trials, function(i) M_est_cpp(p, n, k, q) %^% hours))/trials
  return(M)
}

# simulation of multi gw, takes mean of trials simulations
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

### Section: plots for MGW

# convert to right format (long instead of wide)
to_long <- function(Z) {
  df = as.data.frame(Z)
  colnames(df) <- 0:(ncol(df)-1)
  df$x <- 0:(nrow(df)-1)
  df_long <- df %>% gather(key = "Type", value = "Size", -x)
  return(df_long)
}

# Calculated proportions of types over time (stable type distribution)
type_frequency <- function(Z) {
  row_sums <- rowSums(Z)
  prop <- t(mapply(function(r) Z[r,]/row_sums[r], 1:nrow(Z)))
  df <- as.data.frame(prop)
  colnames(df) <- 0:(ncol(df)-1)
  return(df)
}

# plot using ggplot2
multi_gw_pl <- function(p, n, k, q, hours, Z_0, trials) {
  
  theme_set(theme_minimal())
  
  df_wide <- multi_gw_mean(p, n, k, q, hours, Z_0, trials)
  df_long <- to_long(df_wide)
  
  M_analytic <- type_mat(p, n, k, q,1)
  stable_dist <- round(pf_eigen(M_analytic)[[2]], digits = 3)
  print(stable_dist)
  
  type_long <- to_long(type_frequency(df_wide))
  pl1 <- ggplot(df_long, aes(x, Size)) + geom_line(aes(color = Type, group = Type), size = 1.2)
  pl2 <- ggplot(type_long, aes(x, Size)) + 
    geom_line(aes(color = Type, group = Type), size = 1.2) + 
    labs(y = "Proportion") + geom_hline(yintercept = stable_dist, color = 'grey',size=0.8) +
    scale_y_continuous(breaks = stable_dist)

  pl1 + pl2
  
}

### Section: Numerical computations

# calculated the reproduction matrix M = qA+B as in theory
type_mat <- function(p,n,k,q,hours) {
  # M = qA + B
  M <- matrix(0, k+1, k+1)
  for(i in 0:k) {
    M[i+1,] <- q*type_row(i, p, n, k) + b_row(i,n,k,q)
  }
  return(M %^% hours)
}

# rows in the A matrix
type_row <- function(i, p, n, k) {
  sapply(0:k, function(j) dbinom(i+n-j, n+i, p) + dbinom(j, n+i, p))
}
# rows in the B matrix
b_row <- function(i, n, k, q) {
  b <- c(rep(0,k+1))
  if(i+n <= k) {
    b[i+n+1] <- 1-q
  }
  return(b)
}

# perron frobenius eigenvalue and eigenvector
pf_eigen <- function(M) {
  eigs = eigen(t(M))
  vals <- eigs$values
  vecs <- eigs$vectors
  ind <- which.max(Re(vals[abs(Im(vals)) < 1e-6]))
  return(list(vals[ind], c(normalize(Re(vecs[,ind])))))
}

# normalize st u*1=1
normalize <- function(vec) {
  return(vec/sum(vec))
}
         
# expected size of population
expected_sz <- function(p, n, k, q, hours, start) {
  M <- type_mat(p,n,k,q,hours)
  sum(start*M)
}

# minimization of rho - 1
critical <- function(p, n, k) {
  qs <- function(q) { abs(pf_eigen(type_mat(p, n, k, q,1))[[1]]-1)}
  optimize(qs, lower = 0, upper = 1)$minimum
}

# results to data frame to be able to plot
critical_df <- function(n,k) {
  p <- seq(0,1,by=0.01)
  qcrit <- sapply(p, function(p) critical(p,n,k))
  df <- as.data.frame(qcrit)
  df$p <- p
  return(df)
}

# plots the critical q value as a function of p
critical_pl <- function(n,k) {
  df <- critical_df(n,k)
  ggplot(df, aes(p, qcrit)) +
    geom_line(size = 1.2) +
    ggtitle("Kritiskt q värde") + 
    xlab("p") + ylab("q") 
}

# plots several curves in the same graph
critical_pls <- function(n,k) {
  
  theme_set(theme_minimal())
  
  df <- critical_df(n, k[1])
  df$k <- k[1]
  for(i in 2:length(k)) {
    df2 <- critical_df(n,k[i])
    df2$k <- k[i]
    df <- rbind(df, df2)
  }
  ggplot(df, aes(p, qcrit, color = factor(k))) +
    geom_line(size = 1.2) +
    ggtitle("Kritiskt q värde") + 
    xlab("p") + ylab("q") +
    coord_fixed(xlim = c(0,1), ylim = c(0,1))
}

## Section: Age dependent process

# find distribution of ages among mother cells
age_sim <- function(p, n, k, q, hours, Z_0, trials) {
  rowSums(replicate(trials, age_prop(p, n, k, q, hours, Z_0)))
}

# convert to distribution (scales)
age_df <- function(Z) {
  scaled <- Z/sum(Z);
  df <- as.data.frame(scaled)
  df$x <- 0:hours
  return(df)
}

# WORK IN PROGRESS: fit gamma disribution to simulated data
# fit gamma disribution to simulated data by maximum likelihoos
gamma_ML <- function(Z) {
  data <- fitdistr(rep(c(0.01, 1:(length(Z)-1)), times=Z), "gamma")
  print(data)
  params <- data$estimate
  df <- as.data.frame(Z/sum(Z))
  colnames(df) <- 'proportion'
  df$x <- 0:(length(Z)-1)
  ggplot(df, aes(x = x, y = proportion)) + geom_bar(stat = 'identity') +
    geom_function(fun = dgamma, args = list(shape = params[1], rate = params[2]), size = 1) 
}
                  
# Frequency polygons
age_pl <- function(p, n, k, q, hours, Z_0, trials) {
  theme_set(theme_minimal())
  
  # run the simulation
  Z <- sapply(p, function(p_) rowSums(replicate(trials, age_prop(p_, n, k, q, hours, Z_0))))
  df <- as.data.frame(sapply(1:ncol(Z), function(i) Z[,i]/sum(Z[,i])))
  colnames(df) <- p
  df$x <- 0:hours
  df_long <- df %>% gather(key = "pvalue", value = "count", -x)
  
  ggplot(df_long, aes(x, count)) + geom_line(aes(color = pvalue, group = pvalue), size = 1) 
}


