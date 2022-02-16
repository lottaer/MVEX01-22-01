# multi type gw + reproduction matrix -------------------------------------

library(expm)
library(ggplot2)
library(tidyverse)

# estimate M trials time and take average
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
      M_i[def[2]+1] <- M_i[def[2]+1] + 1 # survives
    }
  }
  
  # test mother cell 
  if(def_m <= k) {
    M_i[def_m+1] <- M_i[def_m+1] + 1 # survives
  } 
  
  return(M_i)
}

multi_gw <- function(p, n, k, q, hours, Z_0){
  Z <- Z_0
  print(Z_0)
  for(i in 1:hours){
    # extinction Z_t is zero vector
    if(all(Z==0)) break
    # for all types that exist in this timestep run M_row
    Z <- rowSums(sapply(which(Z != 0), function(k_) { 
      rowSums(replicate(Z[k_], M_row(k_-1, p, n, k, q)))
    }))
    print(Z)
  }
  # return(Z)
} 

# repeat multi type gw x nr of times and take mean
multi_sim <- function(p, n, k, q, hours, Z_0, trials) {
  # multi gw returns col vectors
  Z <- rowSums(replicate(trials, multi_gw(p, n, k, q, hours, Z_0)))/trials
  return(Z)
}


# plot test ---------------------------------------------------------------

library(patchwork)

# test plot multi-type galton watson
# multi gw dataframe

multi_gw_df <- function(p, n, k, q, hours, Z_0) {
  Z_mat <- matrix(0,hours+1, k+1)
  Z_mat[1,] <- Z_0
  for(i in 2:(hours+1)){
    # extinction Z_t is zero vector
    if(all(Z_mat==0)) break
    # for all types that exist in this timestep run M_row
    Z_mat[i,] <- rowSums(sapply(which(Z_mat[i-1,] != 0), function(k_) { 
      rowSums(replicate(Z_mat[i-1,k_], M_row(k_-1, p, n, k, q)))
    }))
  }
  return(Z_mat)
}

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
multi_gw_pl <- function(p, n, k, q, hours, Z_0) {
  
  theme_set(theme_minimal())
  
  df_wide <- multi_gw_df(p, n, k, q, hours, Z_0)
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

######## Find critical value for (p,q) (function of p)

type_mat <- function(p,n,k,q,hours) {
  # M = qA + B
  M <- matrix(0, k+1, k+1)
  for(i in 0:k) {
    M[i+1,] <- q*type_row(i, p, n, k) + b_row(i,n,k,q)
  }
  return(M %^% hours)
}

type_row <- function(i, p, n, k) {
  sapply(0:k, function(j) dbinom(i+n-j, n+i, p) + dbinom(j, n+i, p))
}

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

normalize <- function(vec) {
  return(vec/sum(vec))
}


library("pracma")
critical <- function(p, n, k) {
  qs <- function(q) { abs(pf_eigen(type_mat(p, n, k, q,1))[[1]]-1)}
  optimize(qs, lower = 0, upper = 1)$minimum
}

critical_df <- function(n,k) {
  p <- seq(0,1,by=0.01)
  qcrit <- sapply(p, function(p) critical(p,n,k))
  df <- as.data.frame(qcrit)
  df$p <- p
  return(df)
}

critical_pl <- function(n,k) {
  df <- critical_df(n,k)
  ggplot(df, aes(p, qcrit)) +
    geom_line(size = 1.2) +
    ggtitle("Kritiskt q värde") + 
    xlab("p") + ylab("q") 
}

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

## Gamma

# inte färdig 
age_sim <- function(p, n, k, q, hours, Z_0, trials) {
  
  theme_set(theme_minimal())
  
  # multi gw returns col vectors
  Z <- rowSums(replicate(trials, age_prop(p, n, k, q, hours, Z_0)))
  params <- fitdistr(rep(0:hours, times=Z), "gamma")$estimate
  scaled <- Z/sum(Z);
  df <- as.data.frame(scaled)
  df$x <- 0:hours
  
  #curve(dgamma(x, shape = params[1], rate = params[2]), from = 0, to = 10)
  ggplot(data = df, aes(x = x, y = scaled)) + geom_bar(stat="identity") +
    geom_function(fun = dgamma, args = list(shape = params[1], rate = params[2]))
}












