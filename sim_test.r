# simple gw ---------------------------------------------------------------

library(ggplot2)
library(tidyverse)

# PARAMS:
# n = number of generations
# p = probability vector (offspring distribution)
# trials = number of trials in simulation
# sz = grid size

# simulation if simple galton watson process
simple_gw <- function(n, p) {
  # Z_0 = 1
  Z <- c(1, rep(0,n))
  for (i in 2:(n+1)) {
      if (Z[i-1] == 0) break
      Z[i] <- sum(sample(0:(length(p)-1), size = Z[i-1], replace = T, prob = p))
  } 
  return(Z[n+1])
}

# rate of reproduction (m-value)
mu <- function(p) {
  sum(sapply(1:length(p), function(i) (i-1)*p[i]))
}

# simulated m value
sim_mu <- function(p, n, trials) {
  mean(replicate(trials, simple_gw(n,p)))
}

# simulation of probability of extinction
ext_pr <- function(p, n,  trials) {
  sum(replicate(trials, simple_gw(n, p)) == 0)/trials
}

# returns smallest positive (real) root to g(s)=s
pgf_root <- function(p)  {
  p[2] <- p[2]-1
  roots <- polyroot(p)
  real_root <- Re(roots)[abs(Im(roots)) < 1e-6]
  return(real_root[which.min(real_root > 0)])
}

# for what values of p1 and p2 does extinction occur
extinct <- function(sz, n, trials) {
  pr <- seq(0,1, length.out = sz)
  # creates a sz*sz grid och possible p1, p2 values
  gr <- expand.grid(pr, pr)
  # p0 = 1-p1-p2
  p0 <- 1-rowSums(gr)
  # remove all unfeasible combinations
  data <- data.frame(p0, gr) %>% filter_all(all_vars(.>=0 & .<=1))
  
  # simulated population extinction for all possible points
  ext <- sapply(1:nrow(data), function(row) ext_pr(data[row,], n, trials))

  pl <- ggplot(data[,2:3], aes(data[,2], data[,3])) + geom_point(aes(color = ext))
  pl
  # sp3+scale_color_gradientn(colours = rainbow(10))
}
# TODO: graphical, add constraint

### iterated function -------------------------------------------------------

# probability generating function
pgf <- function(s,p) {
  vals <- c(0:(length(p)-1))
  return(sum(s^vals*p))
}

# iterated function
pgf_recurse <- function(x, p, n) {
  e <- pgf(x, p)
  i <- 1
  while(i < n) {
    i <- i+1
    e <- pgf(e, p)
  }
  return(e)
}
                
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
  # TODO: convert to correct data frame
  return(to_dataframe(Z_mat))
}

# convert to right format (long instead of wide)
to_dataframe <- function(Z) {
  df = as.data.frame(Z)
  colnames(df) <- 0:(ncol(df)-1)
  df$x <- 0:(nrow(df)-1)
  df_long <- df %>% gather(key = "Type", value = "Size", -x)
  return(df_long)
}

# plot using ggplot2
multi_gw_pl <- function(p, n, k, q, hours, Z_0) {
  df_long <- multi_gw_df(p, n, k, q, hours, Z_0)
  ggp <- ggplot(df_long, aes(x, Size)) + geom_line(aes(color = Type, group = Type), size = 1.2) + theme_light()
  ggp  + theme(
    aspect.ratio = 1,
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15, face = 'bold')
  )
  ggp
  # TODO: labels, save option
}




























