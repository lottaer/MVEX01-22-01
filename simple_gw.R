
# simple gw ---------------------------------------------------------------

library(ggplot2)
library(tidyverse)

# PARAMS:
# n = number of generations
# p = probability vector (offspring distribution)
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
ext_pr_grid <- function(sz, n, trials) {
  pr <- seq(0,0.5, length.out = sz)
  gr <- expand.grid(pr, pr) # creates a sz*sz grid och possible p1, p2 values
  p0 <- 1-rowSums(gr)  # p0 = 1-p1-p2
  # remove all unfeasible combinations
  data <- data.frame(p0, gr) %>% filter_all(all_vars(.>=0 & .<=1))
  colnames(data) <- c('p0', 'p1', 'p2')
  # simulated population extinction for all possible points
  data$Q <- sapply(1:nrow(data), function(row) ext_pr(data[row,], n, trials))
  return(data)
}

library("RColorBrewer")
gradient_plot <- function(data) {
  eq <- function(x) {0.5*(1-x)}
  pl <- ggplot(data, aes(x = p1, y = p2, color = Q)) + stat_function(fun=eq, size = 1) +  geom_point() + scale_color_gradientn(colours = c("white", "yellow", "orange", "darkred"))
  pl
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

# visualization of convergence 
# TODO: plot multiple lines
pgf_plot <- function(x, p, n) {
  # extinction probability
  ext <- pgf_root(p)
  ll <- sapply(1:n, function(tt) pgf_recurse(x, p, tt))
  data <- data.frame(iterations = c(1:n), probability = ll)
  pl <- ggplot(data, aes(x=iterations, y=probability)) + geom_line(lwd=1.2) 
  pl + geom_hline(yintercept = ext, linetype = "dashed", color = "red", lwd = 1.2)
}
