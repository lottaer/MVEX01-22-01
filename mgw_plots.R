library(expm)
library(ggplot2)
library(tidyverse)
library(pracma)
library(patchwork)
library(MASS)
library(easyGgplot2)
library(nleqslv)
library(tikzDevice)

### ------------------ GENERAL CASE -------------------

# plot stable type distribution
multi_gw_pl <- function(p, n, k, q, hours, Z_0, trials) {
  
  #theme_set(theme_minimal())
  
  df_wide <- multi_gw_mean(p, n, k, q, hours, Z_0, trials)
  df_long <- to_long(rowSums(df_wide))
  
  M_analytic <- type_mat(p, n, k, q,1)
  stable_dist <- round(pf_eigen(M_analytic)[[2]], digits = 3)
  print(stable_dist)
  
  type_long <- to_long(type_frequency(df_wide))
  #pl1 <- ggplot(df_long, aes(x, Size)) + geom_line(aes(color = Type, group = Type), size = 1.2) +
  # xlab("Generation") + ylab("Populationsstorlek")
  pl2 <- ggplot(type_long, aes(x, Size)) + 
    geom_line(aes(color = Type, group = Type), size = 1.2) + xlab("Generation") +
    labs(y = "Proportion") + geom_hline(yintercept = stable_dist, color = 'grey',size=0.8) +
    scale_y_continuous(breaks = stable_dist)
  
  #pl1  + pl2
  pl2
}

# stacked bar plots for stable type distribution
v_plot <- function(p,n,k,q) {
  # check supercriticality
  rhos <- sapply(p, function(p_) rho(p_,n,k,q))
  if(min(rhos) <= 1) {
    print("Subcritical!")
    return(0)
  }    
  v <- st_type(p[1],n,k,q)
  label <- rep(p[1], k+1)
  type <- 0:k
  df <- data.frame(label, v,type)
  for(i in 2:length(p)){
    v <- st_type(p[i],n,k,q)
    label <- rep(p[i], k+1)
    type <- 0:k
    df2 <- data.frame(label, v,type)
    df <- rbind(df, df2)
  }
  ggplot(df, aes(fill=factor(type), y=v,x=factor(label))) + geom_bar(position='fill', stat='identity', alpha = 0.9) +
    scale_fill_brewer(palette = "Greens", direction = -1) + theme_minimal() +
    geom_text(aes(label = round(v,digits=3)),size = 3, hjust = 0.5, vjust = 2, position = "stack", alpha = 0.5)
}

# plot rho for different values to find trend
# q is held fixed
rho_plot <- function(p,n,k) {
  q <- seq(0,1,by=0.05) 
  rho <- sapply(seq(0,1,by=0.05), function(q) rho(p[1],n,k,q))
  df <- data.frame(q, rho)
  df$pvalue <- p[1]
  for(i in 2:length(p)) {
    rho <- sapply(seq(0,1,by=0.05), function(q) rho(p[i],n,k,q))
    df2 <- data.frame(q, rho)
    df2$pvalue <- p[i]
    df <- rbind(df, df2)
  }
  print(df)
  ggplot(df, aes(q,rho)) + geom_line(aes(color=factor(pvalue))) + theme_light() +
    ylim(0,2)
}

## --------------- AGE DIST PLOT ------------------------

# plots one curve of the age distribution
cell_sim <- function(i,p,n,k,q,trials) {
  ages <- cell_ages(i,p,n,k,q,trials)/trials
  df <- as.data.frame(ages)
  df$t <- 1:length(ages)
  pl1 <- ggplot(df, aes(t, ages)) + geom_line()
  pl1
}

# plot several curves in the same plot
# p is a vector e.g c(0.1, 0.3,0.5)
age_plot_p <- function(i,p,n,k,q,trials) {
  theme_set(theme_minimal())
  age <- cell_ages(i,p[1],n,k,q,trials)/trials
  df <- as.data.frame(age)
  df$x <- 1:length(age)
  df$dataset <- c(rep(p[1], length(age)))
  
  #pl <- ggplot() + geom_line(df, mapping = aes(x,age))
  for(i in 2:length(p)) {
    age <- cell_ages(i,p[i],n,k,q,trials)/trials
    df2 <-as.data.frame(age)
    df2$x <- 1:length(age)
    df2$dataset <- c(rep(p[i], length(age)))
    df <- rbind(df, df2)
  }
  print(df)
  pl <- ggplot(df, aes(x, age)) + 
    geom_line(aes(color = factor(dataset), group = dataset), size = 1) 
  pl
}

# plots age distribution for different starting cells
# i should be a vector of integers e.g c(0,1)
age_plot_type <- function(i,p,n,k,q,trials) {
  theme_set(theme_minimal())
  age <- cell_ages(i[1],p,n,k,q,trials)/trials
  df <- as.data.frame(age)
  df$x <- 0:(length(age)-1)
  df$dataset <- c(rep(i[1], length(age)))
  # fult men funkar
  for(ii in 2:length(i)) {
    age <- cell_ages(i[ii],p,n,k,q,trials)/trials
    df2 <-as.data.frame(age)
    df2$x <- 0:(length(age)-1)
    df2$dataset <- c(rep(i[ii], length(age)))
    df <- rbind(df, df2)
  }
  col<- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  pl <- ggplot(df, aes(x, age)) + 
    geom_line(aes(color = factor(dataset),group = dataset), size = 1)
  pl + scale_color_manual(values = col)
  # TODO: fixa x-axeln
}

## ---------------- HELPER FUNCTION ---------------------

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


### ------------------ SPECIAL CASE -------------------

# plots the critical value in the (k,n) plane
# PARAMS: vector of p values ex: c(0.1,0.3,0.5), k = 1:15
critical_pl <- function(p,k) {
  df <- critical_df(p[1],k)
  df$pvalue <- p[1]
  for(i in 2:length(p)) {
    df2 <- critical_df(p[i],k)
    df2$pvalue <- p[i]
    df <- rbind(df, df2)
  }
  ggplot(df, aes(k, crit, color = factor(pvalue))) +
    geom_line(size = 1) +
    ggtitle("(Sub)kritiskt $n$ värde") + 
    xlab("$k$") + ylab("$n$") + theme_light() 
}

## ------------------- HELPER FUNCTIONs ---------------

# finds the first n for every k where the the process is subkritical
critical_n <- function(p, k) {
  rho <- 100
  n <- 0
  while(rho > 1) {
    n <- n+1
    rho <- pf_eigen(M_mat(p,n,k,1,1))[[1]]
  }
  return(n)
}

# converts to dataframe
critical_df <- function(p,k) {
  crit <- sapply(k, function(k_) critical_n(p,k_))
  df <- as.data.frame(crit)
  df$k <- k
  return(df)
}
