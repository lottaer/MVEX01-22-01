library(expm)
library(ggplot2)
library(tidyverse)
library(pracma)
library(patchwork)
library(MASS)
library(easyGgplot2)
library(nleqslv)
library(tikzDevice)

#! In this version the parameter q is ommited
# PARAMS:
# (p,n,k) as always
# trials = number of times the simulation is repeated
# hours = number of timesteps
# Z_0 startvector ex: c(1,0,0) 

# Om något paket saknas laddas det ned genom install.packages(...)
# Ladda paket som du laddat ned genom library(...)
# Alla funktioner och definitoner av variabler görs i terminalen
# Ex skapa väntevärdesmatrisen: M <- M_mat(0.5,2,2,1)
# Ex kör funktion genom terminal: multi_gw_mean(0.5,2,2,20,c(1,0,0),1000)
# Om en variabel är definierad används denna ibland som input t.ex pf_eigen(M)
# där M är väntevärdesmatrisen som vi definierade på rad 21

# För att ladda alla funktioner i R-filen till workspaces så kör ctrl-A enter

# All simulering är skriven i C++ och ligger i filen multi_gw_ver2.cpp
# Denna simulering är bara en direkt översättning av tidigare R-simulering
# men skriven på sådant sätt att minnes-allokeringen är minimerad och således
# är programmet märkbart snabbare

# Klicka på Source i filen för att kompilera. Alla funktioner som det 
# står det export över i cpp filen och kommer dyka upp
# i workspaces som R funktionerna. Dessa funktioner går att använda
# på precis samma sätt som funktionerna i R-filen

# Simulation --------------------------------------------------------------

# simulation of the reproduction matrix
M_est <- function(p, n, k, hours, trials) {
  M <- Reduce('+', lapply(1:trials, function(i) M_est_new(p, n, k) %^% hours))/trials
  return(M)
}

# simulation of multi gw, takes mean of trials simulations
multi_gw_mean <- function(p, n, k, hours, Z_0, trials) {
  if(length(Z_0) != k+1) {
    print("Enter starting vector of correct length")
    return(0)
  }
  Reduce('+', lapply(1:trials, function(i) Z_mat(p, n, k, hours, Z_0)))/trials
}

# Analytic ----------------------------------------------------------------

# Mean matrix. Use hours = 1
M_mat <- function(p, n, k, hours) {
  M <- matrix(0, k+1, k+1)
  for(i in 0:k) {
    M[i+1,] <- M_row(i, p, n, k) 
  }
  return(M %^% hours)
}

# Rows in the mean matrix
M_row <- function(i, p, n, k) {
  M_r <- sapply(0:k, function(j) dbinom(i+n-j, n+i, p) + dbinom(j, n+i, p))
}

# calculates the perron-frobenius eigenvalue and eigenvector
# to get eigenvalue pf_eigen(M)[[1]] or pf_eigen(M)[[2]] for the vector
# M is obtained from M_mat ex: M <- M_mat(p,n,k,1)
pf_eigen <- function(M) {
  eigs = eigen(t(M))
  vals <- eigs$values
  vecs <- eigs$vectors
  ind <- which.max(Re(vals[abs(Im(vals)) < 1e-6]))
  return(list(Re(vals[ind]), c(normalize(Re(vecs[,ind])))))
}

# expected size of population after n hours
expected_sz <- function(p, n, k, hours, start) {
  M <- M_mat(p,n,k, hours)
  sum(start*M)
}

# finds the first n for every k where the the process is subkritical
critical_n <- function(p, k) {
  rho <- 100
  n <- 0
  while(rho > 1) {
    n <- n+1
    rho <- pf_eigen(M_mat(p,n,k,1))[[1]]
  }
  return(n)
}

# plots the critical value in the (k,n) plane
# PARAMS: vector of p values ex: c(0.1,0.3,0.5), k = 1:15
critical_pl <- function(p,k) {
  theme_set(theme_minimal())
  df <- critical_df(p[1],k)
  df$pvalue <- p[1]
  for(i in 2:length(p)) {
    df2 <- critical_df(p[i],k)
    df2$pvalue <- p[i]
    df <- rbind(df, df2)
  }
  ggplot(df, aes(k, crit, color = factor(pvalue))) +
    geom_line(size = 1.2) +
    ggtitle("Kritiskt n värde") + 
    xlab("k") + ylab("n") 
}

### HELPER FUNCTIONS ###

# converts to dataframe
critical_df <- function(p,k) {
  crit <- sapply(k, function(k_) critical_n(p,k_))
  df <- as.data.frame(crit)
  df$k <- k
  return(df)
}

# normalize st u*1 = 1
normalize <- function(vec) {
  return(vec/sum(vec))
}

# Visualization -----------------------------------------------------------

# plot function. plots population size and stable type distribution
multi_gw_pl <- function(p, n, k, hours, Z_0, trials) {
  
  theme_set(theme_minimal())
  
  df_wide <- multi_gw_mean(p, n, k, hours, Z_0, trials)
  print(df_wide)
  print(sum(df_wide[hours+1,]))
  df_long <- to_long(rowSums(df_wide))

  M_analytic <- M_mat(p, n, k, 1)
  stable_dist <- round(pf_eigen(M_analytic)[[2]], digits = 3)
  print(stable_dist)
  
  type_long <- to_long(type_frequency(df_wide))
  
  # geom_line(aes(color = Type, group = Type), size = 1.2)
  pl1 <- ggplot(df_long, aes(x, Size)) + geom_line(size = 0.8) +
    xlab("t") + ylab("Populationsstorlek")
  pl2 <- ggplot(type_long, aes(x, Size)) + 
    geom_line(aes(color = Type, group = Type), size = 1.2) + xlab("t") +
    labs(y = "Proportion") + geom_hline(yintercept = stable_dist, color = 'grey',size=0.8) +
    scale_y_continuous(breaks = stable_dist)
  
  pl1  + pl2
}

### HELPER FUNCTIONS ###
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

#tikz(file="plot_test.tex", width = 5, height =5)
#plot <-  critical_pl(c(0.1,0.3,0.5), 1:15)
#print(plot)
#dev.off()

# Age distribution --------------------------------------------------------

# Age distribution --------------------------------------------------------

# plots one curve of the age distribution
# mother cell until it dies
cell_sim <- function(i,p,n,k,trials) {
  ages <- cell_ages(i,p,n,k,trials)/trials
  #print(E(ages))
  df <- as.data.frame(ages)
  df$t <- 1:length(ages)
  pl1 <- ggplot(df, aes(t, ages)) + geom_line()
  pl1
}

# Expected value from observations
E <- function(probs) {
  vals <- 0:(length(probs)-1)
  sum(vals*probs)
}

# plot several curves in the same plot
# p is a vector e.g c(0.1, 0.3,0.5)
age_plot <- function(i,p,n,k,trials) {
  theme_set(theme_minimal())
  age <- cell_ages(i,p[1],n,k,trials)/trials
  df <- as.data.frame(age)
  df$x <- 1:length(age)
  df$dataset <- c(rep(p[1], length(age)))
  
  #pl <- ggplot() + geom_line(df, mapping = aes(x,age))
  for(i in 2:length(p)) {
    age <- cell_ages(i,p[i],n,k,trials)/trials
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
age_plot_type <- function(i,p,n,k,trials) {
  theme_set(theme_minimal())
  age <- cell_ages(i[1],p,n,k,trials)/trials
  df <- as.data.frame(age)
  df$x <- 1:length(age)
  df$dataset <- c(rep(i[1], length(age)))
  for(ii in 2:length(i)) {
    age <- cell_ages(i[ii],p,n,k,trials)/trials
    df2 <-as.data.frame(age)
    df2$x <- 1:length(age)
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

gamma_par <- function(Z) {
  fitdistr(rep(1:length(Z), times=Z), "gamma")$estimate
}

# plot gamma curve and abserved data
gamma_ML <- function(Z) {
  params <- gamma_par(Z)
  print(params)
  df <- as.data.frame(Z/sum(Z))
  colnames(df) <- 'proportion'
  df$x <- 1:length(Z)
  #print(exp_life(Z))
  ggplot(df, aes(x = x, y = proportion)) + geom_bar(stat="identity") +
    geom_function(fun = dgamma, args = list(shape = params[1], rate = params[2]), size = 1) 
}

exp_life <- function(Z) {
  params <- gamma_par(Z)
  params[1]/params[2] # expected value
}


