library(expm)
library(ggplot2)
library(tidyverse)
library(pracma)
library(patchwork)
library(MASS)
library(easyGgplot2)
library(tikzDevice)
library(EnvStats)

rejuv_mean_plot <- function(p,n,k,q,trials) {
  init <- 0:k
  df <- lapply(p, function(p_) {
    data.frame(
    mean = sapply(init, function(i_) mean(df_drls(i_,p_,n,k,trials))),
    type = init,
    p_value = p_
    )
  }
  )
  df <- do.call("rbind", df)
  print(df)
  ggplot(df, aes(x = type, y = mean, group = p_value)) + 
    geom_line(aes(color=factor(p_value))) + geom_point(aes(color=factor(p_value))) +
    scale_color_brewer(palette="Paired") + labs(color = 'p')
}

### ------------------ GENERAL CASE -------------------

# plot stable type distribution
multi_gw_pl <- function(p, n, k, q, hours, Z_0, trials) {
  
  theme_set(theme_minimal())
  
  df_wide <- multi_gw_mean(p, n, k, q, hours, Z_0, trials)
  df_long <- to_long(rowSums(df_wide))
  
  M_analytic <- M_mat(p, n, k, q,1)
  stable_dist <- round(pf_eigen(M_analytic)[[2]], digits = 3)
  print(stable_dist)
  
  type_long <- to_long(type_frequency(df_wide))
  pl1 <- ggplot(df_long, aes(x, Size)) + geom_line() # +
    # xlab("Tidsenhet") + ylab("Populationsstorlek")
  pl2 <- ggplot(type_long, aes(x, Size)) + 
    geom_line(aes(color = Type, group = Type)) + 
     geom_hline(yintercept = stable_dist, color = 'grey',size=0.6) +
     scale_y_continuous(breaks = stable_dist)
  pl1  + pl2 
}

# Plot the population size over time
multi_gw_size <- function(p,n,k,q,hours,Z_0,trials){
  df <- data.frame(
    sapply(p, function(p_) rowSums(multi_gw_mean(p_, n, k, q, hours, Z_0, trials))),
    t = 0:hours
  )
  colnames(df)[1:length(p)] <- p 
  df <- df %>% gather(key = 'p', value = 'size', -t)
  ggplot(df, aes(x = t, y = size, color = factor(p))) + geom_line() +
    theme_light()
}

# stacked bar plots for stable type distribution
v_plot <- function(p,n,k,q) {
  # check supercriticality
  rhos <- sapply(p, function(p_) rho(p_,n,k,q))
  if(min(rhos) <= 1) {
    print("Subcritical!")
    return(0)
  } 
  df <- data.frame(
    sapply(p, function(p_) st_type(p_,n,k,q)),
    type = 0:k
  )
  colnames(df)[1:length(p)] <- p 
  df <- df %>% gather(key = 'p', value = 'v', -type)
  ggplot(df, aes(fill=factor(type), y=v,x=factor(p))) + geom_bar(position='fill', stat='identity', alpha = 0.9) +
    scale_fill_brewer(palette = "Greens", direction = -1) + theme_minimal() 
}

# plot rho for different values to find trend
# q is held fixed
rho_plot <- function(p,n,k) {
  q <- seq(0,1,by=0.05) 
  df <- data.frame(
    sapply(p, function(p_) sapply(q, function(q_) rho(p_,n,k,q_))),
    q = q
  )
  colnames(df)[1:length(p)] <- p 
  df <- df %>% gather(key = 'p', value = 'rho', -q)
  ggplot(df, aes(q, rho)) + geom_line(aes(color=factor(p))) + theme_light() + ylim(0,2)
}

rho_plot_sc <- function(n,k) {
  p <- seq(0.01,0.5,0.025)
  df <- data.frame(
    rho = sapply(n, function(n_) sapply(p, function(p_) rho(p_,n_,k,1))),
    p = p
  )
  colnames(df)[1:length(n)] <- n 
  df <- df %>% gather(key = 'n', value = 'rho', -p)
  ggplot(df, aes(p, rho)) + geom_line(aes(color=factor(n))) + theme_light() + ylim(1,2)
}

## _______________ REJUVENATION ________________________

rejuv_plot <- function(p,n,k,q,trials) {
  init <- 0:k
  df <- lapply(p, function(p_) {
    data.frame(
    mean = sapply(init, function(i_) mean(df_drls(i_,p_,n,k,trials))),
    type = init,
    p_value = p_
    )
  }
  )
  df <- do.call("rbind", df)
  print(df)
  ggplot(df, aes(x = type, y = mean, group = p_value)) + 
    geom_line(aes(color=factor(p_value))) + geom_point(aes(color=factor(p_value))) +
    scale_color_brewer(palette="Paired") + labs(color = 'p')
}

## --------------- AGE DIST PLOT ------------------------

# Plot several geometric distributions in the same plot
# Parameter g independent of initial biological age
# --> plot for different values of p


# plots one curve of the age distribution
cell_sim <- function(i,p,n,k,q,trials) {
  ages <- cell_ages(i,p,n,k,q,trials)/trials
  df <- as.data.frame(ages)
  df$t <- 1:length(ages)
  pl1 <- ggplot(df, aes(t,ages)) + geom_bar(stat="identity") 
  pl1
}

# plot several curves in the same plot
age_plot_p <- function(i,p,n,k,q,trials) {
  make_df <- function(i,p,n,k,q,trials) {
    prop = cell_ages(i,p,n,k,q,trials)/trials
    df <- data.frame(
      prop = prop,
      p_value = p,
      x = 1:length(prop)
    )
    return(df)
  }
  df <- do.call("rbind", lapply(p, function(p_) make_df(i,p_,n,k,q,trials)))
  ggplot(df, aes(x = x, y = prop)) + 
    geom_line(aes(color = factor(p_value), group = p_value), size = 1) 
}

# plots age distribution for different starting cells
# i should be a vector of integers e.g c(0,1)
age_plot_type <- function(i,p,n,k,q,trials) {
  rep_vec <- function(i,p,n,k,q,trials) {
    data <- cell_ages(i,p,n,k,q,trials)
    rep(1:length(data), times = data)
  }
  make_df <- function(i,p,n,k,q,trials) {
    df <- data.frame(age = rep_vec(i,p,n,k,q,trials), type = i)
    return(df)
  }
  df <- do.call("rbind", lapply(i, function(i_) make_df(i_,p,n,k,q,trials)))
  #ggplot(df, aes(x = type, y = age))
  ggplot(df, aes(x = type, y = age, group=type)) + geom_boxplot() + coord_flip() +
    stat_summary(fun=mean, geom="point", shape=4, size=5, color="red", fill="red") 
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

tikz(file="plot_test.tex", width = 5, height = 4)
#plot <-  critical_pl(c(0.1,0.3,0.5), 1:15)
print(plot)
dev.off()

### ------------------ SPECIAL CASE -------------------

# Idea maximal nr of proteins a multi-GW process can accumulate
# without being subcritical (rho \leq 1)
critical_n <- function(p,k) {
  df <- lapply(p, function(p_)
    data.frame(
      k = k,
      n = sapply(k, function(k_) find_max_n(p_,k_)),
      p = p_
    )
  )
  df <- do.call("rbind", df)
  print(df)
  ggplot(df, aes(x = k,y = n,color = factor(p))) + geom_point() +
    scale_color_brewer(palette="Paired") + geom_line(aes(color=factor(p))) +
    theme_light()
}

## ------------------- HELPER FUNCTIONs ---------------

# find maximal n where rho > 1
find_max_n <- function(p,k) {
  rho <- Inf
  n <- 0
  while(rho > 1) {
    n <- n+1
    rho <- rho(p,n,k,1)
  }
  return(n-1)
}
