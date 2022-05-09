  library(expm)
  library(ggplot2)
  library(tidyverse)
  library(pracma)
  library(patchwork)
  library(MASS)
  library(easyGgplot2)
  library(tikzDevice)
  library(EnvStats)
  
### ------------------ GENERAL CASE --------------------------
  
# plot stable type distribution an population size
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
  
  # Plot the population size over time for different p
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
      print("Subcritical!") # test for subcriticality
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
  
  # plot rho as a function of q \in [0,1]
  rho_plot <- function(p,n,k) {
    q <- seq(0,1,by=0.05) 
    df <- data.frame(
      sapply(p, function(p_) sapply(q, function(q_) rho(p_,n,k,q_))),
      q = q
    )
    colnames(df)[1:length(p)] <- p 
    df <- df %>% gather(key = 'p', value = 'rho', -q)
    ggplot(df, aes(q, rho)) + geom_line(aes(color=factor(p))) + theme_light() + ylim(0,2) +
      scale_color_brewer(palette="Paired")
  }
  
  ## ----------- REJUVENATION ________________________

  #Using E(d-m)/E(m) make sure to change line in def_drls so it doesnt divide by mothers age there.              
  rejuv_mean_plot2 <- function(p,n,k,q,trials) {
    init <- 0:k
    df <- lapply(p, function(p_) {
      mean_age <- function(i_,p_,n,k,trials) {
        ages <- cell_ages(i_,p_,n,k,1,trials)
        mean(rep(0:(length(ages)-1),times = ages))
      }
      data.frame(
        mean = sapply(init, function(i_) mean(df_drls(i_,p_,n,k,trials))/mean_age(i_,p_,n,k,trials)),
        type = init,
        p_value = p_
      )
    }
    )
    df <- do.call("rbind", df)
    print(df)
    ggplot(df, aes(x = type, y = mean, color=factor(p_value))) + 
      geom_line() + geom_point() +
      scale_color_brewer(palette="Paired", direction = -1) + labs(color = 'p') +
      scale_x_continuous(breaks=seq(0,k,1))+
      theme_minimal() + ylim(-1.4,1.2)
  }  
  
  # mean rejuvenation index weighted with stable type distribution
  rejuv_population <- function(p,n,k,trials) {
    init <- 0:k
    st <- st_type(p,n,k,1)
    mean_ages <- sapply(init, function(i_) {
      ages <- cell_ages(i_,p,n,k,1,trials)
      mean(rep(0:(length(ages)-1), times = ages))
    })
    # access element i+1 since R is indexed from 1
    EX <- sapply(init, function(i) mean(df_drls(i,p,n,k,trials))/mean_ages[i+1])
    sum(st*EX) # weighted mean
  }
  
  # visualisering rejuvenation
  rejuv_pop_plot <- function(n,k,q,trials) {
    df <- lapply(n, function(n_) 
      data.frame(
        n = n_,
        mean = sapply(seq(0.1,0.5,by=0.1), function(p_) rejuv_population(p_,n_,k,trials)),
        p = seq(0.1,0.5,by=0.1)
      )
    )
    df <- do.call("rbind", df)
    ggplot(df, aes(x = p, y = mean, color = factor(n))) + geom_point() + geom_line() +
      theme_light()
  }
  
  # fraction of rejuventated cells 
  prop_rejuv_cells <- function(p,n,k,trials) {
    init <- 0:k
    df <- lapply(p, function(p_) {
        cells <- lapply(init, function(i_) df_drls(i_,p_,n,k,trials))
        data.frame(
          p = p_,
          prop = sapply(init, function(i) sum(cells[[i+1]] > 0)/length(cells[[i+1]])),
          type = init
        )
        }
      )
    df <- do.call("rbind", df)
    print(df)
    ggplot(df, aes(x = type, y = prop, color=factor(p))) + geom_point() + geom_line() +
      ylim(0,1) + theme_light() + ggtitle('Fraction of rejuvenated cells') + 
      scale_x_discrete(limits=init)
  }
  
  # proportion of rejuvenated cell weighted with stable type
  weighted_rej_prop <- function(p,n,k,q,trials) {
    stable_dist <- st_type(p,n,k,q)
    EX <- sapply(0:k, function(i_) {
        cells <- df_drls(i_,p,n,k,trials)
        sum(cells < 0)/length(cells)
      }
    ) 
    sum(EX*stable_dist)
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
    make_df <- function(i,p,n,k,q,trials) {
      prop = cell_ages(i,p,n,k,q,trials)/trials
      df <- data.frame(
        prop = prop[2:length(prop)],
        type = i,
        x = 1:(length(prop)-1)
      )
      return(df)
    }
    df <- do.call("rbind", lapply(i, function(i_) make_df(i_,p,n,k,q,trials)))
    #ggplot(df, aes(x = type, y = age))
    ggplot(df, aes(x = x, y = prop,color=factor(type))) + geom_line() + xlim(1,NA) +
      theme_light() + scale_color_brewer(palette="Paired", direction = -1) + labs(color = 'p') 
  }
  
  mean_age <- function(i,p,n,k,q,trials) {
    data <- cell_ages(i,p,n,k,q,trials)
    x <- 0:(length(data)-1)
    mean(rep(x, times = data))
  }
  
  weighted_age <- function(p,n,k,q,trials) {
    stable_dist <- st_type(p,n,k,q)
    EX <- sapply(0:k, function(i_) mean_age(i_,p,n,k,q,trials))
    sum(EX*stable_dist)
  }
  
  # build the data frame to plot weighted curves
  type_ages_df <- function(p,n,k,q,trials){
    st <- st_type(p,n,k,q)
    ls <- lapply(0:k, function(i) st[i+1]*cell_ages(i,p,n,k,q,trials))
    max_length <- max(sapply(ls, length))
    ls <- do.call("cbind", lapply(ls, "length<-", max_length))
    ls[is.na(ls)] <- 0
    row_sums <- rowSums(ls)
    #print(mean(rep(1:length(row_sums),row_sums)))
    data.frame(t = 1:nrow(ls), freq = row_sums/((k+1)*trials), var=q, mean = rep_mean(row_sums)) 
  }
  
  # plot the weighted age curves 
  type_ages_plot <- function(p,n,k,q,trials) {
    df <- do.call("rbind", lapply(q, function(q_) type_ages_df(p,n,k,q_,trials)))
    print(df)
    ggplot(df, aes(t,freq,color=factor(var)))  +
      geom_line() +
      scale_color_brewer(direction = -1, palette="Paired") + theme_minimal() #+
    #geom_vline(aes(xintercept=mean,color=factor(var)),linetype='dashed')
  }
  
  ## ---------------- GEOMETRIC DIST ----------------------
  
  # approximate the parameter g in the geometric distribution
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
    df_dist <- data.frame(prop = data/sum(data), x = 1:length(data)) 
    df_tail <- data.frame(prop = tail/sum(tail), x = 0:(length(tail)-1))
    
    p1 <- ggplot(df_dist, aes(x, prop)) + geom_line() +
      geom_vline(xintercept = strtoi(start), linetype="dashed", color = "red") 
    p2 <- ggplot(df_tail, aes(x,prop)) + geom_line() +
      geom_line(aes(x, dgeom(x,mle)), color = "red") + xlab("Försök")
    p1+p2
  }
  
  MLE_data <- function(i,p,n,k,q,trials) {
    data <- cell_ages(i,p,n,k,q,trials)
    plot(0:(length(data)-1), data/sum(data), type = 'l', ylab = '', xlab = 'Ålder')
    start <- readline(prompt = "Enter starting value: ")
    tail <- data[-(1:start)]
    MLE_geom(tail)
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
  
  tikz(file="plot_test.tex", width = 6, height = 5)
  #plot <-  critical_pl(c(0.1,0.3,0.5), 1:15)
  print(plot)
  dev.off()
  
