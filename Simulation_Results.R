
# First run the file "Functions.R" to load the packages and functions needed in this File


#### 4.8 Determining the Number of Iterations ####
# Function to vary the number of iterations
  vary_niter <- function(p_C, p_E, n_pilot = 1000, r = 1, niter_vec = c(100, 500, 1000, 2500, 5000, 7500, 10000, 15000, 20000)){
    results <- vector(mode = "list", length = length(niter_vec))
    # Iterate over the desired numbers of iterations to vary and save simulation result
    for (i in seq_along(niter_vec)) {
      sim <- simulation(p_C, p_E, r = r, n_pilot = n_pilot, niter = niter_vec[i])
      results[[i]] <- sim
    }
    return(results)
  }
  
# Apply function
  sim_niter <- vary_niter(p_C, p_E, n_pilot = 1000)
  
# Plot the Simulations
  sim_plots(sim_niter, c(100, 500, 1000, 2500, 5000, 7500, 10000, 15000, 20000), same_scale = FALSE, plots=1)
  
# Summary statistics
  sum_iter <- summary_statistics(sim_niter, group = c(100, 500, 1000, 2500, 5000, 7500, 10000, 15000, 20000), d_label = "iterations")
  print(xtable(sum_iter, caption=1, digits=4), include.rownames = FALSE)
  
# Scatterplot for different pilot study sample sizes
  # Function to vary pilot study sample size
  vary_npilot <- function(p_C, p_E, niter = 10000, r = 1, npilot_vec = c(50, 100, 200, 500, 1000, 2500, 5000, 10000, 15000)){
    results <- vector(mode = "list", length = length(npilot_vec))
    # Iterate over the desired pilot study sample sizes to vary and save simulation result
    for (i in seq_along(npilot_vec)) {
      sim <- simulation(p_C, p_E, r = r, niter = niter, n_pilot = npilot_vec[i])
      results[[i]] <- sim
    }
    return(results)
  }
  
# Scatter plot with small number of iterations  
  sim_iter_small <- vary_npilot(p_C, p_E, r=1, niter = 1000,  npilot_vec = c(seq(150,1000,50), seq(1500, 12000, 500)))
  sim_scatter_plot(data=sim_iter_small, group = c(seq(150,1000,50), seq(1500, 12000, 500)), 
                   xlabel="Pilot Study Sample Size")
  
  
# Scatter plot with small number of iterations
  sim_iter_large <- vary_npilot(p_C, p_E, r=1, niter = 10000,  npilot_vec = c(seq(150,1000,50), seq(1500, 12000, 500)))
  sim_scatter_plot(data=sim_iter_large, group = c(seq(150,1000,50), seq(1500, 12000, 500)), 
                   xlabel="Pilot Study Sample Size")
  
# Plot them together
  grid.arrange(sim_scatter_plot(data=sim_iter_small, group = c(seq(150,1000,50), seq(1500, 12000, 500)), 
                                xlabel="Pilot Study Sample Size"),
               sim_scatter_plot(data=sim_iter_large, group = c(seq(150,1000,50), seq(1500, 12000, 500)), 
                                xlabel="Pilot Study Sample Size"), 
               ncol=2)
  
#### 5.1 Effect of the pilot study sample size ####
  #### Results ####
  # Function to vary pilot study sample size introduced in the Code to chapter 4.8
  
  # Apply function
  sim_npilot <- vary_npilot(p_C, p_E)
  
  # Plot the simulations
  sim_plots(sim_npilot, c(50, 100, 200, 500, 1000, 2500, 5000, 10000, 15000), plots=0)
  
  # Summary Statistics
  sum_npilot <- summary_statistics(sim_npilot, group =  c(50, 100, 200, 500, 1000, 2500, 5000, 10000, 15000), d_label = "n_pilot")
  print(xtable(sum_npilot, caption=1, digits=4), include.rownames = FALSE)
  
  # Scatter plot of mean actual powers
  sim_scatter_plot(data=sim_npilot, group = c(50, 100, 200, 500, 1000, 2500, 5000, 10000, 15000), 
                   xlabel="Pilot Study Sample Size")
  
    # Median difference between smallest and biggest sample size in every iteration
  diff_table_npilot <- median_diff(sim_npilot, c(50, 100, 200, 500, 1000, 2500, 5000, 10000, 15000))
  print(xtable(diff_table_npilot), include.rownames = FALSE)
  
  # Plot the mean minimum sample size and the median min sample size with quartiles
  sim_ss_plots(data=sim_npilot, 
               group = c(50, 100, 200, 500, 1000, 2500, 5000, 10000, 15000),
               xlabel="Pilot Study Sample Size", plots=2)
  
  ##### Explanation for over- and underpowering in small pilot studies ####
  
  ### leave out study when effect is very small
  # Same as normal simulation function, but with step to calculate the difference between the true and estimated treatment effect
  simulation_diff_effect<-function(p_C, p_E, r=1, n_pilot, niter=1000, alpha=0.05,beta=0.2, effect_cut){
    n_needed<-matrix(NA,niter,4)
    colnames(n_needed)<-c("WMW","ttestord", "po", "index")
    actual_power2<-matrix(NA,niter,4)
    colnames(actual_power2)<-c("WMW","ttestord", "po", "index")
    theta_p <- calculate_theta_A(p_C, p_E, r)[1,1] 
    for (i in 1:niter){
      set.seed(i)
      print(i)
      phat<-generate.pEpC.estimate(p_C=p_C, p_E=p_E, r=r, n_pilot=n_pilot)
      
      # Calculate log odds ratio of simulated pilot probability vectors
      theta_phat <- calculate_theta_A(p_C=phat$C, p_E=phat$E, r)[1,1]
      
      # Calculate difference between simulated pilot theta and true theta
      theta_diff <- abs(theta_phat-theta_p)
      
      # If the difference is too big, then skip the iteration
      if(theta_diff>effect_cut){    print ("effect too small") ; next }
      
      # Calculate sample sizes
      result_WMW<-samplesize_WMW(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
      result_ttestord<-samplesize_ttestord(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
      result_po <- samplesize_po_NN(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
      n_needed[i,1]<- result_WMW$n_total
      n_needed[i,2]<- result_ttestord$n_total
      n_needed[i,3]<- result_po$n_total
      n_needed[i,4]<- i
      actual_power2[i,1]<-result_WMW$actual_power2
      actual_power2[i,2]<-result_ttestord$actual_power2
      actual_power2[i,3]<-result_po$actual_power2
      actual_power2[i,4]<-i
    }
    n_needed <- as.data.frame(na.omit(n_needed))
    actual_power2 <- as.data.frame(na.omit(actual_power2))
    nmin<-apply(n_needed,FUN=min,MARGIN=1)
    whichnmin<-apply(n_needed[,c(1,2,3)],FUN=which.min,MARGIN=1)
    actual_power_nmin<-actual_power2[cbind(1:length(whichnmin),whichnmin)]
    list(n_needed=n_needed,actual_power = actual_power2, 
         nmin=nmin,actual_power_nmin=actual_power_nmin,
         method = c("WMW", "ttest", "PO")[whichnmin])
  }
  
  # Apply Simulation function
  test_eff_cut <- simulation_diff_effect(p_C, p_E, n_pilot=100, niter=10000, effect_cut = 0.15)
  length(test_eff_cut$actual_power_nmin) # 3171 iterations estimated accurate enough treatment effects
  mean(test_eff_cut$actual_power_nmin)

  # Save in Data Frame
  df_eff_cut <- data.frame(power = test_eff_cut$actual_power_nmin) %>%
    mutate(mean = mean(power),
           type = "Cut") 
  
  ### Simulation for Comparison
  sim100 <- simulation(p_C,p_E,n_pilot=100,niter=10000)
  df_no_cut <- data.frame(power_no_cut = sim100$actual_power_nmin) %>%
    mutate(mean = mean(power_no_cut),
           type = "No Cut") %>%        
    rename(power = power_no_cut)            
  
  
  # Combine Datasets
  df_eff <- bind_rows(df_no_cut, df_eff_cut)
  df_eff$type <- as.factor(df_eff$type)
  df_eff$type <- relevel(df_eff$type, "No Cut")

  # Plot 
  ggplot(df_eff) + 
    geom_histogram(aes(x = power, fill = type, colour = type, alpha = type), position = "identity")+
    geom_vline(aes(xintercept = mean, alpha = type, colour = type), alpha = 0.5)+
    geom_vline(aes(xintercept = 0.8, linetype = "80% Power", colour = "No Cut"), colour = "red")+
    scale_alpha_manual(name = "Method", values = c("Cut"=1, "No Cut"=0.3))+
    scale_colour_manual(name= "Method", values = c("Cut" = "black", "No Cut" = "red")) +
    scale_fill_manual(name= "Method", values = c("Cut" = "grey", "No Cut" = "red")) +
    scale_linetype_manual(name = "", values = c("80% Power" = "dotted")) +
    xlab("Actual Power of Minimum Sample Size")+
    ylab("Count")+
    theme_bw()+
    theme(legend.spacing.y = unit(0, "cm"), legend.margin=margin(t=0, r=0.5, b=0, l=0.5, unit="cm")) +
    guides(fill = guide_legend(order = 1), colour = guide_legend(order = 1), alpha = guide_legend(order = 1),
           linetype = guide_legend(order = 2))

  

  ### leave out pilot study when calculated sample size is too small or too big to try and get rid of under and overpowering
  # Same as normal simulation function, but with step to skip iteration if the minimum sample size is too small or large
  simulation_samp_range <- function(p_C, p_E, r=1, n_pilot, niter=1000, alpha=0.05,beta=0.2, min_samp = 100, max_samp=1000){
    n_needed<-matrix(NA,niter,4)
    colnames(n_needed)<-c("WMW","ttestord", "po", "index")
    actual_power2<-matrix(NA,niter,4)
    colnames(actual_power2)<-c("WMW","ttestord", "po", "index")
    for (i in 1:niter){
      set.seed(i)
      print(i)
      # Generate a sampled pilot study and the estimated probability vectors
      phat<-generate.pEpC.estimate(p_C=p_C, p_E=p_E, r=r, n_pilot=n_pilot)
      
      # Calculate sample sizes
      result_WMW<-samplesize_WMW(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
      result_ttestord<-samplesize_ttestord(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
      result_po <- samplesize_po_NN(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
      n_min_samp <- min(c(result_WMW$n_total, result_ttestord$n_total, result_po$n_total))
      
      # Stop iteration when the calculated sample size is out of range
      if((n_min_samp<min_samp) | (n_min_samp>max_samp)){    print ("sample size not in range") ; next }
      
      n_needed[i,1]<- result_WMW$n_total
      n_needed[i,2]<- result_ttestord$n_total
      n_needed[i,3]<- result_po$n_total
      n_needed[i,4]<- i
      actual_power2[i,1]<-result_WMW$actual_power2
      actual_power2[i,2]<-result_ttestord$actual_power2
      actual_power2[i,3]<-result_po$actual_power2
      actual_power2[i,4]<-i
    }
    n_needed <- as.data.frame(na.omit(n_needed))
    actual_power2 <- as.data.frame(na.omit(actual_power2))
    actual_power2 <- actual_power2[n_needed$index,]
    nmin<-apply(n_needed,FUN=min,MARGIN=1)
    whichnmin<-apply(n_needed[,c(1,2,3)],FUN=which.min,MARGIN=1)
    actual_power_nmin<-actual_power2[cbind(1:length(whichnmin),whichnmin)]
    list(n_needed=n_needed,actual_power = actual_power2, 
         nmin=nmin,actual_power_nmin=actual_power_nmin,
         method = c("WMW", "ttest", "PO")[whichnmin])
  }
  
  sim_samp_range <- simulation_samp_range(p_C, p_E, r=1, niter = 10000, n_pilot = 100, min_samp = 100, max_samp=1000)
  mean(sim_samp_range$actual_power_nmin, na.rm = TRUE)
  length(sim_samp_range$actual_power_nmin)
  
  # Save in Data Frame
  df_samp_cut <- data.frame(power = sim_samp_range$actual_power_nmin) %>%
    mutate(mean = mean(power, na.rm = TRUE),
           type = "Cut") 
  
  # Combine Datasets
  df_samp <- bind_rows(df_no_cut, df_samp_cut)
  df_samp$type <- as.factor(df_samp$type)
  df_samp$type <- relevel(df_samp$type, "No Cut")

  #  Plot the distribution of the cut off simulation and of the un-cut simulation 
  ggplot(df_samp) + 
    geom_histogram(aes(x = power, fill = type, colour = type, alpha = type), position = "identity")+
    geom_vline(aes(xintercept = mean, alpha = type, colour = type), alpha = 0.5)+
    geom_vline(aes(xintercept = 0.8, linetype = "80% Power", colour = "No Cut"), colour = "red")+
    scale_alpha_manual(name = "Method", values = c("Cut"=1, "No Cut"=0.3))+
    scale_colour_manual(name= "Method", values = c("Cut" = "black", "No Cut" = "red")) +
    scale_fill_manual(name= "Method", values = c("Cut" = "grey", "No Cut" = "red")) +
    scale_linetype_manual(name = "", values = c("80% Power" = "dotted")) +
    xlab("Actual Power of Minimum Sample Size")+
    ylab("Count")+
    theme_bw()+
    theme(legend.spacing.y = unit(0, "cm"), legend.margin=margin(t=0, r=0.5, b=0, l=0.5, unit="cm")) +
    guides(fill = guide_legend(order = 1), colour = guide_legend(order = 1), alpha = guide_legend(order = 1),
           linetype = guide_legend(order = 2))
  
  
  ### Do this with milder limits
  sim_samp_range_mild <- simulation_samp_range(p_C, p_E, r=1, niter = 10000, n_pilot = 100, min_samp = 50, max_samp=1500)
  mean(sim_samp_range_mild$actual_power_nmin, na.rm = TRUE)
  
  # Save in Data Frame
  df_samp_cut_mild <- data.frame(power = sim_samp_range_mild$actual_power_nmin) %>%
    mutate(mean = mean(power, na.rm = TRUE),
           type = "Cut") 
  
  # Combine Datasets
  df_samp_mild <- bind_rows(df_no_cut, df_samp_cut_mild)
  df_samp_mild$type <- as.factor(df_samp_mild$type)
  df_samp_mild$type <- relevel(df_samp_mild$type, "No Cut")
  
  #  Plot the distribution of the cut off simulation and of the un-cut simulation 
  ggplot(df_samp_mild) + 
    geom_histogram(aes(x = power, fill = type, colour = type, alpha = type), position = "identity")+
    geom_vline(aes(xintercept = mean, alpha = type, colour = type), alpha = 0.5)+
    geom_vline(aes(xintercept = 0.8, linetype = "80% Power", colour = "No Cut"), colour = "red")+
    scale_alpha_manual(name = "Method", values = c("Cut"=1, "No Cut"=0.3))+
    scale_colour_manual(name= "Method", values = c("Cut" = "black", "No Cut" = "red")) +
    scale_fill_manual(name= "Method", values = c("Cut" = "grey", "No Cut" = "red")) +
    scale_linetype_manual(name = "", values = c("80% Power" = "dotted")) +
    xlab("Actual Power of Minimum Sample Size")+
    ylab("Count")+
    theme_bw()+
    theme(legend.spacing.y = unit(0, "cm"), legend.margin=margin(t=0, r=0.5, b=0, l=0.5, unit="cm")) +
    guides(fill = guide_legend(order = 1), colour = guide_legend(order = 1), alpha = guide_legend(order = 1),
           linetype = guide_legend(order = 2))
  
  ### Do this with stricter limits
  sim_samp_range_strict <- simulation_samp_range(p_C, p_E, r=1, niter = 10000, n_pilot = 100, min_samp = 150, max_samp=750)
  mean(sim_samp_range_strict$actual_power_nmin, na.rm = TRUE)
  
  # Save in Data Frame
  df_samp_cut_strict <- data.frame(power = sim_samp_range_strict$actual_power_nmin) %>%
    mutate(mean = mean(power, na.rm = TRUE),
           type = "Cut") 
  
  # Combine Datasets
  df_samp_strict <- bind_rows(df_no_cut, df_samp_cut_strict)
  df_samp_strict$type <- as.factor(df_samp_strict$type)
  df_samp_strict$type <- relevel(df_samp_strict$type, "No Cut")
  
  #  Plot the distribution of the cut off simulation and of the un-cut simulation 
  ggplot(df_samp_strict) + 
    geom_histogram(aes(x = power, fill = type, colour = type, alpha = type), position = "identity")+
    geom_vline(aes(xintercept = mean, alpha = type, colour = type), alpha = 0.5)+
    geom_vline(aes(xintercept = 0.8, linetype = "80% Power", colour = "No Cut"), colour = "red")+
    scale_alpha_manual(name = "Method", values = c("Cut"=1, "No Cut"=0.3))+
    scale_colour_manual(name= "Method", values = c("Cut" = "black", "No Cut" = "red")) +
    scale_fill_manual(name= "Method", values = c("Cut" = "grey", "No Cut" = "red")) +
    scale_linetype_manual(name = "", values = c("80% Power" = "dotted")) +
    xlab("Actual Power of Minimum Sample Size")+
    ylab("Count")+
    theme_bw()+
    theme(legend.spacing.y = unit(0, "cm"), legend.margin=margin(t=0, r=0.5, b=0, l=0.5, unit="cm")) +
    guides(fill = guide_legend(order = 1), colour = guide_legend(order = 1), alpha = guide_legend(order = 1),
           linetype = guide_legend(order = 2))

  
#### 5.2 Effect of the allocation ratio ####
  #### Visualise the sample sizes for different allocations ####
  # Function to calculate the group sample sizes of a total sample size for different allocation ratios
  all_rat <- function(n=100, r_vec=c(seq(0.1, 0.9, 0.2), 1, seq(2,7,1))){
    # Create empty vectors
    n_C_vec <- c()
    n_E_vec <- c()
    
    # Loop over the different allocation ratios
    for (r in seq_along(1:length(r_vec))) {
      # Calculate the group sample sizes and round them up
      n_E_vec[r] <- ceiling(r_vec[r]/(1+r_vec[r]) * n)
      n_C_vec[r] <- ceiling(1/(1+r_vec[r]) * n)
    }
    return(data.frame("n" = c(n_E_vec, n_C_vec), 
                      "r" = rep(r_vec, 2),
                      "group"=c(rep("n_E", length(n_E_vec)),rep("n_C", length(n_C_vec))) ) )
  }
  
  # Apply function
  r_vec <- all_rat()
  
  # Plot the sample sizes and corresponding allocation ratios
  ggplot(r_vec, aes(x=r, y=n, color=group))+
    geom_line( )+
    geom_point( )+
    scale_color_manual(values = c("n_C" = "#669933", "n_E" = "#336699")) +
    ylab("Sample Size")+
    xlab("Allocation Ratio")+
    theme_bw()
  

  #### Results ####
  # Function to vary the allocation ratio
  vary_r <- function(p_C, p_E, n_pilot = 10000, niter = 10000, r_vec){
    results <- vector(mode = "list", length = length(r_vec))
    # Iterate over the desired allocation ratios to vary and save simulation result
    for (i in seq_along(r_vec)) {
      sim <- simulation(p_C, p_E, n_pilot=n_pilot, niter=niter, r=r_vec[i])
      results[[i]] <- sim
    }
    return(results)
  }
  
  # Apply simulation function
  sim_r <- vary_r(p_C, p_E, n_pilot = 1000, r_vec=c(seq(0.1, 0.9, 0.2), 1, seq(2,7,1)))
  
  # Plot the simulation
  sim_plots(sim_r, c(seq(0.1, 0.9, 0.2), 1, seq(2,7,1)), plots=0)
  
  # Study the confidence intervals, to try and find out the reason for the higher power under r=0.5
  t.test(sim_r[[3]]$actual_power_nmin) # CI [0.7899;0.7944]
  t.test(sim_r[[4]]$actual_power_nmin) # CI [0.7877;0.7921]
  t.test(sim_r[[5]]$actual_power_nmin) # CI [0.7873;0.7917]
  t.test(sim_r[[6]]$actual_power_nmin) # CI [0.7880;0.7923]

  # Summary
  sum_r <- summary_statistics(sim_r, group = c(seq(0.1, 0.9, 0.2), 1, seq(2,7,1)), d_label = "r")
  print(xtable(sum_r, caption=1, digits=4), include.rownames = FALSE)
  
  # Scatter plot of mean actual powers
  sim_scatter_plot(data=sim_r, group = c(seq(0.1, 0.9, 0.2), 1, seq(2,7,1)), 
                   xlabel="Allocation Ratio")
  
  # Median difference between smallest and biggest sample size in every iteration
  diff_table_r <- median_diff(sim_r, c(seq(0.1, 0.9, 0.2), 1, seq(2,7,1)), power = FALSE)
  xtable(diff_table_r, caption=1)

  # Plot the mean minimum sample size and the median min sample size with quartiles
  sim_ss_plots(data=sim_r, group = c(seq(0.1, 0.9, 0.2), 1, seq(2,7,1)), 
               xlabel="Allocation Ratio", plots=2)

#### 5.3 Effect of the treatment effect size ####
  ## Function to vary effect size theta
  vary_p_theta <- function(n_pilot = 1000, niter = 10000, r_val = 1, cat = 6, theta_vec=log(seq(0.1, 2.5, 0.4))){
    p_results <- data.frame()
    results <- vector(mode = "list", length = length(theta_vec))
    # Iterate over the desired treatment effects to vary and save simulation result
    for (i in seq_along(theta_vec)) {
      set.seed(1)
      # Sample uniformly distributed control probability vector p_C
      p_C1 <- uniform_simplex(cat)[[1]] 
      # Calculate experimental vector p_E to have the specific treatment effect
      p_E1 <- calc_p_E(p_C1, theta_A = theta_vec[i])
      # Add random noise, so PO assumption isnt fulfilled and PO method doesnt have advantage
      noise <- runif(cat, min = 0.95, max = 1.05) 
      # Apply noise and normalize
      p_noisy <- p_E1 * noise  
      p_E1 <- p_noisy/sum(p_noisy)
      # Apply the simulation to the new probability vectors
      sim <- simulation(p_C=p_C1, p_E=p_E1, n_pilot=n_pilot, niter=niter, r=r_val)
      results[[i]] <- sim
    }
    return(results)
  }
  
  # What is the odds ratio for OR=1 after applying noise 
    set.seed(1)
    p_C1 <- uniform_simplex(4)[[1]] 
    p_E1 <- calc_p_E(p_C1, theta_A = log(1))
    set.seed(1)
    # Apply noise
    noise <- runif(4, min = 0.95, max = 1.05) 
    p_noisy <- p_E1 * noise  
    p_E1 <- p_noisy/sum(p_noisy)
    # calculate the odds ratio
    exp(calculate_theta_A(p_C1, p_E1)[1,1])
  
  # Apply simulation function
  sim_p_theta <- vary_p_theta(n_pilot = 1000, theta_vec = log(c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4)), cat=4)
 
  # Plot the simulations
  sim_plots(sim_p_theta, c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4), plots = 0)
  
  # Summary statistics
  sum_theta <- summary_statistics(sim_p_theta, group = c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4), d_label = "theta")
  print(xtable(sum_theta, caption=1, digits = 4), include.rownames = FALSE)
  
  # Scatter plot of mean actual powers
  sim_scatter_plot(data=sim_p_theta, group = c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4), 
                   xlabel="Treatment Effect")
  
  # Median difference between smallest and biggest sample size in every iteration
  sim_theta_diff <- median_diff(sim_p_theta, c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4), power=FALSE)
  xtable(sim_theta_diff, caption=1)
  
  # Plot the mean minimum sample size and the median min sample size with quartiles
  sim_ss_plots(sim_p_theta, c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4), "Odds Ratio", plots = 2)

#### 5.4 Effect of the number of categories ####
  # Simulation function
  vary_prob_length_PO <- function(n_pilot = 1000, niter = 10000, r = 1, seed = 1, theta=log(1.8)){
    vec_length<- c(3:14)
    results <- vector(mode = "list", length = length(vec_length))
    # Iterate over the desired numbers of categories to vary and save simulation result
    for (i in seq_along(vec_length)) {
      set.seed(seed)
      # Create uniformly distributed control probability vector p_C
      p_C_u <- uniform_simplex(vec_length[i])[[1]]
      # Calculate experimental probability with specific treatment effect
      p_E_u <- calc_p_E(p_C_u, theta_A = theta) # no noise applied this time, to make all simulations comparable
      sim_u <- simulation(p_C_u, p_E_u, n_pilot=n_pilot, niter=niter, r=r)
      results[[i]] <- sim_u
    }
    return(results)
  }
  
  # Odds ratio from p_C and p_E
    exp(calculate_theta_A(p_C, p_E)[1,1])
    
  # Apply simulation function with the same odds ratio as p_C and p_E have
  sim_cat_PO <- vary_prob_length_PO(n_pilot = 1000, theta = log(0.4821251))
  
  # Plot simulations
  sim_plots(sim_cat_PO, c(3:14), plots=0)
  
  # Summary statistics
  sum_cat <- summary_statistics(sim_cat_PO, group = c(3:14), d_label = "categories")
  print(xtable(sum_cat, caption=1, digits=4), include.rownames = FALSE)
  
  # Scatter plot of mean actual powers
  sim_scatter_plot(data=sim_cat_PO, group = c(3:14), xlabel="Number of Categories")
  
  # Median difference between smallest and biggest sample size in every iteration
  diff_table_cat <- median_diff(sim_cat_PO, c(3:14), power = FALSE)
  xtable(diff_table_cat, caption=1)
  
  # Plot the mean minimum sample size and the median min sample size with quartiles
  sim_ss_plots(sim_cat_PO, c(3:14), xlabel = "Number of Categories", plots = 2)
  
  

#### 5.5 Effect of normally distributed probabilities ####
  # Function to sample a vector with a normal distribution
  normal_vector <- function(cat, niter=10000, theta_A , seed = 1){
    set.seed(seed)
    # sample normally distributed values
    samp_norm <- rnorm(n = niter, mean = 0, sd = 1)
    samp_df <- data.frame(value = samp_norm) 
    # calculate the breaking points to categorize the values
    breaks_norm <- (max(samp_norm)-min(samp_norm))/(cat)
    # cut the values up into categories
    samp_df$category <- cut(samp_df$value, 
                            breaks=c(round(seq(min(samp_norm), max(samp_norm), breaks_norm),6)), 
                            labels=c(1:cat)) 
    p <- as.vector(prop.table(table(samp_df$category)))
      p2 <- calc_p_E(p, theta_A = theta_A)
      # Add random noise
      noise <- runif(cat, min = 0.95, max = 1.05) 
      # Apply noise and normalize
      p_noisy <- p2 * noise  
      p2 <- p_noisy/sum(p_noisy)
      p2 <- p2/sum(p2)

    return(list(p, p2))
  }
  
  # Simulation function
  vary_norm <- function(n_pilot = 1000, niter = 10000, r = 1, rep = 10, vec_length = 6, theta = NA){
    rep_vec <- c(1:rep)
    results <- vector(mode = "list", length = rep)
    # calculate the normally distributed probability vectors
    p_C_norm <- normal_vector(cat = vec_length, niter = niter, seed = 1, theta_A = theta)[[1]]
    p_E_norm <- normal_vector(cat = vec_length, niter = niter, seed = 1, theta_A = theta)[[2]]
    # Iterate over the desired number of repetitions to vary and save simulation result
    for (i in seq_along(1:rep)) {
      set.seed(i)
      sim <- simulation(p_C_norm, p_E_norm, r = r, niter = niter, n_pilot = n_pilot)
      results[[i]] <- sim
    }
    return(results)
  }
  
  # Apply simulation function 
  sim_norm <- vary_norm(n_pilot = 1000, rep=1, theta = log(0.4821251), vec_length = 4)
  
  # Plot simulation
  sim_plots(sim_norm, c(""), plots=0)
  
  # Summary statistics
  sum_norm <- summary_statistics(sim_norm, group =  c("norm"), d_label = "norm")
  print(xtable(sum_norm, caption=1, digits = 4), include.rownames = FALSE)
  
  # Median difference between smallest and biggest sample size in every iteration
  median_diff(sim_norm, c(""), power = FALSE)
  
  ### Comparison to a not-normally distributed simulation
  # Simulation with n_pilot=1000
  sim_npilot_1000 <- list(sim_npilot[[5]])
  
  # Plot simulation
  sim_plots(sim_npilot_1000, c(""), plots=0)
  
  # Summary Statistics
  sum_norm2 <- summary_statistics(sim_npilot_1000, group =  c("not-norm"), d_label = "not-norm") # summary statistics
  print(xtable(sum_norm2, caption=1, digits = 4), include.rownames = FALSE) 
  
  # Median difference between smallest and biggest sample size in every iteration
  median_diff(sim_npilot_1000, c(""), power = FALSE)
  
#### 5.6 Combined Factors ####
  ### Combination of treatment effect size and pilot study sample size
  # small pilot study
  sim_p_theta_npilot_small <- vary_p_theta(n_pilot = 100, theta_vec = log(c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4)), cat=4, r=1)
  sim_plots(sim_p_theta_npilot_small, c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4), plots = 0)
  
  # big pilot study
  sim_p_theta_npilot_big <- vary_p_theta(n_pilot = 10000, theta_vec = log(c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4)), cat=4, r=1)
  sim_plots(sim_p_theta_npilot_big, c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4), plots = 0)
  
  ### Combination of allocation ratio with other factors
  ## varying effect size
  # r = 4
  sim_p_theta_r4 <- vary_p_theta(n_pilot = 1000, theta_vec = log(c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4)), cat=4, r_val=4)
  sim_plots(sim_p_theta_r4, c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4), plots = 0)
  
  # r = 0.5
  sim_p_theta_r05 <- vary_p_theta(n_pilot = 1000, theta_vec = log(c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4)), cat=4, r_val =0.5)
  sim_plots(sim_p_theta_r05, c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4), plots = 0)
  
  ## varying pilot study effect size
  # r = 4
  sim_p_npilot_r4 <- vary_npilot(p_C, p_E, r=4)
  sim_plots(sim_p_npilot_r4, c(50, 100, 200, 500, 1000, 2500, 5000, 10000, 15000), plots=0)
  
  # r = 0.5
  sim_p_npilot_r05 <- vary_npilot(p_C, p_E, r=0.5)
  sim_plots(sim_p_npilot_r05, c(50, 100, 200, 500, 1000, 2500, 5000, 10000, 15000), plots=0)
  
  