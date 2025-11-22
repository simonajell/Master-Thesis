# load these packages
library(MASS)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(xtable)

#### Basis values ####
p_C <- c(0.2, 0.3, 0.3, 0.2)
p_E <- c(0.1, 0.2, 0.4, 0.3)
#### Functions used for the t-test method ####
  # Function to calculate delta_A and sigma 
  calculate_delta_A_sigma<-function(p_C,p_E){
    # Set number of outcome levels
    K<-length(p_C)
    
    # Calculate the means
    mu_E<-sum((1:K)*p_E)
    mu_C<-sum((1:K)*p_C)
    
    # Calculate delta_A, which is the difference of means
    delta_A<-mu_E-mu_C  
    
    # Calculate sigma
    sigma_C<- sqrt(sum(p_C*((1:K)-mu_C)^2))
    sigma_E<- sqrt(sum(p_E*((1:K)-mu_E)^2))
    sigma<-0.5*(sigma_C+sigma_E)
    
    list(delta_A=delta_A,sigma=sigma)
  }
  
  # Function to calculate the needed sample size for the two-sample t-test using the Guenther & Schouten correction 
  # adapted from Kieser (2020)
  samplesize_ttestord<-  function(p_C, p_E, alpha, beta, r, p_C2=NULL, p_E2=NULL){
    # Calculate delta_A and sigma
    param<-calculate_delta_A_sigma(p_C,p_E)
    
    # Sample size formula adapted from Formula 3.11a to 3.11c from Kieser (2020)
    n_C <- ceiling((1+r)/r * (qnorm(1-alpha/2) + qnorm(1-beta))^2 *
                     (param$sigma/param$delta_A)^2 + (qnorm(1-alpha/2)^2) / (2*(1+r)))
    n_E <- ceiling((1+r) * (qnorm(1-alpha/2) + qnorm(1-beta))^2 *
                     (param$sigma/param$delta_A)^2 + r * (qnorm(1-alpha/2)^2) / (2*(1+r)))
    
    # Calculate the actual allocation ratio resulting from the calculated sample sizes 
    actual_r <- n_E/n_C
    
    # Calculate the total sample size
    n_total <- n_E + n_C
    
    # Calculate the power
    actual_power <- pnorm(1/((1+actual_r) * abs(param$sigma/param$delta_A)) *
                            sqrt(actual_r*(n_total-(qnorm(1-alpha/2)^2)/
                                             (2))) - qnorm(1-alpha/2))
    
    # Calculate the actual power (power using the true probability vector p_C2 and p_E2)
    actual_power2<-NA
    if (!is.null(p_C2)&!is.null(p_E2))
    {
      param2<-calculate_delta_A_sigma(p_C2,p_E2)
      actual_power2 <- pnorm(1/((1+actual_r) * abs(param2$sigma/param2$delta_A)) *
                               sqrt(actual_r*(n_total-(qnorm(1-alpha/2)^2)/
                                                (2))) - qnorm(1-alpha/2))
    }
    return(data.frame(n_total = n_total, n_E = n_E, n_C = n_C,
                      actual_r = actual_r, actual_power = actual_power, actual_power2=actual_power2))
  }
  

#### Functions used for the PO method ####
  # Function to estimate theta_A and V_A
  calculate_theta_A <- function(p_C, p_E, r = 1) {
    # Validate inputs
    if (length(p_C) != length(p_E)) {
      stop("p_C and p_E must have the same length (same number of outcome levels).")
    }
    if (abs(sum(p_C) - 1) > 1e-6 || abs(sum(p_E) - 1) > 1e-6) {
      stop("Probabilities p_C and p_E must each sum to 1.")
    }
    
    # Set number of outcome levels
    K <- length(p_C)  
    
    # Construct outcome levels and group labels
    outcome <- rep(c(1:K), times = 2)
    group <- c(rep(c("control"), times = K), rep(c("treatment"), times = K))
    
    # Calculate weights under the alternative hypothesis
    weights_alt <- c(((r * p_C) / (r + 1)), (p_E / (r + 1)))
    
    # Create synthetic dataset
    data_alt <- data.frame(
      outcome = factor(outcome, ordered = TRUE),
      group = factor(group, levels = c("control", "treatment")),
      weight = weights_alt
    )
    
    # Fit proportional odds model under the alternative
    model_alt <- suppressWarnings(
      polr(outcome ~ group, data = data_alt, weights = weight, method = "logistic", Hess = TRUE)
    )
    
    # Extract coefficient of treatment effect
    coef_alt <- -coef(summary(model_alt))["grouptreatment", "Value"]
    var_alt <- coef(summary(model_alt))["grouptreatment", "Std. Error"]
    alt <- data.frame("coef" = coef_alt, "Var" = var_alt)
    return(alt)
  }

  # Function to estimate theta_N and V_N
  calculate_theta_N <- function(p_C, p_E, r = 1) {
    # Validate inputs
    if (length(p_C) != length(p_E)) {
      stop("p_C and p_E must have the same length (same number of outcome levels).")
    }
    if (abs(sum(p_C) - 1) > 1e-6 || abs(sum(p_E) - 1) > 1e-6) {
      stop("Probabilities p_C and p_E must each sum to 1.")
    }
    
    # Set number of outcome levels
    K <- length(p_C)  
    
    # Calculate the expected probability vector
    p_i <- (r*p_C+p_E)/(r+1) 
    
    # Construct outcome levels and group labels
    outcome <- rep(c(1:K), times = 2)
    group <- c(rep(c("control"), times = K), rep(c("treatment"), times = K))
    
    # Calculate weights under the null hypothesis
    weights_null <- c(((r * p_i) / (r + 1)), (p_i / (r + 1)))
    
    # Create synthetic dataset
    data_null <- data.frame(
      outcome = factor(outcome, ordered = TRUE),
      group = factor(group, levels = c("control", "treatment")),
      weight = weights_null
    )
    
    # Fit proportional odds model under the alternative
    model_null <- suppressWarnings(
      polr(outcome ~ group, data = data_null, weights = weight, method = "logistic", Hess = TRUE)
    )
    
    # Extract coefficient of treatment effect
    coef_null <- -coef(summary(model_null))["grouptreatment", "Value"]
    var_null <- coef(summary(model_null))["grouptreatment", "Std. Error"]
    null <- data.frame("coef" = coef_null, "Var" = var_null)
    return(null)
  }

  # Function to calculate sample size for "NN" Method of White (2023)
  samplesize_po_NN <- function(p_C, p_E, alpha=0.05, beta=0.2, r, p_C2=NULL, p_E2=NULL){
    
    # Estimate theta_A and V_N
    theta_A <- calculate_theta_A(p_C, p_E, r)[1,1]
    Var_N <- calculate_theta_N(p_C, p_E, r)[1,2]^2
    
    # Sample size formula (Formula 3 from White)
    n <- (Var_N*((qnorm(1-alpha/2)+qnorm(1-beta))^2))/(theta_A^2)
    
    # Calculate sample size in the two groups
    n_E <- ceiling(r/(1+r) * n)
    n_C <- ceiling(1/(1+r) * n)
    
    # Calculate the total sample size
    n_total <- n_E + n_C
    
    # Calculate the actual allocation ratio resulting from the calculated sample sizes
    actual_r <- n_E/n_C
    
    # Calculate power
    actual_power <- pnorm(sqrt((n_total*(theta_A^2))/Var_N)-qnorm(1-alpha/2))
    
    # Calculate actual power (power using the true probability vector p_C2 and p_E2)
    actual_power2<-NA
    if (!is.null(p_C2)&!is.null(p_E2)){
      theta_A2 <- calculate_theta_A(p_C2, p_E2, r)[1,1]
      Var_N2 <- calculate_theta_N(p_C2, p_E2, r)[1,2]^2
      actual_power2 <- pnorm(sqrt((n_total*(theta_A2^2))/Var_N2)-qnorm(1-alpha/2))
    }
    return(data.frame(n_total = n_total, n_E = n_E, n_C = n_C,
                      actual_r = actual_r, actual_power = actual_power, actual_power2=actual_power2))
  }
  
  # Additionally the other two approaches were also implemented ("NA" and "AA")
  samplesize_po_NA <- function(p_C, p_E, alpha, beta, r, p_C2=NULL, p_E2=NULL){
    
    d <- calculate_theta_A(p_C, p_E, r)[1,1]
    V_A <- calculate_theta_A(p_C, p_E, r)[1,2]^2
    V_N <- calculate_theta_N(p_C, p_E, r)[1,2]^2
    
    n <- ((sqrt(V_N)*qnorm(1-alpha/2)+sqrt(V_A)*qnorm(1-beta))^2)/(d^2) # White, 2023 Formel 4 (NA)
    n_E <- ceiling(r/(1+r) * n)
    n_C <- ceiling(1/(1+r) * n)
    n_total <- n_E + n_C
    actual_r <- n_E/n_C
    actual_power <- pnorm(((sqrt(n_total*d^2))-(sqrt(V_N)*qnorm(1-alpha/2)))/(sqrt(V_A)))
    actual_power2<-NA
    if (!is.null(p_C2)&!is.null(p_E2)){
      d2 <- calculate_theta_A(p_C2, p_E2, r)[1,1]
      V_A2 <- calculate_theta_A(p_C2, p_E2, r)[1,2]^2
      V_N2 <- calculate_theta_N(p_C2, p_E2, r)[1,2]^2
      actual_power2 <- pnorm(((sqrt(n_total*d2^2))-(sqrt(V_N2)*qnorm(1-alpha/2)))/(sqrt(V_A2)))
    }
    return(data.frame(n_total = n_total, n_E = n_E, n_C = n_C,
                      actual_r = actual_r, actual_power = actual_power, actual_power2=actual_power2))
  }
  samplesize_po_AA <- function(p_C, p_E, alpha, beta, r, p_C2=NULL, p_E2=NULL){
    
    theta_A <- calculate_theta_A(p_C, p_E, r)[1,1]
    Var_A <- calculate_theta_A(p_C, p_E, r)[1,2]^2
    
    n <- (Var_A*((qnorm(1-alpha/2)+qnorm(1-beta))^2))/(theta_A^2) # White, 2023 Formel 5 (AA)
    n_E <- ceiling(r/(1+r) * n)
    n_C <- ceiling(1/(1+r) * n)
    n_total <- n_E + n_C
    actual_r <- n_E/n_C
    actual_power <- pnorm(sqrt((n_total*(theta_A^2))/Var_A)-qnorm(1-alpha/2))
    actual_power2<-NA
    if (!is.null(p_C2)&!is.null(p_E2)){
      theta_A2 <- calculate_theta_A(p_C2, p_E2, r)[1,1]
      Var_A2 <- calculate_theta_A(p_C2, p_E2, r)[1,2]^2
      actual_power2 <- pnorm(sqrt((n_total*(theta_A2^2))/Var_A2)-qnorm(1-alpha/2))
    }
    return(data.frame(n_total = n_total, n_E = n_E, n_C = n_C,
                      actual_r = actual_r, actual_power = actual_power, actual_power2=actual_power2))
  }
  
  # Sample size calculation by Whitehead (1993) using the formulas from Kieser (2020)
  samplesize_po_kieser <- function(p_C, p_E, alpha, beta, r, p_C2=NULL, p_E2=NULL){
    
    theta_A <- calculate_theta_A(p_C, p_E, r)[1,1]
    x = 0
    for (i in 1:length(p_C)){
      x = x + ((p_C[i] + r*p_E[i]) / (1 + r))^3
    }
    n <- (1+r)^2/r * 3*(qnorm(1-alpha/2) + qnorm(1-beta))^2 /
      (theta_A^2 * (1-x)) # Kieser Formel 4.8
    n_E <- ceiling(r/(1+r) * n)
    n_C <- ceiling(1/(1+r) * n)
    n_total <- n_E + n_C
    actual_r <- n_E/n_C
    actual_power <- pnorm(sqrt(n_total * actual_r * theta_A^2 * (1-x)/
                                 (3 * (1 + actual_r)^2)) - qnorm(1-alpha/2))
    actual_power2<-NA
    if (!is.null(p_C2)&!is.null(p_E2)){
      theta_A2 <- calculate_theta_A(p_C2, p_E2, r)[1,1]
      actual_power2 <- pnorm(sqrt(n_total * actual_r * theta_A2^2 * (1-x)/
                                    (3 * (1 + actual_r)^2)) - qnorm(1-alpha/2))
    }
    return(data.frame(n_total = n_total, n_E = n_E, n_C = n_C,
                      actual_r = actual_r, actual_power = actual_power, actual_power2=actual_power2))
  }
  
  
#### Functions used for the WMW method ####  
  # Function to calculate pi_A from Formula 4.4 from Kieser (2020)
  calculate_pi_A<-function(p_C,p_E){
    x <- 0
    p_C_cum <- cumsum(p_C)
    for (i in 2:length(p_E)){
      x <- x + p_E[i] * p_C_cum[i-1]
    }
    y <- 0
    for (i in 1:length(p_E)){
      y <- y + 0.5 * p_E[i] * p_C[i]
    }
    return(x+y)  
  }
  
  # Function to calculate the needed sample size for a form of Mann-Whitney test taking ties into account
  # adapted from Kieser (2020)
  samplesize_WMW <- function(p_C, p_E, alpha, beta, r, p_C2=NULL, p_E2=NULL){
    # Calculate pi_A
    pi_A <- calculate_pi_A(p_C,p_E)
    
    # Calculate sample size using Formula 4.7 from Kieser (2020)
    x <- 0
    for (i in 1:length(p_C)){
      x <- x + ((r*p_C[i] + p_E[i]) / (1 + r))^3
    } 
    n <- ((1+r)^2/r) * ((qnorm(1-(alpha/2)) + qnorm(1-beta))^2) / 12 *
      (1-x)/((pi_A - 0.5)^2) 
    
    # Calculate sample size in the two groups
    n_E <- ceiling(r/(1+r) * n)
    n_C <- ceiling(1/(1+r) * n)
    
    # Calculate the total sample size
    n_total <- n_E + n_C
    
    # Calculate the actual allocation ratio resulting from the calculated sample sizes
    actual_r <- n_E/n_C
    
    # Calculate the power
    actual_power <- pnorm(abs(pi_A-0.5) / (1 + actual_r)*
                            sqrt(n_total * actual_r * 12/(1-x)) - qnorm(1-alpha/2))
    
    # Calculate actual power (power using the true probability vector p_C2 and p_E2)
    actual_power2<-NA
    if (!is.null(p_C2)&!is.null(p_E2))
    {
      pi_A2<-calculate_pi_A(p_C2,p_E2)
      actual_power2 <- pnorm(abs(pi_A2-0.5) / (1 + actual_r)*
                               sqrt(n_total * actual_r * 12/(1-x)) - qnorm(1-alpha/2))
    }
    return(data.frame(n_total = n_total, n_E = n_E, n_C = n_C,
                      actual_r = actual_r, actual_power = actual_power, actual_power2=actual_power2))
  }

#### Simulation Functions ####
  # This function generates RCT data of total size n_pilot with allocation ratio r and 
  # probability vectors p_C and p_E in the control and experimental arms, respectively, and returns the 
  # resulting estimates for p_C and p_E.
  generate.pEpC.estimate<-function(p_C, p_E, r=1, n_pilot){
    # Allocation to the groups is sampled using the pilot study sample size and the allocation ratio
    y <- rbinom(n_pilot,1,r/(r+1))  
    
    # Sample Sizes in the groups are calculated
    n_E <- sum(y)
    n_C <- n_pilot - n_E
    
    # Probability vectors are drawn using the probability vectors p_C and p_E and the group sample sizes
    phat_C <- as.vector(rmultinom(1, size=n_C, prob=p_C)/n_C)
    phat_E <- as.vector(rmultinom(1, size=n_E, prob=p_E)/n_E)
    list(C=phat_C,E=phat_E)
  }
  
  # This function returns the distribution of the needed sample sizes calculated using methods ttest_ord, WMW and PO over niter pilot datasets
  simulation<-function(p_C, p_E, r=1, n_pilot, niter=1000, alpha=0.05,beta=0.2){
    # Create empty vectors and matrices
    n_needed<-matrix(NA,niter,3)
    colnames(n_needed)<-c("WMW","ttestord", "po")
    actual_power2<-matrix(NA,niter,3)
    colnames(actual_power2)<-c("WMW","ttestord", "po")
    
    # Loop over the number of iterations
    for (i in 1:niter){
      # Set a seed for reproducability and print it
      set.seed(i)
      print(i)
      
      # Generate a sampled pilot study and the estimated probability vectors
      phat<-generate.pEpC.estimate(p_C=p_C, p_E=p_E, r=r, n_pilot=n_pilot)
      
      # Calculate sample sizes and powers using the three methods and the samples probability vectors
        # the probability vectors p_C and p_E are the assumed true probability vectors
      result_WMW<-samplesize_WMW(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
      result_ttestord<-samplesize_ttestord(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
      result_po <- samplesize_po_NN(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E) # newly added
      
      # Save the sample size calculations
      n_needed[i,1]<- result_WMW$n_total
      n_needed[i,2]<- result_ttestord$n_total
      n_needed[i,3]<- result_po$n_total
      
      # Save the actual power calculations
      actual_power2[i,1]<-result_WMW$actual_power2
      actual_power2[i,2]<-result_ttestord$actual_power2
      actual_power2[i,3]<-result_po$actual_power2
    }
    
    # Save the minimum sample size from every iteration
    nmin<-apply(n_needed,FUN=min,MARGIN=1)
    
    # Save the method thats calculated the smallest sample size in every iteration
    whichnmin<-apply(n_needed,FUN=which.min,MARGIN=1)
    
    # Save the actual power of the smallest sample size in every iteration
    actual_power_nmin<-actual_power2[cbind(1:niter,whichnmin)]
    
    list(n_needed=n_needed,actual_power = actual_power2, 
         nmin=nmin,actual_power_nmin=actual_power_nmin,
         method = c("WMW", "ttest", "PO")[whichnmin])
  }
  
#### Visualization and Summary Functions ####
  # Function to plot histograms for multiple simulations (also separated by methods)
    # "data": should be a list of simulations from the "simulation" function
  sim_plots <- function(data, group, plots=0, same_scale=TRUE){
    df_data <- data.frame()
    
    # Loop over the number of simulation in "data" to create a dataset including all simulations
    for (i in seq_along(1:length(data))) {
      # Create a data frame including the actual powers of the minimum sample sizes, the corresponding methos and the index of the simulation
      df <- data.frame("power_nmin" = data[[i]]$actual_power_nmin, 
                       "method" = data[[i]]$method,
                       "group" = group[i])
      df_data <- rbind(df_data,df)
    }
    
    # Transform the data frame to include a column with the mean actual power
    df_data <- df_data %>% group_by(group) %>%  mutate(mean = mean(power_nmin))
    df_data$method <- as.factor(df_data$method)
    df_data$method <- relevel(df_data$method, "ttest")
    
    # Plot the histograms
    # plot1: histogram of the actual power distribution for every simulation
    # plot2: plot1 separated by the method used to calculate the smallest sample size
    # here one has the choice to have the same scales in every histogram or free scales in every simulation
    if(same_scale == TRUE){
      plot1 <- ggplot(df_data, aes(x = power_nmin))+ 
        geom_histogram(fill = "grey", color = "black", position = "identity")+
        geom_vline(aes(xintercept = mean, group = group), color = "red")+
        geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
        facet_wrap(group ~.) +
        xlab("Actual Power")+
        ylab("Count")+
        theme_bw()
      plot2 <- ggplot(df_data, aes(x = power_nmin, fill = method, colour = method))+ 
        facet_wrap(group ~.) +
        geom_histogram(alpha = 0.3, position = "identity")+
        scale_colour_manual(values = c("PO" = "#7CAE00", "ttest" = "#00008B", "WMW" = "#C77CFF")) +
        scale_fill_manual(values = c("PO" = "#7CAE00", "ttest" = "#00008B", "WMW" = "#C77CFF")) +
        geom_vline(aes(xintercept = mean, group = group), color = "red")+
        geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
        xlab("Actual Power")+
        ylab("Count")+
        labs(colour="Method", fill="Method")+
        theme_bw()
    } else if(same_scale == FALSE){
      plot1 <- ggplot(df_data, aes(x = power_nmin))+ 
        geom_histogram(fill = "grey", color = "black", position = "identity")+
        geom_vline(aes(xintercept = mean, group = group), color = "red")+
        geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
        facet_wrap(group ~., scales="free_y") +
        xlab("Actual Power")+
        ylab("Count")+
        theme_bw()
      plot2 <- ggplot(df_data, aes(x = power_nmin, fill = method, colour = method))+ 
        facet_wrap(group ~., scales="free_y") +
        geom_histogram(alpha = 0.3, position = "identity")+
        scale_colour_manual(values = c("PO" = "#7CAE00", "ttest" = "#00008B", "WMW" = "#C77CFF")) +
        scale_fill_manual(values = c("PO" = "#7CAE00", "ttest" = "#00008B", "WMW" = "#C77CFF")) +
        geom_vline(aes(xintercept = mean, group = group), color = "red")+
        geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
        xlab("Actual Power")+
        ylab("Count")+
        labs(colour="Method", fill="Method")+
        theme_bw()
    }
    if(plots==0){
      grid.arrange(plot1, plot2, ncol=2)
    } else if(plots==1){
      return(plot1)
    } else if(plots == 2){
      return(plot2)
    }
  }
  
  # Function to plot the mean actual powers in a scatter plot
  sim_scatter_plot <- function(data, group, xlabel="Varied Component"){
    mean_col <- c()
    
    # Loop over the number of simulations in "data" to create a vector with the mean actual powers of every simulation
    for(i in seq(1:length(data))){
      mean_col[i] <- mean(data[[i]]$actual_power_nmin)
    }
    
    # Store the mean actual powers and the varied component values in a data frame
    sim_df <- data.frame(xcol =group, "nmin_power"=mean_col)
    
    # Plot the mean actual powers on the y-axis and the varied component on the x-axis
    ggplot(sim_df, aes(x=xcol, y=mean_col)) +
      geom_point()+
      geom_hline(yintercept=0.8, color = "red", linetype = "dotted")+
      xlab(xlabel)+
      ylab("Mean Actual Power")+
      theme_bw()
  }
  
  # Function to plot the mean minimum sample size of every simulation
    # and the median minimum sample size with the quartiles for every simulation
  sim_ss_plots <- function(data, group, xlabel="Varied Component", plots=0){
    df_samp <- data.frame()
    
    # Loop over the number of simulations in "data" to store the minimum sample sizes of every simulation 
      # and the values of the varied component
    for (i in seq_along(1:length(data))) {
      df <- data.frame("nmin" = data[[i]]$nmin, 
                       "group" = group[i])
      df_samp <- rbind(df_samp,df)
    }
    
    # Plot the mean minimum sample size
    plot1 <- ggplot(df_samp) +
      geom_pointrange(mapping = aes(x=group,y=nmin),
                      stat = "summary",
                      fun = mean, size=0.3)+
      xlab(xlabel)+
      ylab("Mean Minimum Sample Size")+
      theme_bw()
    
    # Plot the median minimum sample sizes with the respective IQR
    plot2 <- ggplot(df_samp) +
      geom_pointrange(mapping = aes(x=group,y=nmin),
                      stat = "summary",
                      fun.min = function(z) {quantile(z,0.25)},
                      fun.max = function(z) {quantile(z,0.75)},
                      fun = median, size=0.3)+
      xlab(xlabel)+
      ylab("Median Minimum Sample Size with Quartiles")+
      theme_bw()
    
    # Choose to plot both or only one of the two plots
    if(plots==0){
      grid.arrange(plot1, plot2, ncol=2)
    } else if(plots==1){
      return(plot1)
    } else if(plots == 2){
      return(plot2)
    }
    
  }
  
  # Function to create data frame with the summary statistics including:
    # Mean and standard deviation of the actual power distribution
    # Proportions of methods
    # Mean, median and IQR of minimum sample size
  summary_statistics <- function(data, group, d_label="Varied Component"){
    # Create an empty data frame with all summary statistics
    result <- data.frame(d_label = group,
                         "mean"=c(1:length(data)),
                         "sd" = c(1:length(data)),
                         "t-test" =  c(1:length(data)),
                         "PO" =  c(1:length(data)),
                         "WMW" =  c(1:length(data)),
                         "mean minimum sample size" = c(1:length(data)),
                         "median minimum sample size" = c(1:length(data)),
                         "IQR: 25" = c(1:length(data)),
                         "IQR: 75" = c(1:length(data)))
    colnames(result)[[1]] <- d_label
    
    # Calculate mean actual power
    for (i in seq_along(1:length(data))) {
      result[i, 2] <- mean(data[[i]]$actual_power_nmin)
    }
    
    # Calculate standard deviation of actual power
    for (i in seq_along(1:length(data))) {
      result[i, 3] <- sd(data[[i]]$actual_power_nmin)
    }
    # Calculate the proportion of methods
    for (i in seq_along(1:length(data))) {
      result[i,4] <- length(which(data[[i]]$method=="ttest"))/100
      result[i,5] <- length(which(data[[i]]$method=="PO"))/100
      result[i,6] <- length(which(data[[i]]$method=="WMW"))/100
    }
    # Calculate mean minimum sample size
    for (i in seq_along(1:length(data))) {
      result[i,7] <- mean(data[[i]]$nmin)
    }
    # Calculate median minimum sample size
    for (i in seq_along(1:length(data))) {
      result[i,8] <- median(data[[i]]$nmin)
    }
    # Calculate IQR
    for (i in seq_along(1:length(data))) {
      result[i,9] <-   quantile(data[[i]]$nmin, c(0.25))
      result[i,10] <-   quantile(data[[i]]$nmin, c(0.75))
    }
    
    return(result)
  }
  
  # Function to calculate median difference of smallest vs biggest sample size in a iteration for every simulation
  median_diff <- function(data, group, power = TRUE){
    # Create empty data frame
    diff_df <- data.frame( "sim" = c(1:length(data)),
                           "median sample size diff"=c(1:length(data)),
                           "mean sample size diff" = c(1:length(data)),
                           "median actual power diff"=c(1:length(data)),
                           "mean actual power diff"=c(1:length(data))
    )
    # Loop over the number of simulations contained in "data
    for (i in seq_along(1:length(data))) {
      # Calculate the difference between the maximum and minimum sample size of every iteration
      diff_vec_samp <- c()
      for(j in seq_along(1:length(data[[1]]$nmin))){
        diff_vec_samp[j] <-max(data[[i]]$n_needed[j,]) - min(data[[i]]$n_needed[j,])
      }
      
      # Calculate the difference in the actual power of the minimum and maximum sample sizes
      diff_vec_power <- c()
      for(k in seq_along(1:length(data[[1]]$nmin))){
        diff_vec_power[k] <-abs(max(data[[i]]$actual_power[k,]) - min(data[[i]]$actual_power[k,]))
      }
      # Save the varied component index and the mean and median of the difference in actual powers
      diff_df[i,1] <- group[i]
      diff_df[i,2] <- median(diff_vec_samp)
      diff_df[i,3] <- mean(diff_vec_samp[is.finite(diff_vec_samp)])
      
      # Delete or keep the difference in actual power columns
      if(power == TRUE){
        diff_df[i,4] <- median(diff_vec_power, na.rm=TRUE)
        diff_df[i,5] <- mean(diff_vec_power, na.rm=TRUE)
      } else if(power == FALSE){
        diff_df <- diff_df[,-5]
        diff_df <- diff_df[,-4]
      }
    }
    return(diff_df)
  }
  
  
  
#### Function to calculate p_E with a treatment effect #####
  # Function to calculate p_E when p_C and theta_A are known
  calc_p_E <- function(p_C, theta_A){
  
    # Use Formula 8 in Section 2.5 to calculate the cumulative experimental probability vector
    p_e_1 <- rep(0, length(p_C))
    for (i in seq_along(1:length(p_C))) {
      p_e_1[i] <- sum(p_C[1:i])/(sum(p_C[1:i])+((1-sum(p_C[1:i]))*exp(-theta_A)))
    }
    
    # Subtract the elements of the vector
    p_e_2 <- rep(0, length(p_C))
    for (i in (2:length(p_C))) {
      p_e_2[1]<- p_e_1[1]
      p_e_2[i] <- (p_e_1[i]-p_e_1[(i-1)])
    }
    return(p_e_2)
  }
#### Function that generates a simplex vector  ####
  # Calculate a random simplex vector without zero categories
  random_simplex <- function(n) {
    success <- FALSE
    # a random vector is sampled, as long as all vector elements are not 0
    while (!success) {
      # Gamma Distribution
      x <- rgamma(n, shape = 1, rate = 1)  
      
      # Normalize vector
      x <- x / sum(x)
      
      # Check for success
      success <- all(x > 0)
    }
    return(x)
  }
  
#### Function that generates a simplex vector with a uniform distribution ####
  uniform_simplex <- function(n) {
    # Create the uniform vector
    uni <- rep(1/n, n)
    
    # Add random noise
    noise <- runif(n, min = 0.7, max = 0.9) 
    
    # Apply noise and normalize
    uni_noisy <- uni * noise  
    uni_noisy <- uni_noisy/sum(uni_noisy)
    return(list(uni, uni_noisy))
  }  
  
#### Function that generates two simplex vectors with a specific log odd ratio ####
  generate_two_simplex_vectors <- function(n, theta = log(1.8)) {
    # Vector A: Gamma Distribution
    a <- rgamma(n, shape = 1, rate = 1)
    a <- a / sum(a)
    
    # Vector B: calculated with log odds ratio and random noise
    b <- calc_p_E(a, theta_A = theta)
    
    # Add random noise, so PO assumption isnt fulfilled 
    noise <- runif(n, min = 0.9, max = 1.1)
    
    # Apply noise and normalize
    p_noisy <- b * noise  
    b <- p_noisy/sum(p_noisy)
    
    list(a = a, b = b)
  }
  