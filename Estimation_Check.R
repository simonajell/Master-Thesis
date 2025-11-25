
# First run the file "Functions.R" to load the packages and functions needed in this File

#### Estimation check for the t-test ####
### Function to compare power calculation from this t-test implementation and R base t-test
  comp_ttest <- function(alpha = 0.05, beta = 0.2, r = 1, iter = 1000){
    results_R <-  c()
    results_me <-  c()
    for (i in seq_along(1:iter)) {
      set.seed(i)
      # Sample a number of categories
      prob_length <- sample(c(3:14), 1)
      # Sample a treatment effect
      theta <- runif(1, min = 0.1, max = 5)
      
      print(i)
      # Generate probability vectors
      p <- generate_two_simplex_vectors(prob_length, log(theta))
      p_C <- p[[1]]
      p_E <- p[[2]]
      # Apply t-test function
      ttest <- samplesize_ttestord(p_C=p_C, p_E=p_E, alpha, beta, r)
      results_me <- rbind(results_me, ttest$actual_power)
      # save the difference of means
      del_sig <- calculate_delta_A_sigma(p_C, p_E)
      results_R <- c(results_R, 
                     # calculate the t-test with the base r function
                     power.t.test(n=ttest$n_E, delta = del_sig$delta_A, power = NULL,
                                  sd = del_sig$sigma, alternative = "two.sided", type = "two.sample")$power)
    }
    return(data.frame("ttestord"=results_me, "R"=results_R)) 
  }
  
  # Apply comparison function
  sim_comp_ttest <- comp_ttest(iter = 10000)
  
  # Mean absolute difference between power calculated from this t-test and R base t-test
  mean(abs(sim_comp_ttest$ttestord-sim_comp_ttest$R)) # 0.0002397253
  
  # Plot comparison
  ggplot(data=sim_comp_ttest, aes(x=ttestord, y=R))+
    geom_point(alpha=0.6)+
    geom_abline(slope=1, color="red")+
    xlab("t-test")+
    ylab("base R")+
    theme_bw()

### Function to compare sample size calculation from this t-test implementation and R base t-test
  comp_ttest_samp <- function(alpha = 0.05, beta = 0.2, r = 1, iter = 1000, theta){
    results_R <-  c()
    results_me <-  c()
    for (i in seq_along(1:iter)) {
      set.seed(i)
      # Sample a number of categories
      prob_length <- sample(c(3:14), 1)
      # Sample a treatment effect
      theta <-  runif(1, min = 0.1, max = 5)
      
      print(i)
      # Generate probability vectors
      p <- generate_two_simplex_vectors(prob_length, log(theta))
      p_C <- p[[1]]
      p_E <- p[[2]]
      # Apply t-test function
      ttest <- samplesize_ttestord(p_C=p_C, p_E=p_E, alpha, beta, r)
      results_me <- rbind(results_me, ttest$n_total)
      # save the difference of means
      del_sig <- calculate_delta_A_sigma(p_C, p_E)
      results_R <- c(results_R, 
                     # Calculate the sample size with the base R function
                     ceiling(power.t.test(n=NULL, delta = del_sig$delta_A, power = 0.8, sig.level = 0.05,
                                          sd = del_sig$sigma, alternative = "two.sided", type = "two.sample")$n)*2)
    }
    return(data.frame("ttestord"=results_me, "R"=results_R)) 
  }
  
  # Apply comparison function
  sim_comp_ttest_samp <- comp_ttest_samp(iter=10000)
  
  # Mean absolute difference between sample size calculated from this t-test and R base t-test
  mean(abs(sim_comp_ttest_samp$ttestord-sim_comp_ttest_samp$R)) # 0.0386
  
  # Percentage of unequal sample size calculations 
  length(which(sim_comp_ttest_samp$ttestord != sim_comp_ttest_samp$R, arr.ind = TRUE))/nrow(sim_comp_ttest_samp)
  
  # Plot comparison  
  ggplot(data=sim_comp_ttest_samp, aes(x=ttestord, y=R))+
    geom_point(alpha=0.6)+
    geom_abline(slope=1, color="red")+
    xlab("t-test")+
    ylab("base R")+
    ylim(0,100000)+
    xlim(0,100000)+
    theme_bw()
  
#### Estimation check for the PO method ####
### Function to compare theta estimation by White (2023) and true theta
  theta_comparison <- function(comp_iter, r){
    theta_results <- data.frame("theta_A" = c(NA), "est_theta"=c(NA), "diff" = c(NA))
    p_c <- 0
    p_e <- 0
    for (i in seq_along(1:comp_iter)) {
      set.seed(i)
      print(i)
      # Sample a number of categories
      vec_length <- sample(c(3:14), 1)
      # Sample a treatment effect
      theta_random <- runif(1, min = 0.1, max = 5)
      # Generate probability vectors (not generate_two_simplex_vectors because no noise should be applied)
      p_c <- random_simplex(n=vec_length)
      p_e <- calc_p_E(p_C = p_c, theta_A = log(theta_random))
      # Save true log odds ratio
      theta_results[i,1] <- log(theta_random)
      # Calculate log odds ratio White method
      theta_results[i,2] <- calculate_theta_A(p_c, p_e, r)[1,1]
    }
    # Calculate the absolute difference between the true and estimated log odds ratio
    theta_results$diff <- abs(theta_results$theta_A - theta_results$est_theta)
    return(theta_results)
  }
  
  # Apply comparison function
  sim_comp_PO_theta <- theta_comparison(comp_iter = 10000, r=1)
  
  # Mean absolute difference between estimated theta and true theta
  mean(sim_comp_PO_theta$diff) # 0.0005697382
  
  # Plot comparison
  ggplot(data=sim_comp_PO_theta, aes(x=theta_A, y=est_theta))+
    geom_point(alpha=0.6)+
    geom_abline(slope=1, color="red")+
    xlab("True Log Odds Ratio")+
    ylab("Estimated Theta")+
    theme_bw()

### Function to compare variance estimation by White (2023) and by Whitehead (1993) with formula from Kieser (2020)
  var_comparison <- function(comp_iter, r){
    var_results <- data.frame("theta_A" = c(NA), "var_w"=c(NA), "var_k" = c(NA), "diff" = c(NA))
    for (i in seq_along(1:comp_iter)) {
      set.seed(i)
      print(i)
      # Sample the number of categories
      vec_length <- sample(3:14, 1)
      # Sample the odds ratio
      theta_random <- runif(1, min = 0.1, max = 5)
      # Generate probability vectors (not generate_two_simplex_vectors because no noise should be applied)
      p_c_po <- random_simplex(vec_length)
      p_e_po <- calc_p_E(p_c_po, theta_A = log(theta_random))
      # Save the true log odds ratio
      var_results[i,1] <- theta_random
      # Calculate variance with White method
      var_results[i,2] <- calculate_theta_N(p_c_po, p_e_po, r)[1,2]^2
      # Calculate variance with Kieser method
      x = 0
      for (j in seq_along(1:vec_length)){
        x = x + ((r*p_c_po[j] + p_e_po[j]) / (1 + r))^3
      }
      var_results[i,3] <- ((1+r)^2/r)* (3/(1-x))
    }
    # Calculate difference between the two variances
    var_results$diff <- abs(var_results$var_w - var_results$var_k)
    return(var_results)
  }
  
  # Apply comparison function
  sim_comp_PO_var <- var_comparison(comp_iter = 10000, r=1)
  
  # Mean absolute difference between estimated variance by White and Whitehead
  mean(sim_comp_PO_var$diff) #0.000818111

  # Plot comparison
  ggplot(data=sim_comp_PO_var, aes(x=var_w, y=var_k))+
    geom_point(alpha=0.6)+
    geom_abline(slope=1, color="red")+
    xlab("Variance White")+
    ylab("Variance Whitehead")+
    theme_bw()
  
  
### Function to compare sample size calculation by White (2023) and by Whitehead (1993) with formula from Kieser (2020)
  comp_PO_samp <- function(alpha=0.05, beta=0.2, r=1, iter=1000){
    results_NN <-  data.frame()
    results_AA <-  data.frame()
    results_NA <-  data.frame()
    results_k <- data.frame()
    results_k_known <- data.frame()
    theta_vec <- rep(0,iter)
    true_theta_vec <- rep(0,iter)
    for (i in seq_along(1:iter)) {
      set.seed(i)
      # Sample the number of categories
      prob_length <- sample(3:14, 1)
      # Sample the treatment effect
      theta_A <- runif(1, 0.1, 5)
      # Save the true odds ratio
      true_theta_vec[i] <- theta_A
      print(i)
      # Generate probability vectors (not generate_two_simplex_vectors because no noise should be applied)
      p_C <- random_simplex(prob_length)
      p_E <- calc_p_E(p_C, theta_A = log(theta_A))
      # Estimate treatment effect with White method and save
      theta_vec[i] <- calculate_theta_A(p_C, p_E)[1,1]
      # Calculate the sample sizes with the NN method
      results_NN <- rbind(results_NN, as.data.frame(samplesize_po_NN(p_C=p_C, p_E=p_E, alpha, beta, r)))
      # Calculate the sample size with the AA method
      results_AA <- rbind(results_AA, as.data.frame(samplesize_po_AA(p_C=p_C, p_E=p_E, alpha, beta, r)))
      # Calculate the sample size with the NA method
      results_NA <- rbind(results_NA, as.data.frame(samplesize_po_NA(p_C=p_C, p_E=p_E, alpha, beta, r)))
      # Calculate the sample size with the Whitehead method with formula from Kieser
      results_k <- rbind(results_k, as.data.frame(samplesize_po_kieser(p_C=p_C, p_E=p_E, alpha, beta, r)))
    }
    return(list(results_NN, results_k, results_AA, results_NA, results_k_known, theta_vec, true_theta_vec))
  }
  
  # Apply comparison function
  sim_comp_PO_samp <- comp_PO_samp(iter = 10000)
  
  # Percentage of equal sample sizes calculated from White and Whitehead
  length(which(sim_comp_PO_samp[[1]]$n_total == sim_comp_PO_samp[[2]]$n_total))/10000 # 96.63% , 149 unequal
  
  # Mean absolute difference between sample size calculated by White and Whitehead
  mean(abs(sim_comp_PO_samp[[1]]$n_total - sim_comp_PO_samp[[2]]$n_total)) # 16.3688
  
  # Plot comparison
  SS_comp_PO <- data.frame("White" = sim_comp_PO_samp[[1]]$n_total,
                           "Whitehead" = sim_comp_PO_samp[[2]]$n_total)
  ggplot(data=SS_comp_PO, aes(x=White, y=Whitehead))+
    geom_point(alpha=0.6)+
    geom_abline(slope=1, color="red")+
    xlab("Sample Size by White")+
    ylab("Sample Size by Whitehead")+
    ylim(0,100000)+
    xlim(0,100000)+
    theme_bw()

  