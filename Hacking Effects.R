### Sample Size Hacking in different scenarios
sim <- simulation(p_C,p_E,n_pilot=1000,niter=10000, r=1)
mean(sim$actual_power_nmin) # Mean actual power of minimal sample
hist(sim$actual_power_nmin) # Histogram of power of minimal sample
df_sim <- data.frame("power_nmin" = sim$actual_power_nmin, "method" = sim$method)
ggplot(df_sim, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(sim$method == "PO"))/10000
length(which(sim$method == "ttest"))/10000
length(which(sim$method == "AfS"))/10000

# choose smallest sample size and look at method and power of it
colnames(sim$n_needed)[which(as.matrix(sim$n_needed) == min(sim$n_needed), arr.ind = TRUE)[1,2]] # ttest liefert kleinste sample size
sim$actual_power[which(as.matrix(sim$n_needed) == min(sim$n_needed), arr.ind = TRUE)] # min actual power: 0.398 & 0.404

# choose biggest sample size and look at method and power of it
colnames(sim$n_needed)[which(as.matrix(sim$n_needed) == max(sim$n_needed), arr.ind = TRUE)[1,2]] # ttest liefert größte sample size
sim$actual_power[which(as.matrix(sim$n_needed) == max(sim$n_needed), arr.ind = TRUE)] # min actual power: 0.999



# function to show all histograms
sim_plots <- function(data, group, plots=0, same_scale=TRUE){
  df_data <- data.frame()
  for (i in seq_along(1:length(data))) {
    df <- data.frame("power_nmin" = data[[i]]$actual_power_nmin, 
                     "method" = data[[i]]$method,
                     "group" = group[i])
    df_data <- rbind(df_data,df)
  }
  
  df_data <- df_data %>% group_by(group) %>%  mutate(mean = mean(power_nmin))
  df_data$method <- as.factor(df_data$method)
  df_data$method <- relevel(df_data$method, "ttest")
  
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

sim_ss_plots <- function(data, group, xlabel="Varied Component", plots=0){
  df_samp <- data.frame()
  for (i in seq_along(1:length(data))) {
    df <- data.frame("nmin" = data[[i]]$nmin, 
                     "group" = group[i])
    df_samp <- rbind(df_samp,df)
  }
  plot1 <- ggplot(df_samp) +
                   geom_pointrange(mapping = aes(x=group,y=nmin),
                                   stat = "summary",
                                   fun = mean, size=0.3)+
                   xlab(xlabel)+
                   ylab("Mean Minimum Sample Size")+
                   theme_bw()
  plot2 <- ggplot(df_samp) +
                   geom_pointrange(mapping = aes(x=group,y=nmin),
                                   stat = "summary",
                                   fun.min = function(z) {quantile(z,0.25)},
                                   fun.max = function(z) {quantile(z,0.75)},
                                   fun = median, size=0.3)+
                   xlab(xlabel)+
                   ylab("Median Minimum Sample Size with Quartiles")+
                   theme_bw()
  
  if(plots==0){
    grid.arrange(plot1, plot2, ncol=2)
  } else if(plots==1){
    return(plot1)
  } else if(plots == 2){
    return(plot2)
  }
  
}

sim_scatter_plot <- function(data, group, xlabel="Varied Component"){
  mean_col <- c()
  for(i in seq(1:length(data))){
    mean_col[i] <- mean(data[[i]]$actual_power_nmin)
  }
  sim_df <- data.frame(xcol =group, "nmin_power"=mean_col)
  ggplot(sim_df, aes(x=xcol, y=mean_col)) +
    geom_point()+
    geom_hline(yintercept=0.8, color = "red", linetype = "dotted")+
    xlab(xlabel)+
    ylab("Mean Actual Power")+
    theme_bw()
}

# Data frame with mean and standard deviation of the actual power distribution and the proportions of methods
summary_statistics <- function(data, group, d_label="Varied Component"){
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
  
  # calculate mean actual power
  for (i in seq_along(1:length(data))) {
    result[i, 2] <- mean(data[[i]]$actual_power_nmin)
  }
  
  # calculate standard deviation of actual power
  for (i in seq_along(1:length(data))) {
    result[i, 3] <- sd(data[[i]]$actual_power_nmin)
  }
  # calculate the proportion of methods
  for (i in seq_along(1:length(data))) {
    result[i,4] <- length(which(data[[i]]$method=="ttest"))/(length(data[[i]]$nmin)*0.01)
    result[i,5] <- length(which(data[[i]]$method=="PO"))/(length(data[[i]]$nmin)*0.01)
    result[i,6] <- length(which(data[[i]]$method=="WMW"))/(length(data[[i]]$nmin)*0.01)
  }
  # calculate mean minimum sample size
  for (i in seq_along(1:length(data))) {
    result[i,7] <- mean(data[[i]]$nmin)
  }
  # calculate median minimum sample size
  for (i in seq_along(1:length(data))) {
    result[i,8] <- median(data[[i]]$nmin)
  }
  # calculate IQR
  for (i in seq_along(1:length(data))) {
    result[i,9] <-   quantile(data[[i]]$nmin, c(0.25))
    result[i,10] <-   quantile(data[[i]]$nmin, c(0.75))
  }

  return(result)
}


# function to calculate median difference of smallest vs biggest sample size for every iteration
median_diff <- function(data, group, power = TRUE){
  diff_df <- data.frame( "sim" = c(1:length(data)),
                         "median sample size diff"=c(1:length(data)),
                         "mean sample size diff" = c(1:length(data)),
                         "median actual power diff"=c(1:length(data)),
                         "mean actual power diff"=c(1:length(data))
                         )
  for (i in seq_along(1:length(data))) {
    diff_vec_samp <- c()
    for(j in seq_along(1:length(data[[1]]$nmin))){
     diff_vec_samp[j] <-max(data[[i]]$n_needed[j,]) - min(data[[i]]$n_needed[j,])
    }
    diff_vec_power <- c()
    for(k in seq_along(1:length(data[[1]]$nmin))){
      diff_vec_power[k] <-max(data[[i]]$actual_power[k,]) - min(data[[i]]$actual_power[k,])
    }
    diff_df[i,1] <- group[i]
    diff_df[i,2] <- median(diff_vec_samp)
    diff_df[i,3] <- mean(diff_vec_samp[is.finite(diff_vec_samp)])
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

### different r #####
vary_r <- function(p_C, p_E, n_pilot = 10000, niter = 10000, r_vec){
  results <- vector(mode = "list", length = length(r_vec))
  for (i in seq_along(r_vec)) {
    sim <- simulation(p_C, p_E, n_pilot=n_pilot, niter=niter, r=r_vec[i])
    results[[i]] <- sim
  }
  return(results)
}

# thesis:
sim_r_1000 <- vary_r(p_C, p_E, n_pilot = 1000, r_vec=c(seq(0.1, 0.9, 0.2), 1, seq(2,7,1)))
sim_plots(sim_r_1000, c(seq(0.1, 0.9, 0.2), 1, seq(2,7,1)), plots=2)


# plot the mean minimum sample size and the median min sample size with quartiles
sim_ss_plots(data=sim_r_1000, group = c(seq(0.1, 0.9, 0.2), 1, seq(2,7,1)), 
             xlabel="Allocation Ratio", plots=0)

# summary
sum_r <- summary_statistics(sim_r_1000, group = c(seq(0.1, 0.9, 0.2), 1, seq(2,7,1)), d_label = "r")
print(xtable(sum_r, caption=1, digits=4), include.rownames = FALSE)

# scatter plot of mean actual powers
sim_scatter_plot(data=sim_r_1000, group = c(seq(0.1, 0.9, 0.2), 1, seq(2,7,1)), 
             xlabel="Allocation Ratio")

# how big is the median difference between smallest and biggest sample size in every iteration
diff_table_r <- median_diff(sim_r_1000, c(seq(0.1, 0.9, 0.2), 1, seq(2,7,1)), power = FALSE)
xtable(diff_table_r, caption=1)

# look at it for a smaller treatment effect
p_C_theta <- p_C
p_E_theta <- calc_p_E(p_C, theta_A = log(1.2))

sim_r_small_eff <- vary_r(p_C_theta, p_E_theta, n_pilot = 1000)
sim_plots(sim_r_small_eff, c(0.1, 0.3, 0.5, 0.7, 0.9, 1))


### vary p_C and p_E #####
## vary effect size theta
vary_p_theta <- function(n_pilot = 1000, niter = 10000, r_val = 1, cat = 6, theta_vec=log(seq(0.1, 2.5, 0.4))){
  p_results <- data.frame()
  results <- vector(mode = "list", length = length(theta_vec))
  for (i in seq_along(theta_vec)) {
    set.seed(1)
    p_C1 <- uniform_simplex(cat)[[1]] ## try random_simplex here
    p_E1 <- calc_p_E(p_C1, theta_A = theta_vec[i])
    noise <- runif(cat, min = 0.95, max = 1.05) # Add random noise, so PO assumption isnt fulfilled and PO method doesnt have advantage
    p_noisy <- p_E1 * noise  # Apply noise and normalize
    p_E1 <- p_noisy/sum(p_noisy)
    # Apply the simulation to the new probability vectors
    sim <- simulation(p_C=p_C1, p_E=p_E1, n_pilot=n_pilot, niter=niter, r=r_val)
    results[[i]] <- sim
  }
  return(results)
}

# which odds ratio does it then have for OR=1 after applying noise
set.seed(1)
p_C1 <- uniform_simplex(4)[[1]] ## try random_simplex here
p_E1 <- calc_p_E(p_C1, theta_A = log(1))
set.seed(1)
noise <- runif(4, min = 0.95, max = 1.05) # Add random noise, so PO assumption isnt fulfilled and PO method doesnt have advantage
p_noisy <- p_E1 * noise  # Apply noise and normalize
p_E1 <- p_noisy/sum(p_noisy)
exp(calculate_theta_A(p_C1, p_E1)[1,1])



sim_p_theta <- vary_p_theta(n_pilot = 1000, theta_vec = log(c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4)), cat=4)
sim_plots(sim_p_theta, c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4), plots = 2)
sim_ss_plots(sim_p_theta, c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4), "Odds Ratio", plots = 2)
sim_theta_diff <- median_diff(sim_p_theta, c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4), power=FALSE)
xtable(sim_theta_diff, caption=1)

sum_theta <- summary_statistics(sim_p_theta, group = c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4), d_label = "theta")
print(xtable(sum_theta, caption=1, digits = 4), include.rownames = FALSE)

# scatter plot of mean actual powers
sim_scatter_plot(data=sim_p_theta, group = c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4), 
                 xlabel="Treatment Effect")


# plots without outliers (not used in thesis)
  sim_p_theta_no_outlier <- sim_p_theta
  sim_p_theta_no_outlier[[5]] <- NULL
  sim_ss_plots(sim_p_theta_no_outlier, c(seq(0.2, 0.8, 0.2),1.25, 1.5, 2.5, 4), "Odds Ratio", plots = 2)
  
  
  
### example in Section 5.3.2 (not used)
sim_p_theta_exp <- vary_p_theta(n_pilot = 1000, theta_vec = c(0.4, 0.5, 0.6), cat=4)
median(sim_p_theta_exp[[1]]$n_needed)
median(sim_p_theta_exp[[2]]$n_needed)
median(sim_p_theta_exp[[3]]$n_needed)


p_C_exp1 <- uniform_simplex(4)[[1]]
p_E_exp1 <- calc_p_E(p_C_exp1, theta_A = 0.2)
p_C_exp2 <- uniform_simplex(4)[[1]]
p_E_exp2 <- calc_p_E(p_C_exp2, theta_A = 0.15)
p_C_exp3 <- uniform_simplex(4)[[1]] 
p_E_exp3 <- calc_p_E(p_C_exp3, theta_A = 0.25)
samplesize_ttestord(p_C_exp1, p_E_exp1, alpha=0.05, beta=0.2, r=1)$n_total
samplesize_ttestord(p_C_exp2, p_E_exp2, alpha=0.05, beta=0.2, r=1)$n_total
samplesize_ttestord(p_C_exp3, p_E_exp3, alpha=0.05, beta=0.2, r=1)$n_total


### now in other setting
#### big pilot study
sim_p_theta_npilot <- vary_p_theta(n_pilot = 10000, theta_vec = log(c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4)), cat=4, r=1)
sim_plots(sim_p_theta_npilot, c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4), plots = 0)

#### small pilot study
sim_p_theta_npilot2 <- vary_p_theta(n_pilot = 100, theta_vec = log(c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4)), cat=4, r=1)
sim_plots(sim_p_theta_npilot2, c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4), plots = 0)

#### unbalanced groups
  # r=4
sim_p_theta_r <- vary_p_theta(n_pilot = 1000, theta_vec = log(c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4)), cat=4, r_val=4)
sim_plots(sim_p_theta_r, c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4), plots = 0)
  # r=0.5
sim_p_theta_r2 <- vary_p_theta(n_pilot = 1000, theta_vec = log(c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4)), cat=4, r_val =0.5)
sim_plots(sim_p_theta_r2, c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4), plots = 0)
  # r=0.1
sim_p_theta_r3 <- vary_p_theta(n_pilot = 1000, theta_vec = log(c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4)), cat=4, r_val=0.1)
sim_plots(sim_p_theta_r3, c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4), plots = 0)

# r=7
sim_p_theta_r4 <- vary_p_theta(n_pilot = 1000, theta_vec = log(c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4)), cat=4, r_val = 7)
sim_plots(sim_p_theta_r4, c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4), plots = 0)

#### more categories
sim_p_theta_cat <- vary_p_theta(n_pilot = 1000, theta_vec = log(c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4)), cat=10, r=1)
sim_plots(sim_p_theta_cat, c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4), plots = 0)

#### less categories
sim_p_theta_cat2 <- vary_p_theta(n_pilot = 1000, theta_vec = log(c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4)), cat=3, r=1)
sim_plots(sim_p_theta_cat2, c(seq(0.2, 1, 0.2),1.25, 1.5, 2.5, 4), plots = 0)

# without noise 
sim_p_theta2 <- vary_p_theta(n_pilot = 1000, theta_vec = log(seq(0.2, 1.8, 0.2)), cat=4)
sim_plots(sim_p_theta2, seq(0.2, 1.8, 0.2), plots = 2)
sim_ss_plots(sim_p_theta2, seq(0.2, 1.8, 0.2), "Odds Ratio", plots = 2)
sim_theta_diff2 <- median_diff(sim_p_theta2, seq(0.2, 1.8, 0.2), power=FALSE)
xtable(sim_theta_diff2, caption=1)


# make this with more odds ratio values, to make a nice scatter plot
sim_p_scatter <- vary_p_theta(n_pilot = 1000, cat=4, theta_vec = log(c(seq(0.2, 1, 0.1), 1.25, 1.5, 1.75, seq(2,5,0.5))))

# make a scatter plot with the mean actual Power
mean_min_power <- c()
for(i in seq(1:length(sim_p_scatter))){
  mean_min_power[i] <- mean(sim_p_scatter[[i]]$actual_power_nmin)
}

sim_p_df <- data.frame("p"=c(seq(0.2, 1, 0.1), 1.25, 1.5, 1.75, seq(2,5,0.5)), "nmin_power"=mean_min_power)
ggplot(sim_p_df, aes(x=p, y=nmin_power)) +
  geom_point()+
  geom_hline(yintercept=0.8, color = "red", linetype = "dotted")+
  xlab("Odds Ratio")+
  ylab("Mean Actual Power")+
  theme_bw()

# just repeat different p with the same bias and different seeds
vary_p <- function(n_pilot = 1000, niter = 10000, r = 1, cat = 6, rep=10, theta = log(1.8)){
  vec <- c(1:rep)
  p_results <- data.frame()
  results <- vector(mode = "list", length = rep)
  for (i in seq_along(vec)) {
    set.seed(i)
    p_C1 <- random_simplex(cat)
    p_E1 <- calc_p_E(p_C1, theta_A = theta)
    sim <- simulation(p_C1, p_E1, n_pilot=n_pilot, niter=niter, r=r)
    results[[i]] <- sim
  }
  return(results)
}
sim_p <- vary_p(n_pilot = 10000, rep=5)

mean_sim_p <- c()
for (i in seq_along(1:length(sim_p))) {
  mean_sim_p[[i]] <- mean(sim_p[[i]]$actual_power_nmin)
}
mean_sim_p

sim_plots(sim_p, c(1:length(sim_p)))

##### normally distributed categories #####
# create vector that is normally distributed around a category (4)
set.seed(1)
p_C_norm <- normal_simplex(6,4)[[1]]
p_E_norm <- normal_simplex(6,4)[[2]]

sim_n <- simulation(p_C_norm, p_E_norm, n_pilot=10000, niter=10000)
mean(sim_n$actual_power_nmin)
hist(sim_n$actual_power_nmin)

sim_n_power <-  as.data.frame(sim_n$actual_power) %>%
  pivot_longer(everything(), names_to = "method", values_to = "power")
ggplot(sim_n_power, aes(x = power, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity") +
  geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
  facet_wrap(method ~.) 

ggplot(sim_n_power, aes(x = power, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity") +
  geom_vline(xintercept = 0.8, color = "red", linetype="dotted")

# Plot of minimal power per Method
df_sim_n <- data.frame("power_nmin" = sim_n$actual_power_nmin, "method" = sim_n$method)

ggplot(df_sim_n, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")+
  geom_vline(aes(xintercept = mean(power_nmin)), color = "red")+
  geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
  theme_minimal()
length(which(sim_n$method == "PO"))/10000
length(which(sim_n$method == "ttest"))/10000
length(which(sim_n$method == "AfS"))/10000
sim_n$actual_power[which(as.matrix(sim_n$n_needed) == min(sim_n$n_needed), arr.ind = TRUE)] # min actual power: 0.403

### function that does this a few times
vary_norm <- function(n_pilot = 1000, niter = 10000, r = 1, rep = 10, cat=2, vec_length = 6, theta = log(1.8)){
  rep_vec <- c(1:rep)
  results <- vector(mode = "list", length = rep)
  for (i in seq_along(1:rep)) {
    set.seed(i)
    p_C_norm <- normal_simplex_theta(vec_length, cat, theta)[[1]]
    p_E_norm <- normal_simplex_theta(vec_length, cat, theta)[[2]]
    sim <- simulation(p_C_norm, p_E_norm, r = r, niter = niter, n_pilot = n_pilot)
    results[[i]] <- sim
  }
  return(results)
}

sim_norm <- vary_norm(cat=4, n_pilot = 10000, theta = log(2), rep=1)
sim_plots(sim_norm, c(1:length(sim_norm)))

mean_sim_norm <- c()
for (i in seq_along(1:length(sim_norm))) {
  mean_sim_norm[[i]] <- mean(sim_norm[[i]]$actual_power_nmin)
}
mean_sim_norm

# function with different extreme point locations
vary_norm_loc <- function(n_pilot = 1000, niter = 10000, r = 1, vec_length = 6, seed = 1, theta = log(1.8)){
  loc <- c(1:vec_length)
  results <- vector(mode = "list", length = vec_length)
  for (i in seq_along(loc)) {
    set.seed(seed)
    p_C_norm <- normal_simplex_theta(vec_length, i, theta)[[1]]
    p_E_norm <- normal_simplex_theta(vec_length, i, theta)[[2]]
    sim <- simulation(p_C_norm, p_E_norm, r = r, niter = niter, n_pilot = n_pilot)
    results[[i]] <- sim
  }
  return(results)
}

sim_norm_loc <- vary_norm_loc(n_pilot = 1000, theta = log(2))
sim_plots(sim_norm_loc, c(1:length(sim_norm_loc)))

sim_norm_loc2 <- vary_norm_loc(n_pilot = 1000, theta = log(1.2), seed = 1)
sim_plots(sim_norm_loc2, c(1:length(sim_norm_loc2)))


#### Try Boulesteix idea of sampling the p vector with a continuous variable (in Thesis)

vary_norm_2 <- function(n_pilot = 1000, niter = 10000, r = 1, rep = 10, vec_length = 6, theta = NA){
  rep_vec <- c(1:rep)
  results <- vector(mode = "list", length = rep)
  p_C_norm <- normal_vector(cat = vec_length, niter = niter, seed = 1, theta_A = theta)[[1]]
  p_E_norm <- normal_vector(cat = vec_length, niter = niter, seed = 1, theta_A = theta)[[2]]
  for (i in seq_along(1:rep)) {
    set.seed(i)
    sim <- simulation(p_C_norm, p_E_norm, r = r, niter = niter, n_pilot = n_pilot)
    results[[i]] <- sim
  }
  return(results)
}
# random noise
sim_norm_new <- vary_norm_2(n_pilot = 10000, rep=5)
sim_plots(sim_norm_new, c(1:length(sim_norm_new)))


# set effect size
# use same odds ratio as p_C and p_E have, to make it comparable to the other simulations
exp(calculate_theta_A(p_C,p_E)[1,1])
plot(normal_vector(cat = 4, niter = 10000, seed = 1, theta_A = log(0.4821251))[[1]])
sim_norm_new2 <- vary_norm_2(n_pilot = 1000, rep=1, theta = log(0.4821251), vec_length = 4)
sim_plots(sim_norm_new2, c(""), plots=2)

sum_norm <- summary_statistics(sim_norm_new2, group =  c("norm"), d_label = "norm")
print(xtable(sum_norm, caption=1, digits = 4), include.rownames = FALSE)


median_diff(sim_norm_new2, c(""), power = FALSE)

### for comparison the not-normally distributed simulation
sim_npilot_1000 <- list(sim_npilot[[5]])
sim_plots(sim_npilot_1000, c(""), plots=2)
sum_norm2 <- summary_statistics(sim_npilot_1000, group =  c("not-norm"), d_label = "not-norm")
print(xtable(sum_norm2, caption=1, digits = 4), include.rownames = FALSE)


mean(sim_npilot_1000[[1]]$nmin)
median(sim_npilot_1000[[1]]$nmin)
quantile(sim_npilot_1000[[1]]$nmin, c(0.25,0.75))
median_diff(sim_npilot_1000, c(""), power = FALSE)


# with more categories (no change)
sim_norm_new3 <- vary_norm_2(n_pilot = 1000, rep=1, theta = log(0.4821251), vec_length = 10)
sim_plots(sim_norm_new3, c(1:length(sim_norm_new3)))

# with 3 categories (a little more overpowering)
sim_norm_new4 <- vary_norm_2(n_pilot = 1000, rep=1, theta = log(0.4821251), vec_length = 3)
sim_plots(sim_norm_new4, c(1:length(sim_norm_new4)))

# large n_pilot  (no change)
sim_norm_npilot1 <- vary_norm_2(n_pilot = 10000, rep=1, theta = log(0.4821251), vec_length = 4)
sim_plots(sim_norm_npilot1, c(1:length(sim_norm_npilot1)))

# small n_pilot (no change)
sim_norm_npilot2 <- vary_norm_2(n_pilot = 100, rep=1, theta = log(0.4821251), vec_length = 4)
sim_plots(sim_norm_npilot2, c(1:length(sim_norm_npilot2)))

# strong treatment effect (no change)
sim_norm_theta1 <- vary_norm_2(n_pilot = 1000, rep=1, theta = log(2.5), vec_length = 4)
sim_plots(sim_norm_theta1, c(1:length(sim_norm_theta1)))

# weak treatment effect (no change)
sim_norm_theta2 <- vary_norm_2(n_pilot = 1000, rep=1, theta = log(0.8), vec_length = 4)
sim_plots(sim_norm_theta2, c(1:length(sim_norm_theta2)))

### different number of categories #####
#### try it with fixed effect size for p_E
vary_prob_length_PO <- function(n_pilot = 1000, niter = 10000, r = 1, seed = 1, theta=log(1.8)){
  vec_length<- c(3:14)
  results <- vector(mode = "list", length = length(vec_length))
  for (i in seq_along(vec_length)) {
    set.seed(seed)
    p_C_u <- uniform_simplex(vec_length[i])[[1]]
    p_E_u <- calc_p_E(p_C_u, theta_A = theta)
#    noise <- runif(vec_length[i], min = 0.95, max = 1.05) 
#    p_noisy <- p_E_u * noise  
#    p_E_u <- p_noisy/sum(p_noisy)
    sim_u <- simulation(p_C_u, p_E_u, n_pilot=n_pilot, niter=niter, r=r)
    results[[i]] <- sim_u
  }
  return(results)
}

sim_cat_PO <- vary_prob_length_PO(n_pilot = 1000)
sim_plots(sim_cat_PO, c(3:14))

sim_cat_PO2 <- vary_prob_length_PO(n_pilot = 1000, theta = log(1.3))
sim_plots(sim_cat_PO2, c(3:14))

sim_cat_PO3 <- vary_prob_length_PO(n_pilot = 100, theta = log(1.8))
sim_plots(sim_cat_PO3, c(3:14))

# same odds ratio as p_C and p_E
exp(calculate_theta_A(p_C, p_E)[1,1])
sim_cat_PO_same <- vary_prob_length_PO(n_pilot = 1000, theta = log(0.4821251))
sim_plots(sim_cat_PO_same, c(3:14), plots=2)

sum_cat <- summary_statistics(sim_cat_PO_same, group = c(3:14), d_label = "categories")
print(xtable(sum_cat, caption=1, digits=4), include.rownames = FALSE)

# scatter plot of mean actual powers
sim_scatter_plot(data=sim_cat_PO_same, group = c(3:14), 
                 xlabel="Number of Categories")


## visualise the mean minimum sample size
sim_ss_plots(sim_cat_PO_same, c(3:14), xlabel = "Number of Categories", plots = 2)

# how big is the median difference between smallest and biggest sample size in every iteration
diff_table_cat <- median_diff(sim_cat_PO_same, c(3:14), power = FALSE)
xtable(diff_table_cat, caption=1)


# same odds ratio as p_C and p_E but with n_pilot=200
sim_cat_PO_same_small <- vary_prob_length_PO(n_pilot = 200, theta = log(0.4821251))
sim_cat_PO_same_small2 <- vary_prob_length_PO(n_pilot = 100, theta = log(0.4821251))
sim_plots(sim_cat_PO_same_small, c(3:14))



#### try it with ad-hoc method
vary_prob_length_ad_hoc <- function(n_pilot = 1000, niter = 10000, r = 1, seed = 1){
  vec_length<- c(3:15)
  results <- vector(mode = "list", length = length(vec_length))
  for (i in seq_along(vec_length)) {
    set.seed(seed)
    p_C_u <- uniform_simplex(vec_length[i])[[1]]
    shift_vec <- p_C_u * 0 
    shift_vec[1] <- -p_C_u[1]*0.2
    shift_vec[length(p_C_u)] <- p_C_u[length(p_C_u)]*0.2
    p_E_u <- p_C_u + shift_vec
    sim_u <- simulation(p_C_u, p_E_u, n_pilot=n_pilot, niter=niter, r=r)
    results[[i]] <- sim_u
  }
  return(results)
}

sim_cat_ad_hoc <- vary_prob_length_ad_hoc()
sim_plots(sim_cat_ad_hoc, c(3:15))

### different n_pilot #####
vary_npilot <- function(p_C, p_E, niter = 10000, r = 1){
  npilot_vec <- c(50, 100, 200, 500, 1000, 2500, 5000, 10000, 15000)
  results <- vector(mode = "list", length = length(npilot_vec))
  for (i in seq_along(npilot_vec)) {
    sim <- simulation(p_C, p_E, r = r, niter = niter, n_pilot = npilot_vec[i])
    results[[i]] <- sim
  }
  return(results)
}

sim_npilot <- vary_npilot(p_C, p_E)
sim_npilot[[10]]<- NULL

sim_plots(sim_npilot, c(50, 100, 200, 500, 1000, 2500, 5000, 10000, 15000), plots=2)
  
sum_npilot <- summary_statistics(sim_npilot, group =  c(50, 100, 200, 500, 1000, 2500, 5000, 10000, 15000), d_label = "n_pilot")
print(xtable(sum_npilot, caption=1, digits=4), include.rownames = FALSE)

# scatter plot of mean actual powers
sim_scatter_plot(data=sim_npilot, group = c(50, 100, 200, 500, 1000, 2500, 5000, 10000, 15000), 
                 xlabel="Pilot Study Sample Size")

# plot the mean minimum sample size and the median min sample size with quartiles
sim_ss_plots(data=sim_npilot, 
             group = c(50, 100, 200, 500, 1000, 2500, 5000, 10000, 15000),
             xlabel="Pilot Study Sample Size", plots=2)

# how big is the median difference between smallest and biggest sample size in every iteration
diff_table_npilot <- median_diff(sim_npilot, c(50, 100, 200, 500, 1000, 2500, 5000, 10000, 15000))
print(xtable(diff_table_npilot), include.rownames = FALSE)

### now in other setting
#### unbalanced groups
sim_p_npilot_r <- vary_npilot(p_C, p_E, r=4)
sim_plots(sim_p_npilot_r, c(50, 100, 200, 500, 1000, 2500, 5000, 10000, 15000), plots=0)

sim_p_npilot_r2 <- vary_npilot(p_C, p_E, r=0.5)
sim_plots(sim_p_npilot_r2, c(50, 100, 200, 500, 1000, 2500, 5000, 10000, 15000), plots=0)


### different iterations #####

### write function that does all this
vary_niter <- function(p_C, p_E, n_pilot = 1000, r = 1){
  niter_vec <- c(100, 500, 1000, 2500, 5000, 7500, 10000, 15000, 20000)
  results <- vector(mode = "list", length = length(niter_vec))
  for (i in seq_along(niter_vec)) {
    sim <- simulation(p_C, p_E, r = r, n_pilot = n_pilot, niter = niter_vec[i])
    results[[i]] <- sim
  }
  return(results)
}
sim_niter <- vary_niter(p_C, p_E, n_pilot = 1000)
sim_plots(sim_niter, c(100, 500, 1000, 2500, 5000, 7500, 10000, 15000, 20000), same_scale = FALSE)

# Summary statistics
sum_iter <- summary_statistics(sim_niter, group = c(100, 500, 1000, 2500, 5000, 7500, 10000, 15000, 20000), d_label = "iterations")
print(xtable(sum_iter, caption=1, digits=4), include.rownames = FALSE)


### vary alpha ########
vary_alpha <- function(p_C, p_E, n_pilot = 10000, niter = 10000, r=1){
  alpha_vec <- c(0.01, 0.025, 0.05, 0.075, 0.1)
  results <- vector(mode = "list", length = length(alpha_vec))
  for (i in seq_along(alpha_vec)) {
    sim <- simulation(p_C, p_E, n_pilot=n_pilot, niter=niter, alpha=alpha_vec[i], r=r)
    results[[i]] <- sim
  }
  return(results)
}

sim_alpha <- vary_alpha(p_C, p_E, n_pilot = 1000)
sim_plots(sim_alpha, c(0.01, 0.025, 0.05, 0.075, 0.1))

### vary beta ######
vary_beta <- function(p_C, p_E, n_pilot = 10000, niter = 10000, r=1){
  beta_vec <- c(0.7, 0.75, 0.8, 0.85, 0.9, 0.95)
  results <- vector(mode = "list", length = length(beta_vec))
  for (i in seq_along(beta_vec)) {
    sim <- simulation(p_C, p_E, n_pilot=n_pilot, niter=niter, beta=beta_vec[i], r=r)
    results[[i]] <- sim
  }
  return(results)
}

sim_beta <- vary_beta(p_C, p_E, n_pilot = 1000)
sim_plots(sim_beta, c(0.7, 0.75, 0.8, 0.85, 0.9, 0.95))


