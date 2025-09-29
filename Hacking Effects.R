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

### different r #####
vary_r <- function(p_C, p_E, n_pilot = 10000, niter = 10000, r_vec=c(0.1, 0.3, 0.5, 0.7, 0.9, 1)){
  results <- vector(mode = "list", length = length(r_vec))
  for (i in seq_along(r_vec)) {
    sim <- simulation(p_C, p_E, n_pilot=n_pilot, niter=niter, r=r_vec[i])
    results[[i]] <- sim
  }
  return(results)
}

sim_r <- vary_r(p_C, p_E)

mean_sim_r <- c()
for (i in seq_along(1:length(sim_r))) {
  mean_sim_r[[i]] <- mean(sim_r[[i]]$actual_power_nmin)
}
mean_sim_r

sim_plots(sim_r, c(0.1, 0.3, 0.5, 0.7, 0.9, 1))

## for small n_pilot, to see what it looks like in a scenario with sample size hacking
sim_r_small <- vary_r(p_C, p_E, n_pilot = 200)
sim_plots(sim_r_small, c(0.1, 0.3, 0.5, 0.7, 0.9, 1))

sim_r_small_500 <- vary_r(p_C, p_E, n_pilot = 500)
sim_plots(sim_r_small_500, c(0.1, 0.3, 0.5, 0.7, 0.9, 1))

sim_r_1000 <- vary_r(p_C, p_E, n_pilot = 1000)
sim_plots(sim_r_1000, c(0.1, 0.3, 0.5, 0.7, 0.9, 1), plots = 2)

#Proportions:
## 0.1: PO=79.82%  T-Test=6.59%   WMW=13.59%
## 0.3: PO=91.67%  T-Test=5.37%   WMW=2.96%
## 0.5: PO=82.86%  T-Test=16.9%   WMW=0.24%
## 0.7: PO=66.37%  T-Test=33.63%  WMW=0
## 0.9: PO=52.34%  T-Test=47.66%  WMW=0
## 1:   PO=44.99%  T-Test=55.01%  WMW=0
# proportions of methods
length(which(sim_r_1000[[6]]$method == "PO"))/100
length(which(sim_r_1000[[6]]$method == "ttest"))/100
length(which(sim_r_1000[[6]]$method == "WMW"))/100


# look at it for a smaller treatment effect
p_C_theta <- p_C
p_E_theta <- calc_p_E(p_C, theta_A = log(1.2))

sim_r_small_eff <- vary_r(p_C_theta, p_E_theta, n_pilot = 1000)
sim_plots(sim_r_small_eff, c(0.1, 0.3, 0.5, 0.7, 0.9, 1))

### make a scatter plot
### more r values, to more points in the scatter plot
sim_r_scatter <- vary_r(p_C, p_E, n_pilot=1000,r_vec=seq(0.1,1,0.1))
sim_r_scatter2 <- vary_r(p_C, p_E, n_pilot=10000,r_vec=seq(0.1,1,0.1))

# make a scatter plot with the mean actual Power
mean_min_power_r <- c()
for(i in seq(1:length(sim_r_scatter))){
  mean_min_power_r[i] <- mean(sim_r_scatter[[i]]$actual_power_nmin)
}

sim_r_df <- data.frame("p"=seq(0.1,1,0.1), "nmin_power"=mean_min_power_r)
ggplot(sim_r_df, aes(x=p, y=nmin_power)) +
  geom_point()+
  geom_hline(yintercept=0.8, color = "red", linetype = "dotted")+
  xlab("Allocation Ratio")+
  ylab("Mean Actual Power")+
  theme_bw()
### vary p_C and p_E #####
## vary effect size theta
vary_p_theta <- function(n_pilot = 1000, niter = 10000, r = 1, cat = 6, theta_vec=log(seq(0.1, 2.5, 0.4))){
  p_results <- data.frame()
  results <- vector(mode = "list", length = length(theta_vec))
  for (i in seq_along(theta_vec)) {
    set.seed(1)
    p_C1 <- uniform_simplex(cat)[[1]] ## try random_simplex here
    p_E1 <- calc_p_E(p_C1, theta_A = theta_vec[i])
    noise <- runif(cat, min = 0.95, max = 1.05) # Add random noise, so PO assumption isnt fulfilled and PO method doesnt have advantage
    p_noisy <- p_E1 * noise  # Apply noise and normalize
    p_E1 <- p_noisy/sum(p_noisy)
    
    sim <- simulation(p_C1, p_E1, n_pilot=n_pilot, niter=niter, r=r)
    results[[i]] <- sim
  }
  return(results)
}

sim_p_theta_100 <- vary_p_theta(n_pilot = 100) 
sim_p_theta_1000 <- vary_p_theta(n_pilot = 1000, theta_vec = log(seq(0.2, 3.5, 0.4)))
sim_p_theta_10000 <- vary_p_theta(n_pilot = 10000)

sim_plots(sim_p_theta_100, paste0("log(",seq(0.1, 2.5, 0.4), ")", " = ", round(log(seq(0.1, 2.5, 0.4)), 2)))
sim_plots(sim_p_theta_1000, paste0("log(",seq(0.2, 3.5, 0.4), ")", " = ", round(log(seq(0.2, 3.5, 0.4)), 2)))
sim_plots(sim_p_theta_10000, paste0("log(",seq(0.1, 2.5, 0.4), ")", " = ", round(log(seq(0.1, 2.5, 0.4)), 2)))


# make this with more odds ratio values, to make a nice scatter plot
sim_p_scatter <- vary_p_theta(n_pilot = 1000, theta_vec = log(seq(0.1, 2.5, 0.1)))

# make a scatter plot with the mean actual Power
mean_min_power <- c()
for(i in seq(1:length(sim_p_scatter))){
  mean_min_power[i] <- mean(sim_p_scatter[[i]]$actual_power_nmin)
}
which(seq(0.1, 2.5, 0.1)==1, arr.ind = TRUE)
mean_min_power[8]

sim_p_df <- data.frame("p"=seq(0.1, 2.5, 0.1), "nmin_power"=mean_min_power)
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

sim_norm1.1 <- vary_norm(cat=4, n_pilot = 1000, theta = log(2), rep=1)
sim_plots(sim_norm1.1, c(1:length(sim_norm1.1)))

sim_norm2 <- vary_norm(cat=4, n_pilot = 10000, theta = log(1.2), rep=1)
sim_plots(sim_norm2, c(1:length(sim_norm2)))

sim_norm2.2 <- vary_norm(cat=4, n_pilot = 1000, theta = log(1.2), rep=1)
sim_plots(sim_norm2.2, c(1:length(sim_norm2.2)))

sim_norm3 <- vary_norm(cat=4, n_pilot = 10000, theta = log(1.05), rep=1)
sim_plots(sim_norm3, c(1:length(sim_norm3)))

sim_norm3.2 <- vary_norm(cat=4, n_pilot = 1000, theta = log(1.05), rep=1)
sim_plots(sim_norm3.2, c(1:length(sim_norm3.2)))

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


#### Try Boulesteix idea of sampling the p vector with a continuous variable

vary_norm_2 <- function(n_pilot = 1000, niter = 10000, r = 1, rep = 10, vec_length = 6, theta = NA){
  rep_vec <- c(1:rep)
  results <- vector(mode = "list", length = rep)
  for (i in seq_along(1:rep)) {
    set.seed(i)
    p_C_norm <- normal_vector(cat = vec_length, niter = niter, seed = i, theta_A = theta)[[1]]
    p_E_norm <- normal_vector(cat = vec_length, niter = niter, seed = i, theta_A = theta)[[2]]
    sim <- simulation(p_C_norm, p_E_norm, r = r, niter = niter, n_pilot = n_pilot)
    results[[i]] <- sim
  }
  return(results)
}
# random noise
sim_norm_new <- vary_norm_2(n_pilot = 10000, rep=5)
sim_plots(sim_norm_new, c(1:length(sim_norm_new)))

# set effect size
  sim_norm_new2 <- vary_norm_2(n_pilot = 1000, rep=1, theta = log(1.8))
  sim_plots(sim_norm_new2, c(1:length(sim_norm_new2)))
  # look at those high powered WMW sample sizes 
  as.data.frame(sim_norm_new2[[1]]$actual_power[which(sim_norm_new2[[1]]$method == "WMW", arr.ind = TRUE),]) # the power is high for all methods


# with small treatment effect
sim_norm_new3 <- vary_norm_2(n_pilot = 1000, rep=1, theta = log(1.1))
sim_plots(sim_norm_new3, c(1:length(sim_norm_new3)))

# with 3 categories
sim_norm_new4 <- vary_norm_2(n_pilot = 1000, rep=1, theta = log(1.8), vec_length = 3)
sim_plots(sim_norm_new4, c(1:length(sim_norm_new4)))

### different number of categories #####
vary_prob_length <- function(n_pilot = 1000, niter = 10000, r = 1, seed = 1){
  vec_length<- c(3:15)
  results <- vector(mode = "list", length = length(vec_length))
  for (i in seq_along(vec_length)) {
    set.seed(seed)
    p_C_u <- uniform_simplex(vec_length[i])[[1]]
    p_E_u <- uniform_simplex(vec_length[i])[[2]]
    sim_u <- simulation(p_C_u, p_E_u, n_pilot=n_pilot, niter=niter, r=r)
    results[[i]] <- sim_u
  }
  return(results)
}

sim_cat <- vary_prob_length(niter=10000, n_pilot=10000, r = 1)
mean_sim_cat <- c()
for (i in seq_along(1:13)) {
  mean_sim_cat[[i]] <- mean(sim_cat[[i]]$actual_power_nmin)
}
mean_sim_cat

sim_plots(sim_cat, c(3:15))


###
sim_cat2 <- vary_prob_length(niter=10000, n_pilot=10000, r = 1, seed=4)
mean_sim_cat2 <- c()
for (i in seq_along(1:13)) {
  mean_sim_cat2[[i]] <- mean(sim_cat2[[i]]$actual_power_nmin)
}
mean_sim_cat2
sim_plots(sim_cat2, c(3:15))

#### try it with fixed effect size for p_E
vary_prob_length_PO <- function(n_pilot = 1000, niter = 10000, r = 1, seed = 1, theta=log(1.8)){
  vec_length<- c(3:15)
  results <- vector(mode = "list", length = length(vec_length))
  for (i in seq_along(vec_length)) {
    set.seed(seed)
    p_C_u <- uniform_simplex(vec_length[i])[[1]]
    p_E_u <- calc_p_E(p_C_u, theta_A = theta)
    sim_u <- simulation(p_C_u, p_E_u, n_pilot=n_pilot, niter=niter, r=r)
    results[[i]] <- sim_u
  }
  return(results)
}

sim_cat_PO <- vary_prob_length_PO(n_pilot = 1000)
sim_cat_PO[[13]] <- NULL
sim_plots(sim_cat_PO, c(3:15))

sim_cat_PO2 <- vary_prob_length_PO(n_pilot = 1000, theta = log(1.3))
sim_plots(sim_cat_PO2, c(3:15))

sim_cat_PO3 <- vary_prob_length_PO(n_pilot = 100, theta = log(1.8))
sim_plots(sim_cat_PO3, c(3:15))

### make a scatter plot
# make a scatter plot with the mean actual Power
mean_min_power_cat <- c()
for(i in seq(1:length(sim_cat_PO))){
  mean_min_power_cat[i] <- mean(sim_cat_PO[[i]]$actual_power_nmin)
}

sim_cat_df <- data.frame("p"=c(3:15), "nmin_power"=mean_min_power_cat)
ggplot(sim_cat_df, aes(x=p, y=nmin_power)) +
  geom_point()+
  geom_hline(yintercept=0.8, color = "red", linetype = "dotted")+
  xlab("Number of Categories")+
  ylab("Mean Actual Power")+
  theme_bw()

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

sim_plots(sim_npilot, c(50, 100, 200, 500, 1000, 2500, 5000, 10000, 15000))
  
  
mean_sim_npilot <- c()
for (i in seq_along(1:length(sim_npilot))) {
  mean_sim_npilot[i] <- mean(sim_npilot[[i]]$actual_power_nmin)
}
mean_sim_npilot


sd_sim_npilot <- c()
for (i in seq_along(1:length(sim_npilot))) {
  sd_sim_npilot[i] <- sd(sim_npilot[[i]]$actual_power_nmin)
}
sd_sim_npilot

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

# mean and sd 
m <- c()
for(i in seq(1:length(sim_niter))){
  m[i] <- mean(sim_niter[[i]]$actual_power_nmin)
}
mean(m)

for(i in seq(1:length(sim_niter))){
  sd_i <- sd(sim_niter[[i]]$actual_power_nmin)
  print(sd_i)
}
mean(sd_i)

# proportions of the methods
prop <- data.frame("PO"=c(), "ttest"=c(), "WMW"=c())
for (i in seq(1:length(sim_niter))) {
  prop[i,1] <- length(which(sim_niter[[i]]$method == "PO"))/c(100, 500, 1000, 2500, 5000, 7500, 10000, 15000)[i]
  prop[i,2] <- length(which(sim_niter[[i]]$method == "ttest"))/c(100, 500, 1000, 2500, 5000, 7500, 10000, 15000)[i]
  prop[i,3] <- length(which(sim_niter[[i]]$method == "WMW"))/c(100, 500, 1000, 2500, 5000, 7500, 10000, 15000)[i]
}
prop


