
### make simulation that chooses biggest sample size
  # look if the power values are better, or in some cases, the same

simulation_biggest<-function(p_C, p_E, r=1, n_pilot, niter=1000, alpha=0.05,beta=0.2){
  n_needed<-matrix(NA,niter,3)
  colnames(n_needed)<-c("WMW","ttestord", "po")
  actual_power2<-matrix(NA,niter,3)
  colnames(actual_power2)<-c("WMW","ttestord", "po")
  for (i in 1:niter){
    set.seed(i)
    print(i)
    phat<-generate.pEpC.estimate(p_C=p_C, p_E=p_E, r=r, n_pilot=n_pilot)
    result_WMW<-samplesize_AfS(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
    result_ttestord<-samplesize_ttestord(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
    result_po <- samplesize_po_NN(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E) # newly added
    n_needed[i,1]<- result_WMW$n_total
    n_needed[i,2]<- result_ttestord$n_total
    n_needed[i,3]<- result_po$n_total
    actual_power2[i,1]<-result_WMW$actual_power2
    actual_power2[i,2]<-result_ttestord$actual_power2
    actual_power2[i,3]<-result_po$actual_power2
  }
  nmax<-apply(n_needed,FUN=max,MARGIN=1)
  whichnmax<-apply(n_needed,FUN=which.max,MARGIN=1)
  actual_power_nmax<-actual_power2[cbind(1:niter,whichnmax)]
  list(n_needed=n_needed,actual_power = actual_power2, 
       nmax=nmax,actual_power_nmax=actual_power_nmax,
       method = c("WMW", "ttest", "PO")[whichnmax])
}



### 
max1000<-simulation_biggest(p_C, p_E, n_pilot=1000, niter=10000)
mean(max1000$actual_power_nmax)
sd(max1000$actual_power_nmax)
hist(max1000$actual_power_nmax)

df_max1000 <- data.frame()
df_max1000 <- data.frame("power_nmax" = max1000$actual_power_nmax, 
                          "method" = max1000$method) %>%
  mutate(mean = mean(power_nmax))
df_max1000$method <- as.factor(df_max1000$method)
df_max1000$method <- relevel(df_max1000$method, "ttest")


grid.arrange(ggplot(df_max1000, aes(x = power_nmax))+ 
               geom_histogram(fill = "grey", color = "black", position = "identity")+
               geom_vline(aes(xintercept = mean), color = "red")+
               geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
               theme_bw(),
             ggplot(df_max1000, aes(x = power_nmax, fill = method, colour = method))+ 
               geom_histogram(alpha = 0.3, position = "identity")+
               geom_vline(aes(xintercept = mean), color = "red")+
               scale_colour_manual(values = c("PO" = "#7CAE00", "ttest" = "#00BFC4", "WMW" = "#C77CFF")) +
               scale_fill_manual(values = c("PO" = "#7CAE00", "ttest" = "#00BFC4", "WMW" = "#C77CFF")) +
               geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
               theme_bw())

sim_plots_biggest <- function(data, group, plots=0){
  df_data <- data.frame()
  for (i in seq_along(1:length(data))) {
    df <- data.frame("power_nmax" = data[[i]]$actual_power_nmax, 
                     "method" = data[[i]]$method,
                     "group" = group[i])
    df_data <- rbind(df_data,df)
  }
  
  df_data <- df_data %>% group_by(group) %>%  mutate(mean = mean(power_nmax))
  df_data$method <- as.factor(df_data$method)
  df_data$method <- relevel(df_data$method, "ttest")
  
  df_means <- df_data %>%
    group_by(group) %>%
    summarise(mean_val = mean(power_nmax, na.rm = TRUE))
  
  plot1 <- ggplot(df_data, aes(x = power_nmax))+ 
    geom_histogram(fill = "grey", color = "black", position = "identity")+
    geom_vline(data = df_means, aes(xintercept = mean_val), color = "red")+
    geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
    facet_wrap(group ~.) +
    theme_bw()
  plot2 <- ggplot(df_data, aes(x = power_nmax, fill = method, colour = method))+ 
    facet_wrap(group ~.) +
    geom_histogram(alpha = 0.2, position = "identity")+
    scale_colour_manual(values = c("PO" = "#7CAE00", "ttest" = "#00BFC4", "WMW" = "#C77CFF")) +
    scale_fill_manual(values = c("PO" = "#7CAE00", "ttest" = "#00BFC4", "WMW" = "#C77CFF")) +
    geom_vline(data = df_means, aes(xintercept = mean_val), color = "red")+
    geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
    theme_bw()
  if(plots==0){
    grid.arrange(plot1, plot2, ncol=2)
  } else if(plots==1){
    return(plot1)
  } else if(plots == 2){
    return(plot2)
  }
}

### vary n_pilot for biggest sample size
vary_npilot_biggest <- function(p_C, p_E, niter = 10000, r = 1, npilot_vec=c(50, 100, 200, 500, 1000, 2500, 5000, 10000, 15000)){
  results <- vector(mode = "list", length = length(npilot_vec))
  for (i in seq_along(npilot_vec)) {
    sim <- simulation_biggest(p_C, p_E, r = r, niter = niter, n_pilot = npilot_vec[i])
    results[[i]] <- sim
  }
  return(results)
}

sim_npilot_biggest_normal <- vary_npilot_biggest(p_C, p_E)
sim_npilot_biggest <- vary_npilot_biggest(p_C, p_E, npilot_vec = c(seq(150,1000,50), seq(1500, 12000, 500)))

sim_plots_biggest(sim_npilot_biggest_normal, c(50, 100, 200, 500, 1000, 2500, 5000, 10000, 15000, 20000))



### scatter plot for min and max sample sizes
mean_max_power_npilot <- c()
for(i in seq(1:length(sim_npilot_biggest))){
  mean_max_power_npilot[i] <- mean(sim_npilot_biggest[[i]]$actual_power_nmax, na.rm = TRUE)
}


sim_npilot_min_max <- data.frame("npilot"=c(seq(150,1000,50), seq(1500, 12000, 500)), 
                                 "min"=sim_small_test[[1]][,2], 
                                 "max"=mean_max_power_npilot) %>%
  pivot_longer(!npilot, names_to = "min_or_max", values_to = "power")
ggplot(sim_npilot_min_max) +
  geom_point(mapping = aes(x=npilot, y=power, colour = min_or_max))+
  geom_line(mapping = aes(x=npilot, y=power, colour = min_or_max))+
  geom_hline(yintercept=0.8, color = "red", linetype = "dotted")+
  scale_colour_manual(values = c("min" = "#CD661D", "max" = "#191970")) +
  xlab("Pilot Study Sample Size")+
  ylab("Mean Actual Power")+
  theme_bw()
