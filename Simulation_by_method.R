#### AfS #####
simulation_AfS<-function(p_C, p_E, r=1, n_pilot, niter=1000, alpha=0.05,beta=0.2)
{
  n_needed<-c()
  actual_power2<- c()
  for (i in 1:niter){
    set.seed(i)
    print(i)
    phat<-generate.pEpC.estimate(p_C=p_C, p_E=p_E, r=r, n_pilot=n_pilot)
    result_AfS<-samplesize_AfS(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
    n_needed[i]<- result_AfS$n_total
    actual_power2[i]<-result_AfS$actual_power2
  }
  nmin<-min(n_needed)
  whichnmin<-which.min(n_needed)
  actual_power_nmin<-actual_power2[whichnmin]
  list(n_needed=n_needed,actual_power = actual_power2, nmin=nmin,actual_power_nmin=actual_power_nmin)
}
test200_AfS<-simulation_AfS(p_C,p_E,n_pilot=200,niter=1000)
test200_AfS$nmin
hist(test200_AfS$actual_power)
hist(test200_AfS$n_needed)

sim_mean_power <- data.frame("n_pilot"=c(NA), "nmin_power"=c(NA), "nmin"=c(NA), "nmean" =c(NA))
settings_AfS <- function(n_vector, p_C, p_E, r, niter){
  for (i in seq_along(n_vector)) {
    sim <- simulation_AfS(p_C=p_C, p_E=p_E, r=r, n_pilot = n_vector[i], niter=niter)
    sim_mean_power[i,2] <- sim$actual_power_nmin
    sim_mean_power[i,1] <- n_vector[i]
    sim_mean_power[i,3] <- sim$nmin
    sim_mean_power[i,4] <- mean(sim$n_needed)}
  return(sim_mean_power)
}
sim_test_AfS <- settings_AfS(n_vector = c(seq(50,1000,50), seq(1500, 12000, 500)), p_C, p_E, r, niter = 1000)
ggplot(sim_test_AfS, aes(x=n_pilot, y=nmin)) +
  geom_point()
ggplot(sim_test_AfS, aes(x=n_pilot, y=nmin_power)) +
  geom_point()



#### ttest #####

simulation_ttest<-function(p_C, p_E, r=1, n_pilot, niter=1000, alpha=0.05,beta=0.2)
{
  n_needed<-c()
  actual_power2<- c()
  for (i in 1:niter){
    set.seed(i)
    print(i)
    phat<-generate.pEpC.estimate(p_C=p_C, p_E=p_E, r=r, n_pilot=n_pilot)
    result_ttest<-samplesize_ttestord(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
    n_needed[i]<- result_ttest$n_total
    actual_power2[i]<-result_ttest$actual_power2
  }
  nmin<-min(n_needed)
  whichnmin<-which.min(n_needed)
  actual_power_nmin<-actual_power2[whichnmin]
  list(n_needed=n_needed,actual_power = actual_power2, nmin=nmin,actual_power_nmin=actual_power_nmin)
}
test200_ttest<-simulation_ttest(p_C,p_E,n_pilot=200,niter=1000)

settings_ttest <- function(n_vector, p_C, p_E, r, niter){
  for (i in seq_along(n_vector)) {
    sim <- simulation_ttest(p_C=p_C, p_E=p_E, r=r, n_pilot = n_vector[i], niter=niter)
    sim_mean_power[i,2] <- sim$actual_power_nmin
    sim_mean_power[i,1] <- n_vector[i]    
    sim_mean_power[i,3] <- sim$nmin
    sim_mean_power[i,4] <- mean(sim$n_needed)}
  return(sim_mean_power)
}
sim_test_ttest <- settings_ttest(n_vector = c(seq(50,1000,50), seq(1500, 12000, 500)), p_C, p_E, r, niter = 1000)
ggplot(sim_test_ttest, aes(x=n_pilot, y=nmin_power)) +
  geom_point()
ggplot(sim_test_ttest, aes(x=n_pilot, y=nmin)) +
  geom_point()


#### PO NN ####
simulation_po_NN<-function(p_C, p_E, r=1, n_pilot, niter=1000, alpha=0.05,beta=0.2)
{
  n_needed<-c()
  actual_power2<- c()
  for (i in 1:niter){
    set.seed(i)
    print(i)
    phat<-generate.pEpC.estimate(p_C=p_C, p_E=p_E, r=r, n_pilot=n_pilot)
    result_po<-samplesize_po_NN(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
    n_needed[i]<- result_po$n_total
    actual_power2[i]<-result_po$actual_power2
  }
  nmin<-min(n_needed)
  whichnmin<-which.min(n_needed)
  actual_power_nmin<-actual_power2[whichnmin]
  list(n_needed=n_needed,actual_power = actual_power2, nmin=nmin,actual_power_nmin=actual_power_nmin)
}
test200_po<-simulation_po_NN(p_C,p_E,n_pilot=200,niter=1000)

settings_po_NN <- function(n_vector, p_C, p_E, r, niter){
  for (i in seq_along(n_vector)) {
    sim <- simulation_po_NN(p_C=p_C, p_E=p_E, r=r, n_pilot = n_vector[i], niter=niter)
    sim_mean_power[i,2] <- sim$actual_power_nmin
    sim_mean_power[i,1] <- n_vector[i]    
    sim_mean_power[i,3] <- sim$nmin
    sim_mean_power[i,4] <- mean(sim$n_needed)}
  return(sim_mean_power)
}
sim_test_po_NN <- settings_po_NN(c(seq(50,1000,50), seq(1500, 12000, 500)), p_C, p_E, r, niter = 1000)
ggplot(sim_test_po_NN, aes(x=n_pilot, y=nmin_power)) +
  geom_point()
ggplot(sim_test_po_NN, aes(x=n_pilot, y=nmin)) +
  geom_point()


#### PO AA ####
simulation_po_AA<-function(p_C, p_E, r=1, n_pilot, niter=1000, alpha=0.05,beta=0.2)
{
  n_needed<-c()
  actual_power2<- c()
  for (i in 1:niter){
    set.seed(i)
    print(i)
    phat<-generate.pEpC.estimate(p_C=p_C, p_E=p_E, r=r, n_pilot=n_pilot)
    result_po<-samplesize_po_AA(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
    n_needed[i]<- result_po$n_total
    actual_power2[i]<-result_po$actual_power2
  }
  nmin<-min(n_needed)
  whichnmin<-which.min(n_needed)
  actual_power_nmin<-actual_power2[whichnmin]
  list(n_needed=n_needed,actual_power = actual_power2, nmin=nmin,actual_power_nmin=actual_power_nmin)
}
test200_po_AA<-simulation_po_AA(p_C,p_E,n_pilot=200,niter=1000)

settings_po_AA <- function(n_vector, p_C, p_E, r, niter){
  for (i in seq_along(n_vector)) {
    sim <- simulation_po_AA(p_C=p_C, p_E=p_E, r=r, n_pilot = n_vector[i], niter=niter)
    sim_mean_power[i,2] <- sim$actual_power_nmin
    sim_mean_power[i,1] <- n_vector[i]    
    sim_mean_power[i,3] <- sim$nmin
    sim_mean_power[i,4] <- mean(sim$n_needed)}
  return(sim_mean_power)
}
sim_test_po_AA <- settings_po_AA(c(seq(50,1000,50), seq(1500, 12000, 500)), p_C, p_E, r, niter = 1000)
ggplot(sim_test_po_AA, aes(x=n_pilot, y=nmin_power)) +
  geom_point()
ggplot(sim_test_po_AA, aes(x=n_pilot, y=nmin)) +
  geom_point()


#### PO Kieser ####

simulation_po_kieser<-function(p_C, p_E, r=1, n_pilot, niter=1000, alpha=0.05,beta=0.2){
  n_needed<-c()
  actual_power2<- c()
  for (i in 1:niter){
    set.seed(i)
    print(i)
    phat<-generate.pEpC.estimate(p_C=p_C, p_E=p_E, r=r, n_pilot=n_pilot)
    result_po<-samplesize_po_kieser(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
    n_needed[i]<- result_po$n_total
    actual_power2[i]<-result_po$actual_power2
  }
  nmin<-min(n_needed)
  whichnmin<-which.min(n_needed)
  actual_power_nmin<-actual_power2[whichnmin]
  list(n_needed=n_needed,actual_power = actual_power2, nmin=nmin,actual_power_nmin=actual_power_nmin)
}
test200_po_kieser<-simulation_po_kieser(p_C,p_E,n_pilot=200,niter=1000)

settings_po_kieser <- function(n_vector, p_C, p_E, r, niter){
  for (i in seq_along(n_vector)) {
    sim <- simulation_po_kieser(p_C=p_C, p_E=p_E, r=r, n_pilot = n_vector[i], niter=niter)
    sim_mean_power[i,2] <- sim$actual_power_nmin
    sim_mean_power[i,1] <- n_vector[i]    
    sim_mean_power[i,3] <- sim$nmin
    sim_mean_power[i,4] <- mean(sim$n_needed)}
  return(sim_mean_power)
}
sim_test_po_kieser <- settings_po_kieser(c(seq(50,1000,50), seq(1500, 12000, 500)), p_C, p_E, r, niter = 1000)

ggplot(sim_test_po_kieser, aes(x=n_pilot, y=nmin_power)) +
  geom_point()
ggplot(sim_test_po_kieser, aes(x=n_pilot, y=nmin)) +
  geom_point()


#### PO NA ####

simulation_po_NA<-function(p_C, p_E, r=1, n_pilot, niter=1000, alpha=0.05,beta=0.2)
{
  n_needed<-c()
  actual_power2<- c()
  for (i in 1:niter){
    set.seed(i)
    print(i)
    phat<-generate.pEpC.estimate(p_C=p_C, p_E=p_E, r=r, n_pilot=n_pilot)
    result_po<-samplesize_po_NA(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
    n_needed[i]<- result_po$n_total
    actual_power2[i]<-result_po$actual_power2
  }
  nmin<-min(n_needed)
  whichnmin<-which.min(n_needed)
  actual_power_nmin<-actual_power2[whichnmin]
  list(n_needed=n_needed,actual_power = actual_power2, nmin=nmin,actual_power_nmin=actual_power_nmin)
}
test200_po_NA<-simulation_po_NA(p_C,p_E,n_pilot=200,niter=1000)

settings_po_NA <- function(n_vector, p_C, p_E, r, niter){
  for (i in seq_along(n_vector)) {
    sim <- simulation_po_NA(p_C=p_C, p_E=p_E, r=r, n_pilot = n_vector[i], niter=niter)
    sim_mean_power[i,2] <- sim$actual_power_nmin
    sim_mean_power[i,1] <- n_vector[i]    
    sim_mean_power[i,3] <- sim$nmin
    sim_mean_power[i,4] <- mean(sim$n_needed)}
  return(sim_mean_power)
}
sim_test_po_NA <- settings_po_NA(c(seq(50,1000,50), seq(1500, 12000, 500)), p_C, p_E, r, niter = 1000)
ggplot(sim_test_po_NA, aes(x=n_pilot, y=nmin_power)) +
  geom_point()
ggplot(sim_test_po_NA, aes(x=n_pilot, y=nmin)) +
  geom_point()


# Combine All ####


sim_combine <- rbind(cbind(sim_test_AfS, "method" = rep("AfS", nrow(sim_test_AfS))),
                     cbind(sim_test_ttest, "method" = rep("ttest", nrow(sim_test_ttest))),
                     cbind(sim_test_po_NN, "method" = rep("po_NN", nrow(sim_test_po_NN))),
                     cbind(sim_test_po_kieser, "method" = rep("po_kieser", nrow(sim_test_po_kieser))),
                     cbind(sim_test_po_NA, "method" = rep("po_NA", nrow(sim_test_po_NA))),
                     cbind(sim_test_po_AA, "method" = rep("po_AA", nrow(sim_test_po_AA))))
ggplot(sim_combine, aes(x=n_pilot, y=nmin)) +
  geom_point(aes(color = method), alpha = 0.7)

sim_combine_po <- rbind(cbind(sim_test_po_NN, "method" = rep("po_NN", nrow(sim_test_po_NN))),
                     cbind(sim_test_po_kieser, "method" = rep("po_kieser", nrow(sim_test_po_kieser))),
                     cbind(sim_test_po_NA, "method" = rep("po_NA", nrow(sim_test_po_NA))),
                     cbind(sim_test_po_AA, "method" = rep("po_AA", nrow(sim_test_po_AA))))
ggplot(sim_combine_po, aes(x=n_pilot, y=nmin)) +
  geom_point(aes(color = method), alpha = 0.7) # kieser and NN method make smaller sample sizes than NA
ggplot(sim_combine_po, aes(x=n_pilot, y=nmin_power)) +
  geom_point(aes(color = method), alpha = 0.7) 

# How different are the sample sizes from White compared to Kieser?
# NN
mean(sim_test_po_kieser$nmin - sim_test_po_NN$nmin) # no difference in minimum sample size
mean(sim_test_po_kieser$nmean - sim_test_po_NN$nmean) # difference in sample size means
mean(sim_test_po_kieser$nmin_power - sim_test_po_NN$nmin_power) #very small difference in power calculation

# NA
mean(sim_test_po_kieser$nmin - sim_test_po_NA$nmin) # difference : -2.3
mean(sim_test_po_kieser$nmean - sim_test_po_NA$nmean) # big difference in sample size means (-27.36)
mean(sim_test_po_kieser$nmin_power - sim_test_po_NA$nmin_power) # difference: -0.0098

# AA
mean(sim_test_po_kieser$nmin - sim_test_po_AA$nmin) # difference : -7.7
mean(sim_test_po_kieser$nmean - sim_test_po_AA$nmean) # huge difference in sample size means (-95.76)
mean(sim_test_po_kieser$nmin_power - sim_test_po_AA$nmin_power) # difference: -0.5

##### Function that does all in one #####
simulation_method <- function(p_C, p_E, r=1, n_pilot, niter=1000, alpha=0.05,beta=0.2,n_min, n_max, steps) {
  set_AfS <- settings_AfS(n_min, n_max, steps, p_C, p_E, r, niter)
  set_ttest <- settings_ttest(n_min, n_max, steps, p_C, p_E, r, niter)
  set_po <- settings_po(n_min, n_max, steps, p_C, p_E, r, niter)
  sim_combine <- rbind(cbind(set_AfS, "method" = rep("AfS", 20)),
                       cbind(set_ttest, "method" = rep("ttest", 20)),
                       cbind(set_po, "method" = rep("po", 20)))
  plot <- ggplot(sim_combine, aes(x=n_pilot, y=nmin_power)) +
    geom_point(aes(color = method))
  return(plot)
}
sim_all_1 <- simulation_method(n_min = 100, n_max = 10000, steps = 500, p_C, p_E, r, niter = 10000)
sim_all_2 <- simulation_method(n_min = 100, n_max = 10000, steps = 500, p_C_2, p_E_2, r, niter = 10000)
sim_all_3 <- simulation_method(n_min = 100, n_max = 10000, steps = 500, p_C_3, p_E_3, r, niter = 10000)
sim_all_4 <- simulation_method(n_min = 100, n_max = 10000, steps = 500, p_C_4, p_E_4, r, niter = 10000)
sim_all_5 <- simulation_method(n_min = 100, n_max = 10000, steps = 500, p_C_5, p_E_5, r, niter = 10000)


