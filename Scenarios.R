# Some real life scenarios are tested here
  # When one category wasn´t observed (meaning the prob. of this category is 0)
  # When the sample size is very big, because the effect is so small

### do some inspecting first
# mean sample size of sample sizes with power close to 1
mean(npilot_100$n_needed[which(npilot_100$actual_power_nmin > 0.99999),1])


#### Missing Category
# simulation that handles when n_pilot is very small
# either discard the one category or combine it with another neighbor category
simulation_missing_cat<-function(p_C, p_E, r=1, n_pilot, niter=1000, alpha=0.05,beta=0.2){
  n_needed<-matrix(NA,niter,3)
  colnames(n_needed)<-c("AfS","ttestord", "po")
  actual_power2<-matrix(NA,niter,3)
  colnames(actual_power2)<-c("AfS","ttestord", "po")
  for (i in 1:niter){
    set.seed(i)
    print(i)
    phat<-generate.pEpC.estimate(p_C=p_C, p_E=p_E, r=r, n_pilot=n_pilot)
    phat_df <- data.frame(phat$C, phat$E)
    if(any(c(phat$C, phat$E)==0)){
      null_cat <- which(phat_df==0, arr.ind = TRUE)[,1]
      if(any(null_cat>1)){
        combine_cat <- null_cat-1
        phat_df[combine_cat,] <- phat_df[null_cat,] + phat_df[combine_cat,]
        phat_df<-phat_df[-null_cat,]
        phat <- list(C=phat_df$phat.C, E=phat_df$phat.E)
        n_pilot <- n_pilot-1
      }else{
        combine_cat <- null_cat+1
        phat_df[combine_cat,] <- phat_df[null_cat,] + phat_df[combine_cat,]
        phat_df<-phat_df[-null_cat,]
        phat <- list(C=phat_df$phat.C, E=phat_df$phat.E)
        n_pilot <- n_pilot-1
      }
    } else {}
    result_AfS<-samplesize_AfS(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
    result_ttestord<-samplesize_ttestord(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
    result_po <- samplesize_po_NN(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E) # newly added
    n_needed[i,1]<- result_AfS$n_total
    n_needed[i,2]<- result_ttestord$n_total
    n_needed[i,3]<- result_po$n_total
    actual_power2[i,1]<-result_AfS$actual_power2
    actual_power2[i,2]<-result_ttestord$actual_power2
    actual_power2[i,3]<-result_po$actual_power2
  }
  nmin<-apply(n_needed,FUN=min,MARGIN=1)
  whichnmin<-apply(n_needed,FUN=which.min,MARGIN=1)
  actual_power_nmin<-actual_power2[cbind(1:niter,whichnmin)]
  list(n_needed=n_needed,actual_power = actual_power2, 
       nmin=nmin,actual_power_nmin=actual_power_nmin,
       method = c("AfS", "ttest", "PO")[whichnmin])
}
test_miss <- simulation_missing_cat(p_C, p_E, r=1, niter = 10000, n_pilot = 30)


#### Enormous sample size
# simulation that handles when the effect is very small, and the sample size thus very big
# This should also get rid of the many Power=1

### leave out study when effect is very small
simulation_small_effect<-function(p_C, p_E, r=1, n_pilot, niter=1000, alpha=0.05,beta=0.2){
  n_needed<-matrix(NA,niter,4)
  colnames(n_needed)<-c("AfS","ttestord", "po", "index")
  actual_power2<-matrix(NA,niter,4)
  colnames(actual_power2)<-c("AfS","ttestord", "po", "index")
  est_eff <- matrix(NA,niter,2)
  colnames(est_eff)<-c("log-odds-ratio", "index")
  for (i in 1:niter){
    set.seed(i)
    print(i)
    phat<-generate.pEpC.estimate(p_C=p_C, p_E=p_E, r=r, n_pilot=n_pilot)
    est_effect <- calculate_theta_A_old(p_C=phat$C, p_E=phat$E, r=r)[,1]
    est_eff[i,1] <- est_effect
    est_eff[i,2] <- i
    #    if(abs(est_effect)<0.1){    print ("too small effect") ; next }
    if(mean(phat$C-phat$E)==0){    print ("diff is 0") ; next }
    result_AfS<-samplesize_AfS(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
    result_ttestord<-samplesize_ttestord(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
    result_po <- samplesize_po_NN(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
    n_needed[i,1]<- result_AfS$n_total
    n_needed[i,2]<- result_ttestord$n_total
    n_needed[i,3]<- result_po$n_total
    n_needed[i,4]<- i
    actual_power2[i,1]<-result_AfS$actual_power2
    actual_power2[i,2]<-result_ttestord$actual_power2
    actual_power2[i,3]<-result_po$actual_power2
    actual_power2[i,4]<-i
  }
  n_needed <- as.data.frame(na.omit(n_needed))
  est_eff <- as.data.frame(na.omit(est_eff))
  actual_power2 <- as.data.frame(na.omit(actual_power2))
  actual_power2 <- actual_power2[n_needed$index,]
  nmin<-apply(n_needed,FUN=min,MARGIN=1)
  whichnmin<-apply(n_needed[,c(1,2,3)],FUN=which.min,MARGIN=1)
  actual_power_nmin<-actual_power2[cbind(1:length(whichnmin),whichnmin)]
  list(n_needed=n_needed,actual_power = actual_power2, 
       nmin=nmin,actual_power_nmin=actual_power_nmin,
       method = c("AfS", "ttest", "PO")[whichnmin],
       effect = est_eff)
}

test_small_effect <- simulation_small_effect(p_C, p_E, r=1, niter = 10000, n_pilot = 50)
hist(test_small_effect$actual_power_nmin)
hist(npilot_50$actual_power_nmin)
length(which(test_small_effect$actual_power_nmin > 0.96))
length(which(npilot_50$actual_power_nmin > 0.96))
# how many results have a power of 1 
nrow(test_small_effect$effect[which(test_small_effect$actual_power_nmin == 1),])




### leave out study when calculated sample size is very large
simulation_big_samp<-function(p_C, p_E, r=1, n_pilot, niter=1000, alpha=0.05,beta=0.2, max_samp = 1000){
  n_needed<-matrix(NA,niter,4)
  colnames(n_needed)<-c("AfS","ttestord", "po", "index")
  actual_power2<-matrix(NA,niter,4)
  colnames(actual_power2)<-c("AfS","ttestord", "po", "index")
  for (i in 1:niter){
    set.seed(i)
    print(i)
    phat<-generate.pEpC.estimate(p_C=p_C, p_E=p_E, r=r, n_pilot=n_pilot)
  # calc. sample sizes
    result_AfS<-samplesize_AfS(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
    result_ttestord<-samplesize_ttestord(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
    result_po <- samplesize_po_NN(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
    if(any(c(result_AfS$n_total, result_ttestord$n_total, result_po$n_total)>max_samp)){    print ("sample size too big") ; next }
    n_needed[i,1]<- result_AfS$n_total
    n_needed[i,2]<- result_ttestord$n_total
    n_needed[i,3]<- result_po$n_total
    n_needed[i,4]<- i
    actual_power2[i,1]<-result_AfS$actual_power2
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
       method = c("AfS", "ttest", "PO")[whichnmin])
}

test_big_samp <- simulation_big_samp(p_C, p_E, r=1, niter = 10000, n_pilot = 100, max_samp = 1500)
hist(test_big_samp$actual_power_nmin)
hist(npilot_50$actual_power_nmin)
hist(npilot_100$actual_power_nmin)
length(which(test_big_samp$actual_power_nmin > 0.9))
length(which(npilot_100$actual_power_nmin > 0.9))
# mean sample size of sample sizes with power close to 1
mean(test_big_samp$n_needed[which(test_big_samp$actual_power_nmin > 0.99999),1])

# only leaving out the incredibly big ones
test_very_big_samp <- simulation_big_samp(p_C, p_E, r=1, niter = 10000, n_pilot = 100, max_samp = 100000)
mean(test_very_big_samp$n_needed[which(test_very_big_samp$actual_power_nmin > 0.99999),1])
  ### so not only the final sample size is an indicator of power = 1

# do this for different n_pilot
vary_npilot_exclusion <- function(p_C, p_E, niter = 10000, r = 1, max_samp = 1000){
  npilot_vec <- c(50, 100, 200, 500, 1000, 2500, 5000, 10000, 15000, 20000)
  results <- vector(mode = "list", length = length(npilot_vec))
  for (i in seq_along(npilot_vec)) {
    sim <- simulation_big_samp(p_C, p_E, r = r, niter = niter, n_pilot = npilot_vec[i], max_samp = max_samp )
    results[[i]] <- sim
  }
  return(results)
}
sim_ss_exclusion <- vary_npilot_exclusion(p_C, p_E, max_samp=800)
sim_plots(sim_ss_exclusion, c(50, 100, 200, 500, 1000, 2500, 5000, 10000, 15000, 20000))
sim_plots(sim_npilot, c(50, 100, 200, 500, 1000, 2500, 5000, 10000, 15000, 20000))
# way less power = 1 if sample sizes bigger than 1500 are excluded


### try to exclude the iterations, in which phat is very far off from p_C

simulation_off_phat<-function(p_C, p_E, r=1, n_pilot=1000, niter=1000, alpha=0.05,beta=0.2){
  n_needed<-matrix(NA,niter,4)
  colnames(n_needed)<-c("AfS","ttestord", "po", "index")
  actual_power2<-matrix(NA,niter,4)
  colnames(actual_power2)<-c("AfS","ttestord", "po", "index")
  for (i in 1:niter){
    set.seed(i)
    print(i)
    phat<-generate.pEpC.estimate(p_C=p_C, p_E=p_E, r=r, n_pilot=n_pilot)
    if(sd(phat$E-p_E)&sd(phat$C-p_C)>0.1){    print ("phat off") ; next }
    # calc. sample sizes
    result_AfS<-samplesize_AfS(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
    result_ttestord<-samplesize_ttestord(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
    result_po <- samplesize_po_NN(p_C=phat$C, p_E=phat$E, alpha=alpha, beta=beta, r=r, p_C2=p_C, p_E2=p_E)
    n_needed[i,1]<- result_AfS$n_total
    n_needed[i,2]<- result_ttestord$n_total
    n_needed[i,3]<- result_po$n_total
    n_needed[i,4]<- i
    actual_power2[i,1]<-result_AfS$actual_power2
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
       method = c("AfS", "ttest", "PO")[whichnmin])
}


test_big_samp <- simulation_off_phat(p_C, p_E, r=1, niter = 10000, n_pilot = 100)
hist(test_big_samp$actual_power_nmin)
hist(npilot_100$actual_power_nmin)



