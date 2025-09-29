##### Sensitivitätsanalyse, schließe pilotstudien mit großen und kleinen effekten aus

simulation_sens<-function(p_C, p_E, r=1, n_pilot, niter=1000, alpha=0.05,beta=0.2, max_samp = 1000){
  n_needed<-matrix(NA,niter,4)
  colnames(n_needed)<-c("AfS","ttestord", "po", "index")
  actual_power2<-matrix(NA,niter,4)
  colnames(actual_power2)<-c("AfS","ttestord", "po", "index")
  for (i in 1:niter){
    set.seed(i)
    print(i)
    phat<-generate.pEpC.estimate(p_C=p_C, p_E=p_E, r=r, n_pilot=n_pilot)
    theta_A <- calculate_theta_A(p_C=phat$C, p_E=phat$E, r)[1,1]
    if(theta_A<0.2 && theta_A>-0.2){    print ("effect too small") ; next }

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


test_sens <- simulation_sens(p_C, p_E, n_pilot=100, niter=10000)
mean(test_sens$actual_power_nmin)
hist(test_sens$actual_power_nmin)

# without exclusion
hist(test100$actual_power_nmin)
mean(test100$actual_power_nmin)

## exclude pilot studies that are too far off from true probabilities
simulation_diff_effect<-function(p_C, p_E, r=1, n_pilot, niter=1000, alpha=0.05,beta=0.2, effect_cut){
  n_needed<-matrix(NA,niter,4)
  colnames(n_needed)<-c("AfS","ttestord", "po", "index")
  actual_power2<-matrix(NA,niter,4)
  colnames(actual_power2)<-c("AfS","ttestord", "po", "index")
  theta_p <- calculate_theta_A(p_C, p_E, r)[1,1] 
  for (i in 1:niter){
    set.seed(i)
    print(i)
    phat<-generate.pEpC.estimate(p_C=p_C, p_E=p_E, r=r, n_pilot=n_pilot)
    theta_phat <- calculate_theta_A(p_C=phat$C, p_E=phat$E, r)[1,1]
    theta_diff <- abs(theta_phat-theta_p)
    if(theta_diff>effect_cut){    print ("effect too small") ; next }
    
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
  nmin<-apply(n_needed,FUN=min,MARGIN=1)
  whichnmin<-apply(n_needed[,c(1,2,3)],FUN=which.min,MARGIN=1)
  actual_power_nmin<-actual_power2[cbind(1:length(whichnmin),whichnmin)]
  list(n_needed=n_needed,actual_power = actual_power2, 
       nmin=nmin,actual_power_nmin=actual_power_nmin,
       method = c("AfS", "ttest", "PO")[whichnmin])
}
test_eff_cut <- simulation_diff_effect(p_C, p_E, n_pilot=100, niter=10000, effect_cut = 0.15)
length(test_eff_cut$actual_power_nmin)
mean(test_eff_cut$actual_power_nmin)
hist(test_eff_cut$actual_power_nmin)

df_eff_cut <-  data.frame("power"=test_eff_cut$actual_power_nmin)%>% 
  mutate(mean = mean(power))
df_eff_no_cut <-  data.frame("power_no_cut"=test100$actual_power_nmin)%>% 
  mutate(mean = mean(power_no_cut))
ggplot(df_eff_cut, aes(x = power)) + 
  geom_histogram(data=df_eff_no_cut, mapping=aes(x=power_no_cut),alpha=0.15, fill = "red", color = "darkgrey",position = "identity") +
  geom_histogram(fill = "grey", color = "black",position = "identity") +
  geom_vline(aes(xintercept = mean), color = "red")+
  geom_vline(data=df_eff_no_cut, aes(xintercept = mean), color = "red", alpha=0.5)+
  geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
  xlab("Actual Power")+
  ylab("Count")+
  theme_bw()

# without exclusion
hist(test100$actual_power_nmin)
mean(test100$actual_power_nmin)



test_eff_cut2 <- simulation_diff_effect(p_C, p_E, n_pilot=50, niter=10000, effect_cut = 0.15)
length(test_eff_cut2$actual_power_nmin)
mean(test_eff_cut2$actual_power_nmin)
hist(test_eff_cut2$actual_power_nmin)
