### Compare ttest function from R with my implementation
comp_ttest <- function(alpha = 0.05, beta = 0.2, r = 1, iter = 1000){
  results_R <-  c()
  results_me <-  c()
  for (i in seq_along(1:iter)) {
    set.seed(i)
    prob_length <- sample(c(3:14), 1)
    theta <- runif(1, min = 0.1, max = 5)
    print(i)
    p <- generate_two_simplex_vectors(prob_length, log(theta))
    p_C <- p[[1]]
    p_E <- p[[2]]
    ttest <- samplesize_ttestord(p_C=p_C, p_E=p_E, alpha, beta, r)
    results_me <- rbind(results_me, ttest$actual_power)
    
    del_sig <- calculate_delta_A_sigma(p_C, p_E)
    results_R <- c(results_R, power.t.test(n=ttest$n_E, delta = del_sig$delta_A, power = NULL,
                 sd = del_sig$sigma, alternative = "two.sided", type = "two.sample")$power)
        }
 return(data.frame("ttestord"=results_me, "R"=results_R)) 
}

# power.t.test(delta = del_sig$delta_A, power = 0.8, sig.level = 0.05, 
#             sd = del_sig$sigma, alternative = "two.sided", type = "two.sample")


comp_ttest <- comp_ttest(iter = 10000)
mean(abs(comp_ttest$ttestord-comp_ttest2$R)) # 0.004933035
ggplot(data=comp_ttest, aes(x=ttestord, y=R))+
  geom_point(alpha=0.6)+
  geom_abline(slope=1, color="red")+
  xlab("t-test")+
  ylab("base R")+
  theme_bw()
length(which(comp_ttest$R < comp_ttest$ttestord))/10000


# with (1+r) in the Denominator
comp_ttest_kieser <- function(alpha = 0.05, beta = 0.2, r = 1, iter = 1000, bias = 1.2){
  results_R <-  c()
  results_me <-  c()
  for (i in seq_along(1:iter)) {
    set.seed(i)
    prob_length <- sample(c(3,4,5,6,7,8), 1)
    print(i)
    p <- generate_two_simplex_vectors(prob_length, bias)
    p_C <- p[[1]]
    p_E <- p[[2]]
    ttest <- samplesize_ttestord_kieser(p_C=p_C, p_E=p_E, alpha, beta, r)
    results_me <- rbind(results_me, ttest$actual_power)
    
    del_sig <- calculate_delta_A_sigma(p_C, p_E)
    results_R <- c(results_R, power.t.test(n=ttest$n_E, delta = del_sig$delta_A,
                                           sd = del_sig$sigma, alternative = "two.sided", type = "two.sample")$power)
  }
  return(data.frame("ttestord"=results_me, "R"=results_R)) 
}

comp_ttest_k_1 <- comp_ttest_kieser(iter=10000)
mean(comp_ttest_k_1$ttestord -comp_ttest_k_1$R)



### Compare Sample Size Calculation
comp_ttest_samp <- function(alpha = 0.05, beta = 0.2, r = 1, iter = 1000, theta){
  results_R <-  c()
  results_me <-  c()
  for (i in seq_along(1:iter)) {
    set.seed(i)
    prob_length <- sample(c(3:14), 1)
    theta <-  runif(1, min = 0.1, max = 5)
    print(i)
    p <- generate_two_simplex_vectors(prob_length, log(theta))
    p_C <- p[[1]]
    p_E <- p[[2]]
    ttest <- samplesize_ttestord(p_C=p_C, p_E=p_E, alpha, beta, r)
    results_me <- rbind(results_me, ttest$n_total)
    
    del_sig <- calculate_delta_A_sigma(p_C, p_E)
    results_R <- c(results_R, ceiling(power.t.test(n=NULL, delta = del_sig$delta_A, power = 0.8, sig.level = 0.05,
      sd = del_sig$sigma, alternative = "two.sided", type = "two.sample")$n)*2)
  }
  return(data.frame("ttestord"=results_me, "R"=results_R)) 
}

comp_ttest_samp <- comp_ttest_samp(iter=10000)
mean(abs(comp_ttest_samp$ttestord-comp_ttest_samp$R)) # 0.0386

ggplot(data=comp_ttest_samp, aes(x=ttestord, y=R))+
  geom_point(alpha=0.6)+
  geom_abline(slope=1, color="red")+
  xlab("t-test")+
  ylab("base R")+
  ylim(0,100000)+
  xlim(0,100000)+
  theme_bw()
# what percentage of calculations are not equal
length(which(comp_ttest_samp$ttestord != comp_ttest_samp$R, arr.ind = TRUE))/nrow(comp_ttest_samp)
# whats the mean difference of the unequal pairs? -> 2
mean(comp_ttest_samp[which(comp_ttest_samp$ttestord != comp_ttest_samp$R, arr.ind = TRUE),1] -
       comp_ttest_samp[which(comp_ttest_samp$ttestord != comp_ttest_samp$R, arr.ind = TRUE),2])

