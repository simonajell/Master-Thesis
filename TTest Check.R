### Compare ttest function from R with my implementation
comp_ttest <- function(alpha = 0.05, beta = 0.2, r = 1, iter = 1000, bias = 1.2){
  results_R <-  c()
  results_me <-  c()
  for (i in seq_along(1:iter)) {
    set.seed(i)
    prob_length <- sample(c(3,4,5,6,7,8), 1)
    print(i)
    p <- generate_two_simplex_vectors(prob_length, bias)
    p_C <- p[[1]]
    p_E <- p[[2]]
    ttest <- samplesize_ttestord(p_C=p_C, p_E=p_E, alpha, beta, r)
    results_me <- rbind(results_me, ttest$actual_power)
    
    del_sig <- calculate_delta_A_sigma(p_C, p_E)
    results_R <- c(results_R, power.t.test(n=ttest$n_E, delta = del_sig$delta_A,
                 sd = del_sig$sigma, alternative = "two.sided", type = "two.sample")$power)
        }
 return(data.frame("ttestord"=results_me, "R"=results_R)) 
}
comp_ttest1 <- comp_ttest(iter = 1000)
mean(comp_ttest1$ttestord -comp_ttest1$R)
ggplot(data=comp_ttest1, aes(x=ttestord, y=R))+
  geom_point()+
  geom_abline(slope=1, color="red")

comp_ttest2 <- comp_ttest(iter = 10000)
mean(comp_ttest2$ttestord-comp_ttest2$R)
ggplot(data=comp_ttest2, aes(x=ttestord, y=R))+
  geom_point()+
  geom_abline(slope=1, color="red")+
  theme_bw()
length(which(comp_ttest2$R < comp_ttest2$ttestord))/10000


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
