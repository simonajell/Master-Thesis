 
# Function to calculate p_E when p_C and theta_A are known
# p_C is the known control group probabilities known from earlier trials etc
# theta_A is the influence of the treatment (control vs. treatment)(how much better/worse is the treatment)
# with those two one can calculate the probabilities of the experimental group, 
# which is the control group prob. multiplied with the theta (so how much better is the experimental group)
calc_p_E <- function(p_C, theta_A){
  p_e_1 <- rep(0, length(p_C))
  for (i in seq_along(1:length(p_C))) {
    p_e_1[i] <- sum(p_C[1:i])/(sum(p_C[1:i])+((1-sum(p_C[1:i]))*exp(-theta_A)))
  }
  p_e_2 <- rep(0, length(p_C))
  for (i in (2:length(p_C))) {
    p_e_2[1]<- p_e_1[1]
    p_e_2[i] <- (p_e_1[i]-p_e_1[(i-1)])
  }
  return(p_e_2)
}

####### Theta: assuming the PO assumption holds - trying around ####
# Kieser,2020 Example 4.3
# theta_A = ln(1.8)
p_c_po <- c(0.10, 0.20, 0.10, 0.15, 0.20, 0.25)
p_e_po <- c(0.167, 0.269, 0.110, 0.142, 0.156, 0.156)

new_p_e_po <- calc_p_E(p_c_po, theta_A = log(1.8))

calculate_theta_A(p_c_po, new_p_e_po, r=1)
log(1.8)


# Try to get rid of small difference, by simulating a large  integer-weighted sample
N <- 1e6
K <- 6
counts_C <- round(N * (r / (r + 1)) * p_c_po)
counts_E <- round(N * (1 / (r + 1)) * new_p_e_po)

outcome <- c(rep(1:K, counts_C), rep(1:K, counts_E))
group <- c(rep("control", sum(counts_C)), rep("treatment", sum(counts_E)))

data_sim <- data.frame(
  outcome = factor(outcome, ordered = TRUE),
  group = factor(group, levels = c("control", "treatment"))
)
model_sim <- polr(outcome ~ group, data = data_sim, method = "logistic")
model_sim$coefficients


calculate_theta_A2 <- function(p_C, p_E, r = 1) {
  # Validate inputs
  if (length(p_C) != length(p_E)) {
    stop("p_C and p_E must have the same length (same number of outcome levels).")
  }
  if (abs(sum(p_C) - 1) > 1e-6 || abs(sum(p_E) - 1) > 1e-6) {
    stop("Probabilities p_C and p_E must each sum to 1.")
  }
  N <- 1e4
  K <- length(p_C)  # number of outcome levels
  
  counts_C <- round(N * (r / (r + 1)) * p_C)
  counts_E <- round(N * (1 / (r + 1)) * p_E)
  
  # Construct outcome levels and group labels
  outcome <- c(rep(1:K, counts_C), rep(1:K, counts_E))
  group <- c(rep("control", sum(counts_C)), rep("treatment", sum(counts_E)))
  
  # Create synthetic dataset
  data_alt <- data.frame(
    outcome = factor(outcome, ordered = TRUE),
    group = factor(group, levels = c("control", "treatment"))
  )
  
  # Fit proportional odds model under the alternative
  model_alt <- polr(outcome ~ group, data = data_alt, method = "logistic", Hess = TRUE)
  
  # Extract coefficient of treatment effect
  coef_alt <- -coef(summary(model_alt))["grouptreatment", "Value"]
  var_alt <- coef(summary(model_alt))["grouptreatment", "Std. Error"]
  alt <- data.frame("coef" = coef_alt, "Var" = var_alt)
  return(alt)
} # langsame präzise methode

calculate_theta_A(p_c_po, new_p_e_po, r=1)
log(1.8)
calculate_theta_A2(p_c_po, new_p_e_po, r=1)

# a maybe faster way, not as robust

calculate_theta_A3 <- function(p_C, p_E, r = 1) {
  # Validate inputs
  if (length(p_C) != length(p_E)) {
    stop("p_C and p_E must have the same length (same number of outcome levels).")
  }
  if (abs(sum(p_C) - 1) > 1e-6 || abs(sum(p_E) - 1) > 1e-6) {
    stop("Probabilities p_C and p_E must each sum to 1.")
  }
  N <- 1e6
  K <- length(p_C)  # number of outcome levels
  
  counts_C <- round(N * (r / (r + 1)) * p_C)
  counts_E <- round(N * (1 / (r + 1)) * p_E)
  
  # Construct outcome levels and group labels
  outcome <- rep(1:K, 2) 
  group <- c(rep("control", K), rep("treatment", K))
  weight <- c(counts_C, counts_E)
  
  # Create synthetic dataset
  data_alt <- data.frame(
    outcome = factor(outcome, ordered = TRUE),
    group = factor(group, levels = c("control", "treatment")),
    weight = weight
  )
  
  # Fit proportional odds model under the alternative
  model_alt <- polr(outcome ~ group, data = data_alt, weights = weight, method = "logistic", Hess = TRUE)
  
  # Extract coefficient of treatment effect
  coef_alt <- -coef(summary(model_alt))["grouptreatment", "Value"]
  var_alt <- coef(summary(model_alt))["grouptreatment", "Std. Error"]
  alt <- data.frame("coef" = coef_alt, "Var" = var_alt)
  return(alt)
} # schnelle präzise methode
log(1.8)
calculate_theta_A(p_c_po, new_p_e_po, r=1)
calculate_theta_A2(p_c_po, new_p_e_po, r=1)
calculate_theta_A3(p_c_po, new_p_e_po, r=1)


### Function that does this for different probability vectors. and random thetas ###
random_simplex <- function(n) {
  success <- FALSE
  while (!success) {
    x <- rgamma(n, shape = 1, rate = 1)  # oder shape = alpha für andere Verteilungen
    x <- x / sum(x)
    # check for success
    success <- all(x > 0.0001)
  }
  return(x)
}


####### Actual Function to do the estimation check ######
theta_comparison <- function(comp_iter, r){
  theta_results <- data.frame("theta_A" = c(NA), "est_theta"=c(NA), "diff" = c(NA))
  p_c <- 0
  p_e <- 0
  for (i in seq_along(1:comp_iter)) {
    set.seed(i)
    print(i)
    vec_length <- sample(c(3:14), 1)
    theta_random <- runif(1, min = 1, max = 5)
    p_c <- random_simplex(n=vec_length)
    p_e <- calc_p_E(p_C = p_c, theta_A = log(theta_random))
    theta_results[i,1] <- log(theta_random)
    # White method
    theta_results[i,2] <- calculate_theta_A(p_c, p_e, r)[1,1]
  }
  theta_results$diff <- abs(theta_results$theta_A - theta_results$est_theta)
  return(theta_results)
}

t_comp1.2 <- theta_comparison(comp_iter = 10000, r=1)
mean(t_comp1.2$diff) # 0.00053
median(t_comp1.2$diff) # 0.00031
plot(t_comp1.2$theta_A, t_comp1.2$diff) # absolute difference
plot(t_comp1.2$theta_A, t_comp1.2$theta_A-t_comp1.2$est_theta) # difference
ggplot(data=t_comp1.2, aes(x=theta_A, y=est_theta))+
  geom_point(alpha=0.6)+
  geom_abline(slope=1, color="red")+
  xlab("log odds ratio")+
  ylab("estimated theta")+
  theme_bw()


#for r=0.8
t_comp2 <- theta_comparison(comp_iter = 10000, vec_length = 6, r=0.8)
mean(t_comp2$diff) # 0.00054
median(t_comp2$diff) # 0.00034
plot(t_comp2$theta_A, t_comp2$diff) 
plot(t_comp2$theta_A, t_comp2$theta_A-t_comp2$est_theta) 


# for r=0.4
t_comp3 <- theta_comparison(comp_iter = 10000, vec_length = 6, r=0.4)
mean(t_comp3$diff) # 0.00066
median(t_comp3$diff) # 0.00042
plot(t_comp3$theta_A, t_comp3$diff) 
plot(t_comp3$theta_A, t_comp3$theta_A-t_comp3$est_theta) 


# for r=0.2
t_comp8 <- theta_comparison(comp_iter = 10000, vec_length = 6, r=0.2)
mean(t_comp8$diff) # 0.00095
median(t_comp8$diff) # 0.00065
plot(t_comp8$theta_A, t_comp8$diff) 
plot(t_comp8$theta_A, t_comp8$theta_A-t_comp8$est_theta) 


# for vec_length = 4
t_comp4 <- theta_comparison(comp_iter = 10000, vec_length = 4, r=1)
mean(t_comp4$diff) # 0.0005
plot(t_comp4$theta_A, t_comp4$diff) 

# for vec_length = 3
t_comp5 <- theta_comparison(comp_iter = 10000, vec_length = 3, r=1)
mean(t_comp5$diff) # 0.00043
plot(t_comp5$theta_A, t_comp5$diff) 

# for vec_length = 8
t_comp6 <- theta_comparison(comp_iter = 10000, vec_length = 8, r=1) 
mean(t_comp6$diff) # 0.00053
plot(t_comp6$theta_A, t_comp6$diff) 


# with new technique - this is old
theta_comparison_new <- function(comp_iter, vec_length, r){
  theta_results <- data.frame("theta_A" = c(NA), "est_theta"=c(NA), "diff" = c(NA))
  p_c <- 0
  p_e <- 0
  for (i in seq_along(1:comp_iter)) {
    set.seed(i)
    print(i)
    theta_random <- runif(1, min = 1, max = 5)
    p_c <- random_simplex(n=vec_length)
    p_e <- calc_p_E(p_C = p_c, theta_A = log(theta_random))
    theta_results[i,1] <- log(theta_random)
    # White method
    theta_results[i,2] <- calculate_theta_A(p_c, p_e, r)[1,1]
  }
  theta_results$diff <- abs(theta_results$theta_A - theta_results$est_theta)
  return(theta_results)
}

t_comp_new <- theta_comparison_new(comp_iter = 10000, vec_length = 6, r=1)
mean(t_comp_new$diff) # 7.10 e-05
median(t_comp_new$diff) # 1.02 e-05
plot(t_comp_new$theta_A, t_comp_new$diff) 
plot(t_comp_new$theta_A, t_comp_new$theta_A - t_comp_new$est_theta) 



####### Variance Calculation for PO vectors ######

# Compare Kieser variance to White variance:
# White:
var_w <- calculate_theta_N(p_c_po, p_e_po, r)[1,2]^2
calculate_theta_N(p_c_po, p_e_po, r)[1,2]^2
#Kieser:
x = 0
for (i in 1:length(p_c_po)){
  x = x + ((p_c_po[i] + r*p_e_po[i]) / (1 + r))^3
}
var_k <- 12/(1-x)


# make function, that generates multiple different probability vectors that fulfill PO assumption and calculate both variances
var_comparison <- function(comp_iter, vec_length, r){
  var_results <- data.frame("theta_A" = c(NA), "var_w"=c(NA), "var_k" = c(NA), "diff" = c(NA))
  for (i in seq_along(1:comp_iter)) {
    set.seed(i)
    vec_length <- sample(3:14, 1)
    theta_random <- runif(1, min = 1, max = 5)
    p_c_po <- random_simplex(vec_length)
    p_e_po <- calc_p_E(p_c_po, theta_A = log(theta_random))
    var_results[i,1] <- theta_random
    # White method
    var_results[i,2] <- calculate_theta_N(p_c_po, p_e_po, r)[1,2]^2
    
    # Kieser method
    x = 0
    for (j in 1:vec_length){
      x = x + ((r*p_c_po[j] + p_e_po[j]) / (1 + r))^3
    }
    var_results[i,3] <- ((1+r)^2/r)* (3/(1-x))
  }
  var_results$diff <- abs(var_results$var_w - var_results$var_k)
  return(var_results)
}

v_comp1 <- var_comparison(comp_iter = 10000, vec_length = 6, r=1)
mean(v_comp1$diff) # on average a difference of 0.00078
median(v_comp1$diff) # 0.00031
plot(v_comp1$theta_A, v_comp1$diff) # difference between variances gets bigger with bigger thetas
plot(v_comp1$var_w,  abs(v_comp1$var_w - v_comp1$var_k))

ggplot(data=v_comp1, aes(x=var_w, y=var_k))+
  geom_point(alpha=0.6)+
  geom_abline(slope=1, color="red")+
  xlab("V_N")+
  ylab("V_Kieser")+
  theme_bw()

# for vec_length = 4
v_comp2 <- var_comparison(comp_iter = 1000, vec_length = 4, r=1)
mean(v_comp2$diff) 
plot(v_comp2$theta_A, v_comp2$diff) 

####### Variance Calculation for non PO vectors ######
generate_two_simplex_vectors <- function(n, bias_strength = 2) {
  # Vektor A: gleichmäßige Dirichlet-Verteilung
  alpha_a <- rep(1, n)
  a <- rgamma(n, shape = alpha_a, rate = 1)
  a <- a / sum(a)
  
  # Vektor B: größere alpha-Werte an denselben Stellen, wo a hoch ist
  # Das verstärkt die hohen Werte weiter
  alpha_b <- 1 + bias_strength * a  # bias_strength steuert die Verstärkung
  b <- rgamma(n, shape = alpha_b, rate = 1)
  b <- b / sum(b)
  
  list(a = a, b = b)
}

var_comparison_no_PO <- function(comp_iter, vec_length, r, theta_A, bias){
  var_results <- data.frame("theta_A" = c(NA), "var_w"=c(NA), "var_k" = c(NA), "diff" = c(NA))
  for (i in seq_along(1:comp_iter)) {
    set.seed(i)
    theta_random <- runif(1, min = 1, max = 3)
    p <- generate_two_simplex_vectors(vec_length, bias_strength = bias)
    p_c_po <- p[[1]]
    p_e_po <- p[[2]]
    var_results[i,1] <- theta_random
    # White method
    var_results[i,2] <- calculate_theta_N(p_c_po, p_e_po, r)[1,2]^2
    
    # Kieser method
    x = 0
    for (j in 1:vec_length){
      x = x + ((r*p_c_po[j] + p_e_po[j]) / (1 + r))^3
    }
    var_results[i,3] <- ((1+r)^2/r)* (3/(1-x))
  }
  var_results$diff <- var_results$var_w - var_results$var_k
  return(var_results)
}

v_comp_no_PO <- var_comparison_no_PO(comp_iter = 1000, vec_length = 6, r=1, bias = 1)
mean(v_comp_no_PO$diff) # on average a difference of -6.056277e-05
plot(v_comp_no_PO$theta_A, v_comp_no_PO$diff) # difference between variances gets bigger with bigger thetas

# two independent probability vectors
var_comparison_no_PO2 <- function(comp_iter, vec_length, r, theta_A){
  var_results <- data.frame("theta_A" = c(NA), "var_w"=c(NA), "var_k" = c(NA), "diff" = c(NA))
  for (i in seq_along(1:comp_iter)) {
    set.seed(i)
    theta_random <- runif(1, min = 1, max = 3)
    p_c_po <- random_simplex(vec_length)
    p_e_po <- random_simplex(vec_length)
    var_results[i,1] <- theta_random
    # White method
    var_results[i,2] <- calculate_theta_N(p_c_po, p_e_po, r)[1,2]^2
    
    # Kieser method
    x = 0
    for (j in 1:vec_length){
      x = x + ((r*p_c_po[j] + p_e_po[j]) / (1 + r))^3
    }
    var_results[i,3] <- ((1+r)^2/r)* (3/(1-x))
  }
  var_results$diff <- var_results$var_w - var_results$var_k
  return(var_results)
}
v_comp_no_PO2 <- var_comparison_no_PO2(comp_iter = 1000, vec_length = 6, r=1)
mean(v_comp_no_PO2$diff) # on average a difference of -6.589431e-05
plot(v_comp_no_PO2$theta_A, v_comp_no_PO2$diff) 


#### Comparison Plot of PO vs no-PO Varianz Calculation
var_comp_noPO_PO <- data.frame("diff" = c(v_comp1$diff, v_comp_no_PO2$diff), 
                               "theta" = c(v_comp1$theta_A, v_comp_no_PO2$theta_A),
                               "PO"= c(rep("PO", length(v_comp1$diff)), rep("no_PO", length(v_comp_no_PO$diff))))
ggplot(var_comp_noPO_PO, aes(x=theta, y=diff, color=PO))+
  geom_point(shape = 1)+
  theme_bw()



####### Sample Size Calculation when PO fulfilled #####
p_C_PO <- random_simplex(5)
p_E_PO <- calc_p_E(p_C_PO, theta_A = log(1.8))

exp(calculate_theta_A(p_C_PO, p_E_PO))
samplesize_po_NN(p_C_PO, p_E_PO, alpha = 0.05, beta = 0.2, r =1)
samplesize_po_kieser(p_C_PO, p_E_PO, alpha = 0.05, beta = 0.2, r =1)
samplesize_po_kieser_known(p_C_PO, p_E_PO, theta = 1.8, alpha = 0.05, beta = 0.2, r =1)


comp_PO_fullfilled <- function(alpha=0.05, beta=0.2, r=1, iter=1000){
  results_NN <-  data.frame()
  results_AA <-  data.frame()
  results_NA <-  data.frame()
  results_k <- data.frame()
  results_k_known <- data.frame()
  theta_vec <- c(NA)
  for (i in seq_along(1:iter)) {
    set.seed(i)
    prob_length <- sample(3:14, 1)
    theta_A <- runif(1, min = 1.05, max = 5)
    print(i)
    p_C <- random_simplex(prob_length)
    p_E <- calc_p_E(p_C, theta_A = log(theta_A))
    theta_vec[i] <- calculate_theta_A(p_C, p_E)
    results_NN <- rbind(results_NN, as.data.frame(samplesize_po_NN(p_C=p_C, p_E=p_E, alpha, beta, r)))
    results_AA <- rbind(results_AA, as.data.frame(samplesize_po_AA(p_C=p_C, p_E=p_E, alpha, beta, r)))
    results_NA <- rbind(results_NA, as.data.frame(samplesize_po_NA(p_C=p_C, p_E=p_E, alpha, beta, r)))
    results_k <- rbind(results_k, as.data.frame(samplesize_po_kieser(p_C=p_C, p_E=p_E, alpha, beta, r)))
    results_k_known <- rbind(results_k_known, 
                             as.data.frame(samplesize_po_kieser_known(p_C=p_C, p_E=p_E, theta = theta_A, alpha, beta, r)))
    
  }
  return(list(results_NN, results_k, results_AA, results_NA, results_k_known, theta_vec))
}

Samp_Size_PO_fulfilled <- comp_PO_fullfilled(iter = 10000)
# how often is NN the same as Kieser
length(which(Samp_Size_PO_fulfilled[[1]]$n_total == Samp_Size_PO_fulfilled[[2]]$n_total))/10000 # 98.51% , 149 unequal
mean(abs(Samp_Size_PO_fulfilled[[1]]$n_total - Samp_Size_PO_fulfilled[[2]]$n_total))

mean(abs(Samp_Size_PO_fulfilled[[1]]$n_total[which(Samp_Size_PO_fulfilled[[1]]$n_total != Samp_Size_PO_fulfilled[[2]]$n_total, arr.ind = TRUE)] -
       Samp_Size_PO_fulfilled[[2]]$n_total[which(Samp_Size_PO_fulfilled[[1]]$n_total != Samp_Size_PO_fulfilled[[2]]$n_total, arr.ind = TRUE)]))

length(which(Samp_Size_PO_fulfilled[[1]]$n_total!=Samp_Size_PO_fulfilled[[2]]$n_total, arr.ind = FALSE))


# how often is AA the same as Kieser
length(which(Samp_Size_PO_fulfilled[[2]]$n_total == Samp_Size_PO_fulfilled[[3]]$n_total))/10000 # 0%
mean(Samp_Size_PO_fulfilled[[3]]$n_total - Samp_Size_PO_fulfilled[[2]]$n_total)
# how often is NA the same as Kieser
length(which(Samp_Size_PO_fulfilled[[2]]$n_total == Samp_Size_PO_fulfilled[[4]]$n_total))/10000 # 0.73 %
mean(Samp_Size_PO_fulfilled[[4]]$n_total - Samp_Size_PO_fulfilled[[2]]$n_total)

# how often is NN the same as Kieser with true theta
length(which(Samp_Size_PO_fulfilled[[1]]$n_total == Samp_Size_PO_fulfilled[[5]]$n_total))/10000 # 78.22%
# how often is kieser the same as Kieser with true theta
length(which(Samp_Size_PO_fulfilled[[2]]$n_total == Samp_Size_PO_fulfilled[[5]]$n_total))/10000 # 78.25%

mean(Samp_Size_PO_fulfilled[[1]]$n_total - Samp_Size_PO_fulfilled[[2]]$n_total)

# vs. not fulfilled
p <- generate_two_simplex_vectors(5, bias_strength = 1.8)
p_C_no_PO <- p[[1]]
p_E_no_PO <- p[[2]]

samplesize_po_NN(p_C_no_PO, p_E_no_PO, alpha = 0.05, beta = 0.2, r =1)


####### Sample Size Calculation when PO fulfilled and categories are normally distributed #####
comp_PO_fullfilled <- function(alpha=0.05, beta=0.2, r=1, iter=1000, theta_A = 1.8){
  results_NN <-  data.frame()
  results_NN <-  data.frame()
  theta_vec <- c(NA)
  for (i in seq_along(1:iter)) {
    set.seed(i)
    prob_length <- sample(c(3,4,5,6,7,8), 1)
    print(i)
    p_C <- random_simplex(prob_length)
    p_E <- calc_p_E(p_C, theta_A = log(theta_A))
    theta_vec[i] <- calculate_theta_A(p_C, p_E)
    results_NN <- rbind(results_NN, as.data.frame(samplesize_po_NN(p_C=p_C, p_E=p_E, alpha, beta, r)))
    results_AA <- rbind(results_AA, as.data.frame(samplesize_po_AA(p_C=p_C, p_E=p_E, alpha, beta, r)))
    results_NA <- rbind(results_NA, as.data.frame(samplesize_po_NA(p_C=p_C, p_E=p_E, alpha, beta, r)))
    results_k <- rbind(results_k, as.data.frame(samplesize_po_kieser(p_C=p_C, p_E=p_E, alpha, beta, r)))
    results_k_known <- rbind(results_k_known, 
                             as.data.frame(samplesize_po_kieser_known(p_C=p_C, p_E=p_E, theta = theta_A, alpha, beta, r)))
    
  }
  return(list(results_NN, results_k, results_AA, results_NA, results_k_known, theta_vec))
}