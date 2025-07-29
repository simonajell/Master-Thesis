## Try to reproduce ologit in R
library(MASS)
library(ggplot2)
library(tidyr)

p_C2 <- c(0.3,0.2,0.15,0.2, 0.1, 0.05)  # control group outcome probabilities
p_E2 <- c(0.15,0.3,0.1, 0.2, 0.2, 0.05)  # experimental group outcome probabilities
p_C <- c(0.2, 0.3, 0.3, 0.2)
p_E <- c(0.1, 0.2, 0.4, 0.3)

##### Function that generates a simplex vector ####
# Zufallsvektoren auf dem Simplex
random_simplex <- function(n) {
  x <- rgamma(5, shape = 1, rate = 1)  # oder shape = alpha für andere Verteilungen
  x / sum(x)
}

normal_simplex <- function(n, favored_cat) {
  x <- seq(1, n, by = 1)
  y <- dnorm(x, mean =  x[favored_cat], sd = 1.8)
  noise <- runif(n, min = 0.15, max = 0.3) # Add random noise (for example, 10% variation)
  y_noisy <- y * noise  # Apply noise and normalize
  y1 <- y_noisy/sum(y_noisy)
  y <- y/sum(y)
  plot(x,y)
  plot(x,y1)
  return(list(y,y1))
}


# Function that generates two simplex vectors, where the second has higher values on average
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

##### Function that generates a simplex vector with a uniform distribution ####
uniform_simplex <- function(n) {
  uni <- rep(1/n, n)
  noise <- runif(n, min = 0.7, max = 0.9) # Add random noise (for example, 10% variation)
  uni_noisy <- uni * noise  # Apply noise and normalize
  uni_noisy <- uni_noisy/sum(uni_noisy)
  return(list(uni, uni_noisy))
}
uniform_simplex(4)
############# Proportional Odds Method #########################################
# This function (adapted from Kieser 2020) will calculate the needed sample size based on p_C and p_E for proportional odds regression. A piece is still missing. Please ignore at this stage.
# calculate theta_A with artcat (Ian Whites) method
calculate_theta_A_old <- function(p_C, p_E, r = 1) {
  # Validate inputs
  if (length(p_C) != length(p_E)) {
    stop("p_C and p_E must have the same length (same number of outcome levels).")
  }
  if (abs(sum(p_C) - 1) > 1e-6 || abs(sum(p_E) - 1) > 1e-6) {
    stop("Probabilities p_C and p_E must each sum to 1.")
  }
  
  K <- length(p_C)  # number of outcome levels
  
  # Construct outcome levels and group labels
  outcome <- rep(c(1:K), times = 2)
  group <- c(rep(c("control"), times = K), rep(c("treatment"), times = K))
  
  # Calculate weights under the alternative hypothesis
  weights_alt <- c(((r * p_C) / (r + 1)), (p_E / (r + 1)))
  
  # Create synthetic dataset
  data_alt <- data.frame(
    outcome = factor(outcome, ordered = TRUE),
    group = factor(group, levels = c("control", "treatment")),
    weight = weights_alt
  )
  
  # Fit proportional odds model under the alternative
  model_alt <- suppressWarnings(
    polr(outcome ~ group, data = data_alt, weights = weight, method = "logistic", Hess = TRUE)
  )
  
  # Extract coefficient of treatment effect
  coef_alt <- -coef(summary(model_alt))["grouptreatment", "Value"]
  var_alt <- coef(summary(model_alt))["grouptreatment", "Std. Error"]
  alt <- data.frame("coef" = coef_alt, "Var" = var_alt)
  return(alt)
}

calculate_theta_A <- function(p_C, p_E, r = 1) {
  # Validate inputs
  if (length(p_C) != length(p_E)) {
    stop("p_C and p_E must have the same length (same number of outcome levels).")
  }
  if (abs(sum(p_C) - 1) > 1e-6 || abs(sum(p_E) - 1) > 1e-6) {
    stop("Probabilities p_C and p_E must each sum to 1.")
  }
  K <- length(p_C)  # number of outcome levels
  
  # Construct outcome levels and group labels
  outcome <- rep(1:K, 2) 
  group <- c(rep("control", K), rep("treatment", K))
  
  if(length(p_C)==3){
    weight <- c(((r * p_C) / (r + 1)), (p_E / (r + 1)))
  } else {
    N <- 1e6
    weight <- round(N*c(((r * p_C) / (r + 1)), (p_E / (r + 1))))
  }
  
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
}

calculate_theta_A2 <- function(p_C, p_E, r = 1) {
  # Validate inputs
  if (length(p_C) != length(p_E)) {
    stop("p_C and p_E must have the same length (same number of outcome levels).")
  }
  if (abs(sum(p_C) - 1) > 1e-6 || abs(sum(p_E) - 1) > 1e-6) {
    stop("Probabilities p_C and p_E must each sum to 1.")
  }
  if(length(p_C)==3){
    N <- 1
    counts_C <- (r / (r + 1)) * p_C
    counts_E <- (1 / (r + 1)) * p_E
    } else {
      N <- 1e6
      counts_C <- round(N * (r / (r + 1)) * p_C)
      counts_E <- round(N * (1 / (r + 1)) * p_E)
    }
  K <- length(p_C)  # number of outcome levels
  
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
}

calculate_theta_A(p_C = p_C, p_E=p_E, r=1)
calculate_theta_A(p_C=c(0.1,0.2,0.4,0.2, 0.1), p_E=c(0.15,0.35,0.1, 0.2, 0.2), r=1)
calculate_theta_A_old(p_C=c(0.1,0.2,0.4,0.2, 0.1), p_E=c(0.15,0.35,0.1, 0.2, 0.2), r=1)
calculate_theta_A2(p_C=c(0.1,0.2,0.4,0.2, 0.1), p_E=c(0.15,0.35,0.1, 0.2, 0.2), r=1)
calculate_theta_A(p_C=c(0.3,0.2,0.15,0.2, 0.1, 0.05), p_E=c(0.15,0.3,0.1, 0.2, 0.2, 0.05), r=1)
calculate_theta_A(p_C=random_simplex(3), p_E=random_simplex(3), r=1)

calculate_theta_N <- function(p_C, p_E, r = 1) {
  # Null Hypothesis states: theta=0
  # Validate inputs
  if (length(p_C) != length(p_E)) {
    stop("p_C and p_E must have the same length (same number of outcome levels).")
  }
  if (abs(sum(p_C) - 1) > 1e-6 || abs(sum(p_E) - 1) > 1e-6) {
    stop("Probabilities p_C and p_E must each sum to 1.")
  }
  
  K <- length(p_C)  # number of outcome levels
  
  p_i <- (r*p_C+p_E)/(r+1) # expected value of estimand (e.g. risk difference)
                            # make both probabilities equal, so that theta = 0
  
  # Construct outcome levels and group labels
  outcome <- rep(c(1:K), times = 2)
  group <- c(rep(c("control"), times = K), rep(c("treatment"), times = K))
  
  # Calculate weights under the alternative hypothesis
  weights_null <- c(((r * p_i) / (r + 1)), (p_i / (r + 1)))
  
  # Create synthetic dataset
  data_null <- data.frame(
    outcome = factor(outcome, ordered = TRUE),
    group = factor(group, levels = c("control", "treatment")),
    weight = weights_null
  )
  
  # Fit proportional odds model under the alternative
  model_null <- suppressWarnings(
    polr(outcome ~ group, data = data_null, weights = weight, method = "logistic", Hess = TRUE)
  )
  
  # Extract coefficient of treatment effect
  coef_null <- -coef(summary(model_null))["grouptreatment", "Value"]
  var_null <- coef(summary(model_null))["grouptreatment", "Std. Error"]
  null <- data.frame("coef" = coef_null, "Var" = var_null)
  return(null)
}

calculate_theta_N(p_C = p_C, p_E=p_E, r=1)
calculate_theta_N(p_C=c(0.1,0.2,0.4,0.2, 0.1), p_E=c(0.15,0.35,0.1, 0.2, 0.2), r=1)
calculate_theta_N(p_C=c(0.3,0.2,0.15,0.2, 0.1, 0.05), p_E=c(0.15,0.3,0.1, 0.2, 0.2, 0.05), r=1)


samplesize_po_AA <- function(p_C, p_E, alpha, beta, r, p_C2=NULL, p_E2=NULL){
  
  theta_A <- calculate_theta_A_old(p_C, p_E, r)[1,1]
  Var_A <- calculate_theta_A_old(p_C, p_E, r)[1,2]^2
  
  n <- (Var_A*((qnorm(1-alpha/2)+qnorm(1-beta))^2))/(theta_A^2) # White, 2023 Formel 5 (AA)
  n_E <- ceiling(r/(1+r) * n)
  n_C <- ceiling(1/(1+r) * n)
  n_total <- n_E + n_C
  actual_r <- n_E/n_C
  actual_power <- pnorm(sqrt((n_total*(theta_A^2))/Var_A)-qnorm(1-alpha/2))
  actual_power2<-NA
  if (!is.null(p_C2)&!is.null(p_E2)){
    theta_A2 <- calculate_theta_A(p_C2, p_E2, r)[1,1]
    Var_A2 <- calculate_theta_A(p_C2, p_E2, r)[1,2]^2
    actual_power2 <- pnorm(sqrt((n_total*(theta_A2^2))/Var_A2)-qnorm(1-alpha/2))
   }
  return(data.frame(n_total = n_total, n_E = n_E, n_C = n_C,
                    actual_r = actual_r, actual_power = actual_power, actual_power2=actual_power2))
}

samplesize_po_NN <- function(p_C, p_E, alpha, beta, r, p_C2=NULL, p_E2=NULL){
  
  theta_A <- calculate_theta_A_old(p_C, p_E, r)[1,1]
  Var_N <- calculate_theta_N(p_C, p_E, r)[1,2]^2
  
  n <- (Var_N*((qnorm(1-alpha/2)+qnorm(1-beta))^2))/(theta_A^2) # White, 2023 Formel 5 (AA)
  n_E <- ceiling(r/(1+r) * n)
  n_C <- ceiling(1/(1+r) * n)
  n_total <- n_E + n_C
  actual_r <- n_E/n_C
  actual_power <- pnorm(sqrt((n_total*(theta_A^2))/Var_N)-qnorm(1-alpha/2))
  actual_power2<-NA
  if (!is.null(p_C2)&!is.null(p_E2)){
    theta_A2 <- calculate_theta_A_old(p_C2, p_E2, r)[1,1]
    Var_N2 <- calculate_theta_N(p_C2, p_E2, r)[1,2]^2
    actual_power2 <- pnorm(sqrt((n_total*(theta_A2^2))/Var_N2)-qnorm(1-alpha/2))
  }
  return(data.frame(n_total = n_total, n_E = n_E, n_C = n_C,
                    actual_r = actual_r, actual_power = actual_power, actual_power2=actual_power2))
}

samplesize_po_NA <- function(p_C, p_E, alpha, beta, r, p_C2=NULL, p_E2=NULL){
  
  d <- calculate_theta_A_old(p_C, p_E, r)[1,1]
  V_A <- calculate_theta_A_old(p_C, p_E, r)[1,2]^2
  V_N <- calculate_theta_N(p_C, p_E, r)[1,2]^2

  n <- ((sqrt(V_N)*qnorm(1-alpha/2)+sqrt(V_A)*qnorm(1-beta))^2)/(d^2) # White, 2023 Formel 4 (NA)
  n_E <- ceiling(r/(1+r) * n)
  n_C <- ceiling(1/(1+r) * n)
  n_total <- n_E + n_C
  actual_r <- n_E/n_C
  actual_power <- pnorm(((sqrt(n_total*d^2))-(sqrt(V_N)*qnorm(1-alpha/2)))/(sqrt(V_A)))
  actual_power2<-NA
  if (!is.null(p_C2)&!is.null(p_E2)){
    d2 <- calculate_theta_A_old(p_C2, p_E2, r)[1,1]
    V_A2 <- calculate_theta_A_old(p_C2, p_E2, r)[1,2]^2
    V_N2 <- calculate_theta_N(p_C2, p_E2, r)[1,2]^2
    actual_power2 <- pnorm(((sqrt(n_total*d2^2))-(sqrt(V_N2)*qnorm(1-alpha/2)))/(sqrt(V_A2)))
  }
  return(data.frame(n_total = n_total, n_E = n_E, n_C = n_C,
                    actual_r = actual_r, actual_power = actual_power, actual_power2=actual_power2))
}

samplesize_po_kieser <- function(p_C, p_E, alpha, beta, r, p_C2=NULL, p_E2=NULL){
  
  theta_A <- calculate_theta_A_old(p_C, p_E, r)[1,1]
  x = 0
  for (i in 1:length(p_C)){
    x = x + ((p_C[i] + r*p_E[i]) / (1 + r))^3
  }
    n <- (1+r)^2/r * 3*(qnorm(1-alpha/2) + qnorm(1-beta))^2 /
      (theta_A^2 * (1-x)) # Kieser Formel 4.8
  n_E <- ceiling(r/(1+r) * n)
  n_C <- ceiling(1/(1+r) * n)
  n_total <- n_E + n_C
  actual_r <- n_E/n_C
  actual_power <- pnorm(sqrt(n_total * actual_r * theta_A^2 * (1-x)/
                               (3 * (1 + actual_r)^2)) - qnorm(1-alpha/2))
  actual_power2<-NA
  if (!is.null(p_C2)&!is.null(p_E2)){
    theta_A2 <- calculate_theta_A(p_C2, p_E2, r)[1,1]
    actual_power2 <- pnorm(sqrt(n_total * actual_r * theta_A2^2 * (1-x)/
                                  (3 * (1 + actual_r)^2)) - qnorm(1-alpha/2))
  }
  return(data.frame(n_total = n_total, n_E = n_E, n_C = n_C,
                    actual_r = actual_r, actual_power = actual_power, actual_power2=actual_power2))
}
samplesize_po_kieser_known <- function(p_C, p_E, theta, alpha, beta, r, p_C2=NULL, p_E2=NULL){
  
  theta_A <- log(theta)
  x = 0
  for (i in 1:length(p_C)){
    x = x + ((p_C[i] + r*p_E[i]) / (1 + r))^3
  }
  n <- ((1+r)^2/r) * 3*((qnorm(1-alpha/2) + qnorm(1-beta))^2) /
    (theta_A^2 * (1-x)) # Kieser Formel 4.8
  n_E <- ceiling(r/(1+r) * n)
  n_C <- ceiling(1/(1+r) * n)
  n_total <- n_E + n_C
  actual_r <- n_E/n_C
  actual_power <- pnorm(sqrt(n_total * actual_r * theta_A^2 * (1-x)/
                               (3 * (1 + actual_r)^2)) - qnorm(1-alpha/2))
  actual_power2<-NA
  if (!is.null(p_C2)&!is.null(p_E2)){
    theta_A2 <- calculate_theta_A(p_C2, p_E2, r)[1,1]
    actual_power2 <- pnorm(sqrt(n_total * actual_r * theta_A2^2 * (1-x)/
                                  (3 * (1 + actual_r)^2)) - qnorm(1-alpha/2))
  }
  return(data.frame(n_total = n_total, n_E = n_E, n_C = n_C,
                    actual_r = actual_r, actual_power = actual_power, actual_power2=actual_power2))
}

samplesize_po_AA(p_C=p_C, p_E=p_E, alpha = 0.05, beta = 0.2, r=1)
samplesize_po_AA(p_C=p_C, p_E=p_E, alpha = 0.05, beta = 0.2, r=1, p_C2, p_E2)
samplesize_po_AA(p_C=p_C, p_E=p_E, alpha = 0.05, beta = 0.2, r=1, p_C2= c(0.1,0.2,0.4,0.3), p_E2=c(0.25,0.35,0.2, 0.2))
samplesize_po_AA(p_C=c(0.1,0.2,0.4,0.3), p_E=c(0.15,0.35,0.2, 0.3), alpha = 0.05, beta = 0.2, r=1)

t1<- samplesize_po_NN(p_C=p_C, p_E=p_E, alpha = 0.05, beta = 0.2, r=1)
t2<-samplesize_po_NN(p_C=p_C, p_E=p_E, alpha = 0.05, beta = 0.2, r=1, p_C2, p_E2)
t3<-samplesize_po_NN(p_C=p_C, p_E=p_E, alpha = 0.05, beta = 0.2, r=1, p_C2= c(0.1,0.2,0.4,0.3), p_E2=c(0.25,0.35,0.2, 0.2))
samplesize_po_NN(p_C=c(0.1,0.2,0.4,0.3), p_E=c(0.15,0.35,0.2, 0.3), alpha = 0.05, beta = 0.2, r=1)



samplesize_po_NA(p_C=p_C, p_E=p_E, alpha = 0.05, beta = 0.2, r=1)
samplesize_po_NA(p_C=c(0.1,0.2,0.4,0.3), p_E=c(0.15,0.35,0.2, 0.3), alpha = 0.05, beta = 0.2, r=1)

samplesize_po_kieser(p_C=p_C, p_E=p_E, alpha = 0.05, beta = 0.2, r=1)
samplesize_po_kieser(p_C=c(0.1,0.2,0.4,0.3), p_E=c(0.15,0.35,0.2, 0.3), alpha = 0.05, beta = 0.2, r=1)

##### Compare Kieser to White Method #####
comp_PO <- function(alpha=0.05, beta=0.2, r=1, iter=1000, bias = 1.2){
  results_NN <-  data.frame()
  results_AA <-  data.frame()
  results_NA <-  data.frame()
  results_k <- data.frame()
  for (i in seq_along(1:iter)) {
  set.seed(i)
  prob_length <- sample(c(3,4,5,6,7,8), 1)
  print(i)
  p <- generate_two_simplex_vectors(prob_length, bias)
  p_C <- p[[1]]
  p_E <- p[[2]]
  results_NN <- rbind(results_NN, as.data.frame(samplesize_po_NN(p_C=p_C, p_E=p_E, alpha, beta, r)))
  results_AA <- rbind(results_AA, as.data.frame(samplesize_po_AA(p_C=p_C, p_E=p_E, alpha, beta, r)))
  results_NA <- rbind(results_NA, as.data.frame(samplesize_po_NA(p_C=p_C, p_E=p_E, alpha, beta, r)))
  results_k <- rbind(results_k, as.data.frame(samplesize_po_kieser(p_C=p_C, p_E=p_E, alpha, beta, r)))
  }
  return(list(results_NN, results_k, results_AA, results_NA))
}

comp_w_k <- comp_PO(iter=10000)
# how often is NN the same as Kieser
length(which(comp_w_k[[1]]$n_total == comp_w_k[[2]]$n_total))/10000 # 94.8% , 520 unequal
mean(comp_w_k[[1]]$n_total - comp_w_k[[2]]$n_total)
# how often is AA the same as Kieser
length(which(comp_w_k[[2]]$n_total == comp_w_k[[3]]$n_total))/10000 # 0%
mean(comp_w_k[[2]]$n_total - comp_w_k[[3]]$n_total)
# how often is NA the same as Kieser
length(which(comp_w_k[[2]]$n_total == comp_w_k[[4]]$n_total))/10000 # 0.53 %
mean(comp_w_k[[4]]$n_total - comp_w_k[[2]]$n_total)

# look at Power:
mean(comp_w_k[[1]]$actual_power - comp_w_k[[2]]$actual_power)

# what is the mean difference of the unequal sample sizes of Kieser and NN?
neq <- which(comp_w_k[[1]]$n_total != comp_w_k[[2]]$n_total)
mean(abs(comp_w_k[[1]]$n_total[neq] - comp_w_k[[2]]$n_total[neq])) # mean abs difference of unequal sample sizes
mean(abs(comp_w_k[[1]]$n_total - comp_w_k[[2]]$n_total))# mean abs difference: 44.85
median(comp_w_k[[1]]$n_total - comp_w_k[[2]]$n_total)# median abs difference = 0
mean(comp_w_k[[1]]$n_total[neq]) # mean NN sample size of unequal pairs


# Bias = 2
comp_w_k2 <- comp_PO(iter=10000, bias = 2)
length(which(comp_w_k2[[1]]$n_total == comp_w_k2[[2]]$n_total))/10000 # 94.04%
neq <- which(comp_w_k2[[1]]$n_total != comp_w_k2[[2]]$n_total)
mean(comp_w_k2[[1]]$n_total[neq] - comp_w_k2[[2]]$n_total[neq])
plot(comp_w_k2[[1]]$n_total, comp_w_k2[[2]]$n_total)
# look at Power:
mean(comp_w_k2[[1]]$actual_power - comp_w_k2[[2]]$actual_power)

# Bias = 3
comp_w_k3 <- comp_PO(iter=10000, bias = 3)
length(which(comp_w_k3[[1]]$n_total == comp_w_k3[[2]]$n_total))/10000 # 93.93%
neq <- which(comp_w_k3[[1]]$n_total != comp_w_k3[[2]]$n_total)
mean(comp_w_k3[[1]]$n_total[neq] - comp_w_k3[[2]]$n_total[neq])
# look at Power:
mean(comp_w_k3[[1]]$actual_power - comp_w_k3[[2]]$actual_power)

# r = 0.8
comp_w_k4 <- comp_PO(iter=10000, r=0.8)
length(which(comp_w_k4[[1]]$n_total == comp_w_k4[[2]]$n_total))/10000 # 31.82%
# how often is AA the same as Kieser
length(which(comp_w_k4[[2]]$n_total == comp_w_k4[[3]]$n_total))/10000 # 0.88%
# how often is NA the same as Kieser
length(which(comp_w_k4[[2]]$n_total == comp_w_k4[[4]]$n_total))/10000 # 5.36 %
neq <- which(comp_w_k4[[1]]$n_total != comp_w_k4[[2]]$n_total)
mean(comp_w_k4[[1]]$n_total[neq] - comp_w_k4[[2]]$n_total[neq])
plot(comp_w_k4[[1]]$n_total, comp_w_k4[[2]]$n_total)
# look at Power:
mean(comp_w_k4[[1]]$actual_power - comp_w_k4[[2]]$actual_power)

# r = 0.4
comp_w_k5 <- comp_PO(iter=10000, r=0.4)
length(which(comp_w_k5[[1]]$n_total == comp_w_k5[[2]]$n_total))/10000 # 9.94%
# how often is AA the same as Kieser
length(which(comp_w_k5[[2]]$n_total == comp_w_k5[[3]]$n_total))/10000 # 1.75%
# how often is NA the same as Kieser
length(which(comp_w_k5[[2]]$n_total == comp_w_k5[[4]]$n_total))/10000 # 5.47 %
neq <- which(comp_w_k5[[1]]$n_total != comp_w_k5[[2]]$n_total)
mean(comp_w_k5[[1]]$n_total[neq] - comp_w_k5[[2]]$n_total[neq])
plot(comp_w_k5[[1]]$n_total, comp_w_k5[[2]]$n_total)
# look at Power:
mean(comp_w_k5[[1]]$actual_power - comp_w_k5[[2]]$actual_power)

############# Wilcoxon- Mann- Whitney Method ###################################
# This function calculates pi_A based on p_E and p_C
# Kieser Formula 4.4
calculate_pi_A<-function(p_C,p_E)
{
  x <- 0
  p_C_cum <- cumsum(p_C)
  for (i in 2:length(p_E)){
    x <- x + p_E[i] * p_C_cum[i-1]
  }
  y <- 0
  for (i in 1:length(p_E)){
    y <- y + 0.5 * p_E[i] * p_C[i]
  }
  return(x+y)  
}

# Wilcoxon-Mann-Whitney test for ordered categorical data (Kieser, 2020)
# This function (adapted from Kieser 2020) calculates the needed sample size based on p_C and p_E for a form of Mann-Whitney test taking ties into account. 
# p_C2 and p_E2: "true" Probabilites (unknown in real life, but "known" in Simulation)

samplesize_AfS <- function(p_C, p_E, alpha, beta, r, p_C2=NULL, p_E2=NULL){
  # calculate pi_A
  pi_A <- calculate_pi_A(p_C,p_E)
  
  # calculate sample size
  x <- 0
  for (i in 1:length(p_C)){
    x <- x + ((r*p_C[i] + p_E[i]) / (1 + r))^3
  } 
  n <- ((1+r)^2/r) * ((qnorm(1-(alpha/2)) + qnorm(1-beta))^2) / 12 *
    (1-x)/((pi_A - 0.5)^2) # Kieser Formula 4.7
  n_E <- ceiling(r/(1+r) * n)
  n_C <- ceiling(1/(1+r) * n)
  n_total <- n_E + n_C
  actual_r <- n_E/n_C
  actual_power <- pnorm(abs(pi_A-0.5) / (1 + actual_r)*
                          sqrt(n_total * actual_r * 12/(1-x)) - qnorm(1-alpha/2))
  actual_power2<-NA
  if (!is.null(p_C2)&!is.null(p_E2))
  {
  pi_A2<-calculate_pi_A(p_C2,p_E2)
  actual_power2 <- pnorm(abs(pi_A2-0.5) / (1 + actual_r)*
                          sqrt(n_total * actual_r * 12/(1-x)) - qnorm(1-alpha/2))
  }
  return(data.frame(n_total = n_total, n_E = n_E, n_C = n_C,
                    actual_r = actual_r, actual_power = actual_power, actual_power2=actual_power2))
}


############# t-test Method ####################################################

# This function calculates delta_A (difference between control and exp. group) and sigma based on p_C and p_E
# 
calculate_delta_A_sigma<-function(p_C,p_E){
  K<-length(p_C)
  mu_E<-sum((1:K)*p_E)
  mu_C<-sum((1:K)*p_C)
  delta_A<-mu_E-mu_C  
  sigma_C<- sqrt(sum(p_C*((1:K)-mu_C)^2))
  sigma_E<- sqrt(sum(p_E*((1:K)-mu_E)^2))
  sigma<-0.5*(sigma_C+sigma_E)
  list(delta_A=delta_A,sigma=sigma)
}

# This function (adapted from Kieser 2020) calculates the needed sample size based on p_C and p_E for the two-sample t-test using the Guenther & Schouten correction 

samplesize_ttestord<-  function(p_C, p_E, alpha, beta, r, p_C2=NULL, p_E2=NULL){
 
  param<-calculate_delta_A_sigma(p_C,p_E)
  # could be formula 8.4 from kieser 2020
  n_C <- ceiling((1+r)/r * (qnorm(1-alpha/2) + qnorm(1-beta))^2 *
                   (param$sigma/param$delta_A)^2 + (qnorm(1-alpha/2)^2) / (2*(1+r)))
  n_E <- ceiling((1+r) * (qnorm(1-alpha/2) + qnorm(1-beta))^2 *
                   (param$sigma/param$delta_A)^2 + r * (qnorm(1-alpha/2)^2) / (2*(1+r)))
  actual_r <- n_E/n_C
  n_total <- n_E + n_C
  actual_power <- pnorm(1/((1+actual_r) * abs(param$sigma/param$delta_A)) *
                          sqrt(actual_r*(n_total-(qnorm(1-alpha/2)^2)/
                                           (2
                                              #*(1+actual_r) # this should not be in the formula
                                              ))) - qnorm(1-alpha/2))
  actual_power2<-NA
  if (!is.null(p_C2)&!is.null(p_E2))
  {
    param2<-calculate_delta_A_sigma(p_C2,p_E2)
    actual_power2 <- pnorm(1/((1+actual_r) * abs(param2$sigma/param2$delta_A)) *
                             sqrt(actual_r*(n_total-(qnorm(1-alpha/2)^2)/
                                              (2
                                               #*(1+actual_r)
                                               ))) - qnorm(1-alpha/2))
  }
  return(data.frame(n_total = n_total, n_E = n_E, n_C = n_C,
                    actual_r = actual_r, actual_power = actual_power, actual_power2=actual_power2))
}
samplesize_ttestord_kieser<-  function(p_C, p_E, alpha, beta, r, p_C2=NULL, p_E2=NULL){
  
  param<-calculate_delta_A_sigma(p_C,p_E)
  # could be formula 8.4 from kieser 2020
  n_C <- ceiling((1+r)/r * (qnorm(1-alpha/2) + qnorm(1-beta))^2 *
                   (param$sigma/param$delta_A)^2 + (qnorm(1-alpha/2)^2) / (2*(1+r)))
  n_E <- ceiling((1+r) * (qnorm(1-alpha/2) + qnorm(1-beta))^2 *
                   (param$sigma/param$delta_A)^2 + r * (qnorm(1-alpha/2)^2) / (2*(1+r)))
  actual_r <- n_E/n_C
  n_total <- n_E + n_C
  actual_power <- pnorm(1/((1+actual_r) * abs(param$sigma/param$delta_A)) *
                          sqrt(actual_r*(n_total-(qnorm(1-alpha/2)^2)/
                                           (2*(1+actual_r)))) - qnorm(1-alpha/2))
  actual_power2<-NA
  if (!is.null(p_C2)&!is.null(p_E2))
  {
    param2<-calculate_delta_A_sigma(p_C2,p_E2)
    actual_power2 <- pnorm(1/((1+actual_r) * abs(param2$sigma/param2$delta_A)) *
                             sqrt(actual_r*(n_total-(qnorm(1-alpha/2)^2)/
                                              (2*(1+actual_r)))) - qnorm(1-alpha/2))
  }
  return(data.frame(n_total = n_total, n_E = n_E, n_C = n_C,
                    actual_r = actual_r, actual_power = actual_power, actual_power2=actual_power2))
}

samplesize_ttestord(p_C, p_E, r=1, alpha = 0.05, beta = 0.2)

# Compare the two Versions

comp_ttest <- function(alpha=0.05, beta=0.2, r=1, iter=1000, bias = 1.2){
  results_kieser <-  data.frame()
  results_me <-  data.frame()
  for (i in seq_along(1:iter)) {
    set.seed(i)
    prob_length <- sample(c(3,4,5,6,7,8), 1)
    print(i)
    p <- generate_two_simplex_vectors(prob_length, bias)
    p_C <- p[[1]]
    p_E <- p[[2]]
    results_kieser <- rbind(results_kieser, as.data.frame(samplesize_ttestord_kieser(p_C=p_C, p_E=p_E, alpha, beta, r)))
    results_me <- rbind(results_me, as.data.frame(samplesize_ttestord(p_C=p_C, p_E=p_E, alpha, beta, r)))
   }
  return(list(results_kieser, results_me))
}

c_ttest_1 <- comp_ttest(iter = 10000)
length(which(c_ttest_1[[1]]$n_total == c_ttest_1[[2]]$n_total))/10000 # for r=1 the sample size is always the same

# r=0.4
c_ttest_2 <- comp_ttest(iter = 10000, r = 0.4)
length(which(c_ttest_2[[1]]$n_total == c_ttest_2[[2]]$n_total))/10000 # for r=0.8 the sample size is always the same


############# Simulation Functions #############################################
# This function (written by ALB) generates RCT data of total size n_pilot with randomization ratio r and 
# probability vectors p_C and p_E in the control and experimental arms, respectively, and returns the 
# resulting estimates for p_C and p_E.
generate.pEpC.estimate<-function(p_C, p_E, r, n_pilot){
  y<-rbinom(n_pilot,1,r/(r+1))  
  n_E <- sum(y)
  n_C <- n_pilot - n_E
  phat_C <- as.vector(rmultinom(1, size=n_C, prob=p_C)/n_C)
  phat_E <- as.vector(rmultinom(1, size=n_E, prob=p_E)/n_E)
  list(C=phat_C,E=phat_E)
}


# This function (written by ALB) returns the distribution of the needed sample sizes calculated using methods ttest_ord, AfS and PO over niter pilot datasets.

simulation<-function(p_C, p_E, r=1, n_pilot, niter=1000, alpha=0.05,beta=0.2){
  n_needed<-matrix(NA,niter,3)
  colnames(n_needed)<-c("AfS","ttestord", "po")
  actual_power2<-matrix(NA,niter,3)
  colnames(actual_power2)<-c("AfS","ttestord", "po")
  for (i in 1:niter){
    set.seed(i)
    print(i)
    phat<-generate.pEpC.estimate(p_C=p_C, p_E=p_E, r=r, n_pilot=n_pilot)
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


############# Simulation Application ###########################################
# Small demo in an example setting.
p_C<-c(0.2,0.3,0.3,0.2)
p_E<-c(0.1, 0.2, 0.4, 0.3)
r<-1 

# compute needed sample sizes for these values of p_C and p_E:
samplesize_AfS(p_C, p_E, alpha=0.05, beta=0.2, r)
samplesize_ttestord(p_C, p_E, alpha=0.05, beta=0.2, r)
samplesize_po_kieser(p_C, p_E, alpha=0.05, beta=0.2, r)
samplesize_po_NA(p_C, p_E, alpha=0.05, beta=0.2, r)
samplesize_po_NN(p_C, p_E, alpha=0.05, beta=0.2, r)
rbind(samplesize_AfS(p_C, p_E, alpha=0.05, beta=0.2, r), 
      samplesize_ttestord(p_C, p_E, alpha=0.05, beta=0.2, r), 
      samplesize_po_kieser(p_C, p_E, alpha=0.05, beta=0.2, r),
      samplesize_po_NN(p_C, p_E, alpha=0.05, beta=0.2, r))

# perform test simulations for different values of n_pilot and compute how often n calculated by AfS is small than n calculated by ttestord
test50<-simulation(p_C,p_E,n_pilot=50,niter=10000)
mean(test50$n_needed[,1]<test50$n_needed[,2]) #how often n calculated by AfS is smaller than n calculated by ttestord
mean(test50$n_needed[,2]<test50$n_needed[,3])
mean(test50$n_needed[,1]<test50$n_needed[,3])
# which Method calculated the minimal sample size overall
colnames(test50$n_needed)[which(as.matrix(test50$n_needed) == min(test50$n_needed), arr.ind = TRUE)[1,2]]

mean(test50$actual_power_nmin)
median(test50$actual_power_nmin)
hist(test50$actual_power_nmin)
plot(test50$actual_power_nmin, type = 'p', col = factor(test50$method))
legend("topleft", legend = levels(factor(test50$method)), 
       pch = 19, col = factor(levels(factor(test50$method))))
df50 <- data.frame("power_nmin" = test50$actual_power_nmin, "method" = test50$method)
ggplot(df50, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")

test500<-simulation(p_C,p_E,n_pilot=500,niter=10000)
mean(test500$n_needed[,1]<test500$n_needed[,2])     # how often n calculated by AfS is smaller than n calculated by ttestord
mean(test500$n_needed[,2]<test500$n_needed[,3])     # how often n calculated by ttest is smaller than n calculated by PO
mean(test500$n_needed[,1]<test500$n_needed[,3])     # how often n calculated by AfS is smaller than n calculated by PO
c(mean(test500$n_needed[,1]<test500$n_needed[,2]), mean(test500$n_needed[,2]<test500$n_needed[,3]), 
  mean(test500$n_needed[,1]<test500$n_needed[,3]))
mean(test500$actual_power_nmin)
median(test500$actual_power_nmin)
hist(test500$actual_power_nmin)
plot(test500$actual_power_nmin, type = 'p', col = factor(test500$method))
legend("topleft", legend = levels(factor(test500$method)), 
       pch = 19, col = factor(levels(factor(test500$method))))
df500 <- data.frame("power_nmin" = test500$actual_power_nmin, "method" = test500$method)
ggplot(df500, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(test500$method == "PO"))/10000
length(which(test500$method == "ttest"))/10000

test1000<-simulation(p_C, p_E, n_pilot=1000, niter=10000)
mean(test1000$n_needed[,1]<test1000$n_needed[,2])
mean(test1000$actual_power_nmin)
hist(test1000$actual_power_nmin)
plot(test1000$actual_power_nmin, type = 'p', col = factor(test1000$method))
  legend("topleft", legend = levels(factor(test1000$method)), 
         pch = 19, col = factor(levels(factor(test1000$method))))
df1000 <- data.frame("power_nmin" = test1000$actual_power_nmin, "method" = test1000$method)
ggplot(df1000, aes(x = power_nmin, fill = method, colour = method)) + 
    geom_histogram(alpha = 0.3, position = "identity")
length(which(test1000$method == "PO"))/10000
length(which(test1000$method == "ttest"))/10000

#plot the actual_powers of the methods
hist(test1000$actual_power[,1])
hist(test1000$actual_power[,2])
hist(test1000$actual_power[,3])
### now layer them?

### maybe do this with the three po methods???

test10000_1<-simulation(p_C,p_E,n_pilot=10000,niter=1000)
mean(test10000_1$actual_power_nmin)
hist(test10000_1$actual_power_nmin)
plot(test10000_1$actual_power_nmin, type = 'p', col = factor(test10000_1$method))
legend("topleft", legend = levels(factor(test10000_1$method)), 
       pch = 19, col = factor(levels(factor(test10000_1$method))))
df10000_1 <- data.frame("power_nmin" = test10000_1$actual_power_nmin, "method" = test10000_1$method)
ggplot(df10000_1, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
  
test10000<-simulation(p_C,p_E,n_pilot=10000,niter=10000)
mean(test10000$actual_power_nmin)
hist(test10000$actual_power_nmin)
plot(test10000$actual_power_nmin, type = 'p', col = factor(test10000$method))
legend("topleft", legend = levels(factor(test10000$method)), 
       pch = 19, col = factor(levels(factor(test10000$method))))
df10000 <- data.frame("power_nmin" = test10000$actual_power_nmin, "method" = test10000$method)
ggplot(df10000, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(test10000$method == "PO"))/10000
length(which(test10000$method == "ttest"))/10000

test15000<-simulation(p_C,p_E,n_pilot=15000,niter=10000)
mean(test15000$actual_power_nmin)
hist(test15000$actual_power_nmin)
plot(test15000$actual_power_nmin, type = 'p', col = factor(test15000$method))
legend("topleft", legend = levels(factor(test15000$method)), 
       pch = 19, col = factor(levels(factor(test15000$method))))
df15000 <- data.frame("power_nmin" = test15000$actual_power_nmin, "method" = test15000$method)
ggplot(df15000, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")

### for r=0.8
test1000_0.8<-simulation(p_C,p_E,n_pilot=1000,niter=10000, r=0.8)
mean(test1000_0.8$actual_power_nmin)
hist(test1000_0.8$actual_power_nmin)
df1000_0.8 <- data.frame("power_nmin" = test1000_0.8$actual_power_nmin, "method" = test1000_0.8$method)
ggplot(df1000_0.8, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")

test10000_0.8<-simulation(p_C,p_E,n_pilot=10000,niter=10000, r=0.8)
mean(test10000_0.8$actual_power_nmin)
hist(test10000_0.8$actual_power_nmin)
df10000_0.8 <- data.frame("power_nmin" = test10000_0.8$actual_power_nmin, "method" = test10000_0.8$method)
ggplot(df10000_0.8, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")

### for r=0.4
test1000_0.4<-simulation(p_C,p_E,n_pilot=1000,niter=10000, r=0.4)
mean(test1000_0.4$actual_power_nmin)
hist(test1000_0.4$actual_power_nmin)

test10000_0.4<-simulation(p_C,p_E,n_pilot=10000,niter=10000, r=0.4)
mean(test10000_0.4$actual_power_nmin)
hist(test10000_0.4$actual_power_nmin)
df10000_0.4 <- data.frame("power_nmin" = test10000_0.4$actual_power_nmin, "method" = test10000_0.4$method)
ggplot(df10000_0.4, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")

### for r=0.2
test1000_0.2<-simulation(p_C,p_E,n_pilot=1000,niter=10000, r=0.2)
mean(test1000_0.2$actual_power_nmin)
hist(test1000_0.2$actual_power_nmin)
df1000_0.2 <- data.frame("power_nmin" = test1000_0.2$actual_power_nmin, "method" = test1000_0.2$method)
ggplot(df1000_0.2, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")

test10000_0.2<-simulation(p_C,p_E,n_pilot=10000,niter=10000, r=0.2)
mean(test10000_0.2$actual_power_nmin)
hist(test10000_0.2$actual_power_nmin)
df10000_0.2 <- data.frame("power_nmin" = test10000_0.2$actual_power_nmin, "method" = test10000_0.2$method)
ggplot(df10000_0.2, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")


###### Multiple Simulations in one plot #########
### Function to make Simulations in different settings
sim_min_power <- data.frame("n_pilot"=c(NA), "min_n" = c(NA), "min_n_power" = c(NA), "method" = c(NA))
sim_mean_power <- data.frame("n_pilot"=c(NA), "nmin_power"=c(NA))
settings <- function(n_vector, p_C, p_E, r, niter){
  for (i in seq_along(n_vector)) {
    sim <- simulation(p_C, p_E, r, n_pilot = n_vector[i], niter)
    sim_min_power[i,2]<- mean(min(sim$nmin))
    sim_min_power[i,3] <- mean(sim$actual_power[which.min(sim$nmin)])
    sim_min_power[i,4] <- sim$method[which.min(sim$nmin)]
    sim_min_power[i,1] <- n_vector[i]
    sim_mean_power[i,2] <- mean(sim$actual_power_nmin)
    sim_mean_power[i,1] <- n_vector[i]}
  return(list(sim_mean_power, sim_min_power))
}

sim_small_test <- settings(n_vector = c(seq(150,1000,50), seq(1500, 12000, 500)),
                           p_C, p_E, r, niter=10000)
ggplot(sim_small_test[[1]], aes(x=n_pilot, y=nmin_power)) +
  geom_point()+
  geom_hline(yintercept=0.8, color = "red")
ggplot(sim_small_test[[2]], aes(x=n_pilot, y=min_n, color = method)) +
  geom_point()
length(which(sim_small_test[[2]]$method =="PO"))/nrow(sim_small_test[[2]])


ggplot(sim_small_test[[2]], aes(x=n_pilot, y=min_n_power, color = method)) +
  geom_point()+
  geom_point(data = sim_small_test[[1]], mapping = aes(x = n_pilot, y = nmin_power), color="blue")+
  geom_hline(yintercept=0.8, color = "red")


sim_test <- settings(n_vector = c(seq(150,1000,50), seq(1500, 12000, 500)), p_C, p_E, r, niter = 1000)
ggplot(sim_test[[1]], aes(x=n_pilot, y=nmin_power)) +
  geom_point()+
  geom_hline(yintercept=0.8, color = "red")

# now with more iterations, to get rid of the variation
sim_test_iter <- settings(n_vector = c(seq(150,1000,50), seq(1500, 12000, 500)), p_C, p_E, r, niter = 10000)
ggplot(sim_test_iter[[1]], aes(x=n_pilot, y=nmin_power)) +
  geom_point()+
  geom_hline(yintercept=0.8, color = "red")
ggplot(sim_test_iter[[2]], aes(x=n_pilot, y=min_n, color = method)) +
  geom_point()

# r=0.8
sim_test0.8 <- settings(n_vector = c(seq(150,1000,50), seq(1500, 12000, 500)), p_C, p_E, r=0.8, niter = 10000)
ggplot(sim_test0.8[[1]], aes(x=n_pilot, y=nmin_power)) +
  geom_point()+
  geom_hline(yintercept=0.8, color = "red")
mean(sim_test0.8[[1]]$nmin_power)
ggplot(sim_test0.8[[2]], aes(x=n_pilot, y=min_n, color = method)) +
  geom_point()
length(which(sim_test0.8[[2]]$method=="PO"))/nrow(sim_test0.8[[2]])

# r=0.4
sim_test0.4 <- settings(n_vector = c(seq(150,1000,50), seq(1500, 12000, 500)), p_C, p_E, r=0.4, niter = 10000)
ggplot(sim_test0.4[[1]], aes(x=n_pilot, y=nmin_power)) +
  geom_point()+
  geom_hline(yintercept=0.8, color = "red")
mean(sim_test0.4[[1]]$nmin_power)
ggplot(sim_test0.4[[2]], aes(x=n_pilot, y=min_n, color = method)) +
  geom_point()
length(which(sim_test0.4[[2]]$method=="PO"))/nrow(sim_test0.4[[2]])

# r=0.2
sim_test0.2 <- settings(n_vector = c(seq(150,1000,50), seq(1500, 12000, 500)), p_C, p_E, r=0.2, niter = 10000)
ggplot(sim_test0.2[[1]], aes(x=n_pilot, y=nmin_power)) +
  geom_point()+
  geom_hline(yintercept=0.8, color = "red")
mean(sim_test0.2[[1]]$nmin_power)
ggplot(sim_test0.2[[2]], aes(x=n_pilot, y=min_n, color = method)) +
  geom_point()
length(which(sim_test0.2[[2]]$method=="PO"))/nrow(sim_test0.2[[2]])

###
p_C_2 <- c(0.2, 0.1, 0.15, 0.25, 0.05, 0.1, 0.15)
p_E_2 <- c(0.15, 0.15, 0.1, 0.25, 0.1, 0.1, 0.15)
sim_test2 <- settings(n_vector = c(seq(150,1000,50), seq(1500, 12000, 500)), p_C_2, p_E_2, r, niter = 1000)
ggplot(sim_test2, aes(x=nmin_power, y=n_pilot)) +
  geom_point()

p_C_3 <- c(0.2, 0.1, 0.15, 0.25, 0.05, 0.1, 0.15)
p_E_3 <- c(0.26, 0.08, 0.14, 0.24, 0.04, 0.09, 0.15)
sim_test3 <- settings(n_vector = c(seq(150,1000,50), seq(1500, 12000, 500)), p_C_3, p_E_3, r, niter = 1000)
ggplot(sim_test3, aes(x=nmin_power, y=n_pilot)) +
  geom_point()

p_C_4 <- random_simplex(6)
p_E_4 <- random_simplex(6)
sim_test4 <- settings(n_min = 100, n_max = 10000, steps = 100 , p_C_4, p_E_4, r, niter = 1000)
ggplot(sim_test4, aes(x=nmin_power, y=n_pilot)) +
  geom_point()


p_5 <- generate_two_simplex_vectors(6)
p_C_5 <- p_5[[1]]
p_E_5 <- p_5[[2]]

sim_test5 <- settings(n_vector = c(seq(150,1000,50), seq(1500, 12000, 500)), p_C_5, p_E_5, r, niter = 1000)
ggplot(sim_test5, aes(x=n_pilot, y=nmin_power)) +
  geom_point()



sim_test6 <- settings(n_vector = c(seq(150,1000,50), seq(1500, 12000, 500)), p_C, p_E, r, niter = 10000)
ggplot(sim_test6, aes(x=n_pilot, y=nmin_power)) +
  geom_point()
# -> weniger Schwankungen

##### Compare Methods with normally distributed Prob. #####
set.seed(1)
p_C_norm <- normal_simplex(6,4)
p_C_norm
set.seed(2)
p_E_norm <- normal_simplex(6,4)
p_E_norm

test50_norm <-simulation(p_C_norm,p_E_norm,n_pilot=50,niter=10000)
mean(test50_norm$n_needed[,1]<test50_norm$n_needed[,2]) #how often n calculated by AfS is smaller than n calculated by ttestord
mean(test50_norm$n_needed[,2]<test50_norm$n_needed[,3])
mean(test50_norm$n_needed[,1]<test50_norm$n_needed[,3])
# which Method calculated the minimal sample size overall
colnames(test50_norm$n_needed)[which(as.matrix(test50_norm$n_needed) == min(test50_norm$n_needed), arr.ind = TRUE)[1,2]]
mean(test50_norm$actual_power_nmin)
median(test50_norm$actual_power_nmin)
hist(test50_norm$actual_power_nmin)
plot(test50_norm$actual_power_nmin, type = 'p', col = factor(test50_norm$method))
legend("topleft", legend = levels(factor(test50_norm$method)), 
       pch = 19, col = factor(levels(factor(test50_norm$method))))
df50_norm <- data.frame("power_nmin" = test50_norm$actual_power_nmin, "method" = test50_norm$method)
ggplot(df50_norm, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
df50_power_norm <- as.data.frame(test50_norm$actual_power) %>%
  pivot_longer(cols = AfS:po, values_to = "power", names_to = "method")
ggplot(df50_power_norm, aes(x = power, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")

test500_norm<-simulation(p_C_norm,p_E_norm,n_pilot=500,niter=10000)
mean(test500_norm$n_needed[,1]<test500_norm$n_needed[,2])     # how often n calculated by AfS is smaller than n calculated by ttestord
mean(test500_norm$n_needed[,2]<test500_norm$n_needed[,3])     # how often n calculated by ttest is smaller than n calculated by PO
mean(test500_norm$n_needed[,1]<test500_norm$n_needed[,3])     # how often n calculated by AfS is smaller than n calculated by PO
c(mean(test500_norm$n_needed[,1]<test500_norm$n_needed[,2]),
  mean(test500_norm$n_needed[,2]<test500_norm$n_needed[,3]), 
  mean(test500_norm$n_needed[,1]<test500_norm$n_needed[,3]))
mean(test500_norm$actual_power_nmin)
median(test500_norm$actual_power_nmin)
hist(test500_norm$actual_power_nmin)
plot(test500_norm$actual_power_nmin, type = 'p', col = factor(test500_norm$method))
legend("topleft", legend = levels(factor(test500_norm$method)), 
       pch = 19, col = factor(levels(factor(test500_norm$method))))
df500_norm <- data.frame("power_nmin" = test500_norm$actual_power_nmin, "method" = test500_norm$method)
ggplot(df500_norm, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(test500_norm$method == "PO"))/10000
length(which(test500_norm$method == "ttest"))/10000
length(which(test500_norm$method == "AfS"))/10000

test1000_norm<-simulation(p_C_norm, p_E_norm, n_pilot=1000, niter=10000)
mean(test1000_norm$n_needed[,1]<test1000_norm$n_needed[,2])
mean(test1000_norm$actual_power_nmin)
hist(test1000_norm$actual_power_nmin)
# Plot of minimal sample size per Method
df1000_norm_nmin <- data.frame("nmin" = test1000_norm$nmin, "method" = test1000_norm$method)
ggplot(df1000_norm_nmin, aes(x = nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")+
  xlim(0,10000)

plot(test1000_norm$actual_power_nmin, type = 'p', col = factor(test1000_norm$method))
legend("topleft", legend = levels(factor(test1000_norm$method)), 
       pch = 19, col = factor(levels(factor(test1000_norm$method))))
# Plot of minimal power per Method
df1000_norm <- data.frame("power_nmin" = test1000_norm$actual_power_nmin, "method" = test1000_norm$method)
ggplot(df1000_norm, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(test1000_norm$method == "PO"))/10000
length(which(test1000_norm$method == "ttest"))/10000
length(which(test1000_norm$method == "AfS"))/10000

#plot the actual_powers of the methods
hist(test1000_norm$actual_power[,1])
hist(test1000_norm$actual_power[,2])
hist(test1000_norm$actual_power[,3])
### now layer them
df1000_power_norm <- as.data.frame(test1000_norm$actual_power) %>%
  pivot_longer(cols = AfS:po, values_to = "power", names_to = "method")
ggplot(df1000_power_norm, aes(x = power, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")


# n_pilot = 5000
test5000_norm<-simulation(p_C_norm, p_E_norm, n_pilot=5000, niter=10000)
mean(test5000_norm$n_needed[,1]<test5000_norm$n_needed[,2])
mean(test5000_norm$actual_power_nmin)
hist(test5000_norm$actual_power_nmin)
plot(test5000_norm$actual_power_nmin, type = 'p', col = factor(test5000_norm$method))
legend("topleft", legend = levels(factor(test5000_norm$method)), 
       pch = 19, col = factor(levels(factor(test5000_norm$method))))
df5000_norm <- data.frame("power_nmin" = test5000_norm$actual_power_nmin, "method" = test5000_norm$method)
ggplot(df5000_norm, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(test5000_norm$method == "PO"))/10000
length(which(test5000_norm$method == "ttest"))/10000
length(which(test5000_norm$method == "AfS"))/10000

# n_pilot=10000
test10000_norm<-simulation(p_C_norm, p_E_norm, n_pilot=10000, niter=10000)
mean(test10000_norm$n_needed[,1]<test10000_norm$n_needed[,2])
mean(test10000_norm$actual_power_nmin)
hist(test10000_norm$actual_power_nmin)
plot(test10000_norm$actual_power_nmin, type = 'p', col = factor(test10000_norm$method))
legend("topleft", legend = levels(factor(test10000_norm$method)), 
       pch = 19, col = factor(levels(factor(test10000_norm$method))))
df10000_norm <- data.frame("power_nmin" = test10000_norm$actual_power_nmin, "method" = test10000_norm$method)
ggplot(df10000_norm, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(test10000_norm$method == "PO"))/10000
length(which(test10000_norm$method == "ttest"))/10000
length(which(test10000_norm$method == "AfS"))/10000
