# This function will calculate the needed sample size based on p_C and p_E for proportional odds regression

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

calculate_theta_N <- function(p_C, p_E, r = 1) {
  # Validate inputs
  if (length(p_C) != length(p_E)) {
    stop("p_C and p_E must have the same length (same number of outcome levels).")
  }
  if (abs(sum(p_C) - 1) > 1e-6 || abs(sum(p_E) - 1) > 1e-6) {
    stop("Probabilities p_C and p_E must each sum to 1.")
  }
  
  K <- length(p_C)  # number of outcome levels
  
  p_i <- (r*p_C+p_E)/(r+1) # expected value of estimand

  # Construct outcome levels and group labels
  outcome <- rep(c(1:K), times = 2)
  group <- c(rep(c("control"), times = K), rep(c("treatment"), times = K))
  
  # Calculate weights under the null hypothesis
  weights_null <- c(((r * p_i) / (r + 1)), (p_i / (r + 1)))
  
  # Create synthetic dataset
  data_null <- data.frame(
    outcome = factor(outcome, ordered = TRUE),
    group = factor(group, levels = c("control", "treatment")),
    weight = weights_null
  )
  
  # Fit proportional odds model under the null hypothesis
  model_null <- suppressWarnings(
    polr(outcome ~ group, data = data_null, weights = weight, method = "logistic", Hess = TRUE)
  )
  
  # Extract coefficient of treatment effect
  coef_null <- -coef(summary(model_null))["grouptreatment", "Value"]
  var_null <- coef(summary(model_null))["grouptreatment", "Std. Error"]
  null <- data.frame("coef" = coef_null, "Var" = var_null)
  return(null)
}

samp_size_PO_White_2023 <- function(p_C, p_E, alpha, beta, r, method){
  theta_A <- calculate_theta_A(p_C, p_E, r)[1,1]
  Var_A <- calculate_theta_A(p_C, p_E, r)[1,2]^2
  Var_N <- calculate_theta_N(p_C, p_E, r)[1,2]^2
  
  if(method=="NA"){
    n <- ((sqrt(Var_N)*qnorm(1-alpha/2)+sqrt(Var_A)*qnorm(1-beta))^2)/(theta_A^2) # White, 2023 Formel 4 (NA)
    n_E <- ceiling(r/(1+r) * n)
    n_C <- ceiling(1/(1+r) * n)
    n_total <- n_E + n_C
    actual_r <- n_E/n_C
  }else if(method=="NN"){
    n <- (Var_N*((qnorm(1-alpha/2)+qnorm(1-beta))^2))/(theta_A^2) # White, 2023 Formel 5 (NN)
    n_E <- ceiling(r/(1+r) * n)
    n_C <- ceiling(1/(1+r) * n)
    n_total <- n_E + n_C
    actual_r <- n_E/n_C
  }else if(method=="AA"){
    n <- (Var_A*((qnorm(1-alpha/2)+qnorm(1-beta))^2))/(theta_A^2) # White, 2023 Formel 5 (AA)
    n_E <- ceiling(r/(1+r) * n)
    n_C <- ceiling(1/(1+r) * n)
    n_total <- n_E + n_C
    actual_r <- n_E/n_C
  }else print("choose from methods NA, NN, AA")
  return(c("n_total"=n_total, "n_E"=n_E, "n_C"=n_C, "actual_r"=actual_r))
}

##### Estimation check of theta_A #####
# check estimation precision of theta_A estimation
## calculate theta_A on probability vectors that fulfill the PO assumption repeatedly
theta_comparison <- function(comp_iter, vec_length, r){
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

t_comp <- theta_comparison(comp_iter = 10000, vec_length = 6, r=1)
mean(t_comp$diff) # 0.000503
median(t_comp$diff) # 0.00033
plot(t_comp$theta_A, t_comp$diff)  

##### Estimation check of variance N #####
# check estimation of variance under null hypothesis V_N compared to PO variance 
## calculate variance N on probability vectors that fulfill the PO assumption repeatedly and compare to variance of PO method
var_comparison <- function(comp_iter, vec_length, r){
  var_results <- data.frame("theta_A" = c(NA), "var_w"=c(NA), "var_k" = c(NA), "diff" = c(NA))
  for (i in seq_along(1:comp_iter)) {
    set.seed(i)
    print(i)
    theta_random <- runif(1, min = 1, max = 3)
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

v_comp <- var_comparison(comp_iter = 10000, vec_length = 6, r=1)
mean(v_comp$diff) # on average a difference of 0.00065
median(v_comp$diff) # 0.00031
plot(v_comp$theta_A, v_comp$diff) 


##### Compare NN Method to Whitehead method #####
# code to calculate PO sample size from (Kieser, 2020) formula 4.8
samplesize_po <- function(p_C, p_E, alpha, beta, r){
  theta_A <- calculate_theta_A(p_C, p_E, r)[1,1]
  x = 0
  for (i in 1:length(p_C)){
    x = x + ((p_C[i] + r*p_E[i]) / (1 + r))^3
  }
  n <- (1+r)^2/r * 3*(qnorm(1-alpha/2) + qnorm(1-beta))^2 /
    (theta_A^2 * (1-x)) # Kieser formula 4.8
  n_E <- ceiling(r/(1+r) * n)
  n_C <- ceiling(1/(1+r) * n)
  n_total <- n_E + n_C
  actual_r <- n_E/n_C
  return(c("n_total" = n_total, "n_E" = n_E, "n_C" = n_C, "actual_r" = actual_r))
}

n_comparison <- function(alpha=0.05, beta=0.2, r=1, iter=1000, theta_A = log(1.8)){
  results_NN <-  data.frame(n_total = numeric(0), n_E = numeric(0), n_C = numeric(0), actual_r = numeric(0))
  results_po <- data.frame(n_total = numeric(0), n_E = numeric(0), n_C = numeric(0), actual_r = numeric(0))
  theta_vec <- rep(NA, iter)
  for (i in seq_along(1:iter)) {
    set.seed(i)
    prob_length <- sample(c(3,4,5,6,7,8), 1)
    print(i)
    p_C <- random_simplex(prob_length)
    p_E <- calc_p_E(p_C, theta_A = theta_A)
    theta_vec[i] <- calculate_theta_A(p_C, p_E)[1,1]
    results_NN[nrow(results_NN)+1,] <- samp_size_PO_White_2023(p_C=p_C, p_E=p_E, alpha, beta, r, method="NN")
    results_po[nrow(results_po)+1,] <- samplesize_po(p_C=p_C, p_E=p_E, alpha, beta, r)
  }
  return(list(results_NN, results_po, theta_vec))
}
n_comp <- n_comparison(iter=10000)

length(which(n_comp[[1]]$n_total == n_comp[[2]]$n_total))/10000 # 99.08% equal
