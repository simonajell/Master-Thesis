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


### different r #####

### write a function for this
vary_r <- function(p_C, p_E, n_pilot = 10000, niter = 10000){
  r_vec <- c(0.2, 0.4, 0.6, 0.8, 1)
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


# show all histograms
df_sim_r <- data.frame()
for (i in seq_along(1:length(sim_r))) {
  df <- data.frame("power_nmin" = sim_r[[i]]$actual_power_nmin, 
             "method" = sim_r[[i]]$method,
             "r"=c(0.2, 0.4, 0.6, 0.8, 1)[i])
  df_sim_r <- rbind(df_sim_r,df)
}
df_sim_r <- df_sim_r %>% group_by(r) %>%  mutate(mean = mean(power_nmin))

ggplot(df_sim_r, aes(x = power_nmin))+ 
  facet_grid(r ~.) +
  geom_histogram(fill = "grey", color = "black", position = "identity")+
  geom_vline(aes(xintercept = mean, group = r), color = "red")+
  geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
  theme_bw()
# separate by method
ggplot(df_sim_r, aes(x = power_nmin, fill = method, colour = method))+ 
  facet_grid(r ~.) +
  geom_histogram(alpha = 0.3, position = "identity")+
  geom_vline(aes(xintercept = mean, group = r), color = "red")+
  geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
  theme_bw()


### vary p_C and p_E #####
## vary bias
vary_p_bias <- function(n_pilot = 1000, niter = 10000, r = 1, cat = 6){
  bias_vec <- seq(0.1, 2.1, 0.2)
  p_results <- data.frame()
  results <- vector(mode = "list", length = length(bias_vec))
  for (i in seq_along(bias_vec)) {
    set.seed(2)
    p <- generate_two_simplex_vectors(cat, bias = bias_vec[i])
    p_C1 <- p[[1]]
    p_E1 <- p[[2]]
    sim <- simulation(p_C1, p_E1, n_pilot=n_pilot, niter=niter, r=r)
    results[[i]] <- sim
  }
  return(results)
}

sim_p_bias <- vary_p_bias(n_pilot = 10000)

mean_sim_p_bias <- c()
for (i in seq_along(1:length(sim_p_bias))) {
  mean_sim_p_bias[[i]] <- mean(sim_p_bias[[i]]$actual_power_nmin)
}
mean_sim_p_bias

# show all histograms
df_sim_p_bias <- data.frame()
for (i in seq_along(1:length(sim_p_bias))) {
  df <- data.frame("power_nmin" = sim_p_bias[[i]]$actual_power_nmin, 
                   "method" = sim_p_bias[[i]]$method,
                   "p" = c(seq(0.1, 2.1, 0.2))[i])
  df_sim_p_bias <- rbind(df_sim_p_bias,df)
}
df_sim_p_bias <- df_sim_p_bias %>% group_by(p) %>%  mutate(mean = mean(power_nmin))

ggplot(df_sim_p_bias, aes(x = power_nmin))+ 
  geom_histogram(fill = "grey", color = "black", position = "identity")+
  geom_vline(aes(xintercept = mean, group = p), color = "red")+
  geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
  facet_wrap(p ~.) +
  theme_bw()
# separate by method
ggplot(df_sim_p_bias, aes(x = power_nmin, fill = method, colour = method))+ 
  facet_wrap(p ~.) +
  geom_histogram(alpha = 0.3, position = "identity")+
  geom_vline(aes(xintercept = mean, group = p), color = "red")+
  geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
  theme_bw()

# just repeat different p with the same bias and different seeds
vary_p <- function(n_pilot = 1000, niter = 10000, r = 1, cat = 6, rep=10, bias = 2){
  vec <- c(1:rep)
  p_results <- data.frame()
  results <- vector(mode = "list", length = rep)
  for (i in seq_along(vec)) {
    set.seed(i)
    p <- generate_two_simplex_vectors_random(cat)
    p_C1 <- p[[1]]
    p_E1 <- p[[2]]
    sim <- simulation(p_C1, p_E1, n_pilot=n_pilot, niter=niter, r=r)
    results[[i]] <- sim
  }
  return(results)
}
sim_p <- vary_p(n_pilot = 10000)

mean_sim_p <- c()
for (i in seq_along(1:length(sim_p))) {
  mean_sim_p[[i]] <- mean(sim_p[[i]]$actual_power_nmin)
}
mean_sim_p

# show all histograms
df_sim_p <- data.frame()
for (i in seq_along(1:length(sim_p))) {
  df <- data.frame("power_nmin" = sim_p[[i]]$actual_power_nmin, 
                   "method" = sim_p[[i]]$method,
                   "p" = c(1:length(sim_p))[i])
  df_sim_p <- rbind(df_sim_p,df)
}
df_sim_p <- df_sim_p %>% group_by(p) %>%  mutate(mean = mean(power_nmin))

ggplot(df_sim_p, aes(x = power_nmin))+ 
  geom_histogram(fill = "grey", color = "black", position = "identity")+
  geom_vline(aes(xintercept = mean, group = p), color = "red")+
  geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
  facet_wrap(p ~.) +
  theme_bw()
# separate by method
ggplot(df_sim_p, aes(x = power_nmin, fill = method, colour = method))+ 
  facet_wrap(p ~.) +
  geom_histogram(alpha = 0.3, position = "identity")+
  geom_vline(aes(xintercept = mean, group = p), color = "red")+
  geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
  theme_bw()


##### normally distributed categories #####
# create vector that is normally distributed around a category (4)
set.seed(1)
p_C_norm <- normal_simplex(6,4)[[1]]
p_E_norm <- normal_simplex(6,4)[[2]]

sim_n <- simulation(p_C_norm, p_E_norm, n_pilot=10000, niter=10000)
mean(sim_n$actual_power_nmin)
hist(sim_n$actual_power_nmin)
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
vary_norm <- function(n_pilot = 1000, niter = 10000, r = 1, rep = 10, cat=2, vec_length = 6){
  rep_vec <- c(1:rep)
  results <- vector(mode = "list", length = rep)
  for (i in seq_along(1:rep)) {
    set.seed(i)
    p_C_norm <- normal_simplex(vec_length,cat)[[1]]
    p_E_norm <- normal_simplex(vec_length,cat)[[2]]
    sim <- simulation(p_C_norm, p_E_norm, r = r, niter = niter, n_pilot = n_pilot)
    results[[i]] <- sim
  }
  return(results)
}

sim_norm <- vary_norm(cat=4, n_pilot = 10000)

mean_sim_norm <- c()
for (i in seq_along(1:length(sim_norm))) {
  mean_sim_norm[[i]] <- mean(sim_norm[[i]]$actual_power_nmin)
}
mean_sim_norm

# show all histograms
df_sim_norm <- data.frame()
for (i in seq_along(1:length(sim_norm))) {
  df <- data.frame("power_nmin" = sim_norm[[i]]$actual_power_nmin, 
                   "method" = sim_norm[[i]]$method,
                   "seed_norm"=c(1:length(sim_norm))[i])
  df_sim_norm <- rbind(df_sim_norm,df)
}
df_sim_norm <- df_sim_norm %>% group_by(seed_norm) %>%  mutate(mean = mean(power_nmin))

ggplot(df_sim_norm, aes(x = power_nmin))+ 
  geom_histogram(fill = "grey", color = "black", position = "identity")+
  geom_vline(aes(xintercept = mean, group = r), color = "red")+
  geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
  facet_wrap(seed_norm ~.) +
  theme_bw()
# separate by method
ggplot(df_sim_norm, aes(x = power_nmin, fill = method, colour = method))+   
  geom_vline(aes(xintercept = mean, group = r), color = "red")+
  geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
  facet_wrap(seed_norm ~.,scales = "free") +
  geom_histogram(alpha = 0.3, position = "identity")+
  theme_bw()


# function with different extreme point locations
vary_norm_loc <- function(n_pilot = 1000, niter = 10000, r = 1, vec_length = 6, seed = 1){
  loc <- c(1:vec_length)
  results <- vector(mode = "list", length = vec_length)
  for (i in seq_along(loc)) {
    set.seed(seed)
    p_C_norm <- normal_simplex(vec_length,i)[[1]]
    p_E_norm <- normal_simplex(vec_length,i)[[2]]
    sim <- simulation(p_C_norm, p_E_norm, r = r, niter = niter, n_pilot = n_pilot)
    results[[i]] <- sim
  }
  return(results)
}
sim_norm_loc <- vary_norm_loc(n_pilot = 10000)

mean_sim_norm_loc <- c()
for (i in seq_along(1:length(sim_norm_loc))) {
  mean_sim_norm_loc[[i]] <- mean(sim_norm_loc[[i]]$actual_power_nmin)
}
mean_sim_norm_loc

# show all histograms
df_sim_norm_loc <- data.frame()
for (i in seq_along(1:length(sim_norm_loc))) {
  df <- data.frame("power_nmin" = sim_norm_loc[[i]]$actual_power_nmin, 
                   "method" = sim_norm_loc[[i]]$method,
                   "location"=c(1:length(sim_norm_loc))[i])
  df_sim_norm_loc <- rbind(df_sim_norm_loc,df)
}
df_sim_norm_loc <- df_sim_norm_loc %>% group_by(location) %>%  mutate(mean = mean(power_nmin))

ggplot(df_sim_norm_loc, aes(x = power_nmin))+ 
  facet_wrap(location ~.,) +
  geom_histogram(fill = "grey", color = "black", position = "identity")+
  geom_vline(aes(xintercept = mean, group = location), color = "red")+
  geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
  theme_bw()
# separate by method
ggplot(df_sim_norm_loc, aes(x = power_nmin, fill = method, colour = method))+ 
  facet_wrap(location ~.) +
  geom_histogram(alpha = 0.3, position = "identity")+
  geom_vline(aes(xintercept = mean, group = location), color = "red")+
  geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
  theme_bw()

# repeat for different seed
sim_norm_loc2 <- vary_norm_loc(n_pilot = 10000, seed = 10)

mean_sim_norm_loc2 <- c()
for (i in seq_along(1:length(sim_norm_loc2))) {
  mean_sim_norm_loc2[[i]] <- mean(sim_norm_loc2[[i]]$actual_power_nmin)
}
mean_sim_norm_loc2

# show all histograms
df_sim_norm_loc2 <- data.frame()
for (i in seq_along(1:length(sim_norm_loc2))) {
  df <- data.frame("power_nmin" = sim_norm_loc2[[i]]$actual_power_nmin, 
                   "method" = sim_norm_loc2[[i]]$method,
                   "location"=c(1:length(sim_norm_loc2))[i])
  df_sim_norm_loc2 <- rbind(df_sim_norm_loc2,df)
}
df_sim_norm_loc2 <- df_sim_norm_loc2 %>% group_by(location) %>%  mutate(mean = mean(power_nmin))

ggplot(df_sim_norm_loc2, aes(x = power_nmin))+ 
  facet_wrap(location ~.,) +
  geom_histogram(fill = "grey", color = "black", position = "identity")+
  geom_vline(aes(xintercept = mean, group = location), color = "red")+
  geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
  theme_bw()
# separate by method
ggplot(df_sim_norm_loc2, aes(x = power_nmin, fill = method, colour = method))+ 
  facet_wrap(location ~.) +
  geom_histogram(alpha = 0.3, position = "identity")+
  geom_vline(aes(xintercept = mean, group = location), color = "red")+
  geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
  theme_bw()

# repeat for different seed
sim_norm_loc3 <- vary_norm_loc(n_pilot = 10000, seed = 3, vec_length = 8)

mean_sim_norm_loc3 <- c()
for (i in seq_along(1:length(sim_norm_loc3))) {
  mean_sim_norm_loc3[[i]] <- mean(sim_norm_loc3[[i]]$actual_power_nmin)
}
mean_sim_norm_loc3

# show all histograms
df_sim_norm_loc3 <- data.frame()
for (i in seq_along(1:length(sim_norm_loc3))) {
  df <- data.frame("power_nmin" = sim_norm_loc3[[i]]$actual_power_nmin, 
                   "method" = sim_norm_loc3[[i]]$method,
                   "location"=c(1:length(sim_norm_loc3))[i])
  df_sim_norm_loc3 <- rbind(df_sim_norm_loc3,df)
}
df_sim_norm_loc3 <- df_sim_norm_loc3 %>% group_by(location) %>%  mutate(mean = mean(power_nmin))

ggplot(df_sim_norm_loc3, aes(x = power_nmin))+ 
  facet_wrap(location ~.,) +
  geom_histogram(fill = "grey", color = "black", position = "identity")+
  geom_vline(aes(xintercept = mean, group = location), color = "red")+
  geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
  theme_bw()
# separate by method
ggplot(df_sim_norm_loc3, aes(x = power_nmin, fill = method, colour = method))+ 
  facet_wrap(location ~.) +
  geom_histogram(alpha = 0.3, position = "identity")+
  geom_vline(aes(xintercept = mean, group = location), color = "red")+
  geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
  theme_bw()

### different number of categories #####

# function
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
set.seed(1)
uniform_simplex(7)


sim_cat <- vary_prob_length(niter=10000, n_pilot=10000, r = 1)
mean_sim_cat <- c()
for (i in seq_along(1:13)) {
  mean_sim_cat[[i]] <- mean(sim_cat[[i]]$actual_power_nmin)
}
mean_sim_cat

# show all histograms
df_sim_cat <- data.frame()
for (i in seq_along(1:length(sim_cat))) {
  df <- data.frame("power_nmin" = sim_cat[[i]]$actual_power_nmin, 
                   "method" = sim_cat[[i]]$method,
                   "cat" = c(3:15)[i])
  df_sim_cat <- rbind(df_sim_cat,df)
}
df_sim_cat <- df_sim_cat %>% group_by(cat) %>%  mutate(mean = mean(power_nmin))

ggplot(df_sim_cat, aes(x = power_nmin))+ 
  geom_histogram(fill = "grey", color = "black", position = "identity")+
  geom_vline(aes(xintercept = mean, group = cat), color = "red")+
  geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
  facet_wrap(cat ~.) +
  theme_bw()
# separate by method
ggplot(df_sim_cat, aes(x = power_nmin, fill = method, colour = method))+ 
  facet_wrap(cat ~.) +
  geom_vline(aes(xintercept = mean, group = cat), color = "red")+
  geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
  geom_histogram(alpha = 0.3, position = "identity")+
  theme_bw()

###
sim_cat2 <- vary_prob_length(niter=10000, n_pilot=10000, r = 1, seed=2)
mean_sim_cat2 <- c()
for (i in seq_along(1:13)) {
  mean_sim_cat2[[i]] <- mean(sim_cat2[[i]]$actual_power_nmin)
}
mean_sim_cat2

# show all histograms
df_sim_cat2 <- data.frame()
for (i in seq_along(1:length(sim_cat2))) {
  df <- data.frame("power_nmin" = sim_cat2[[i]]$actual_power_nmin, 
                   "method" = sim_cat2[[i]]$method,
                   "cat" = c(3:15)[i])
  df_sim_cat2 <- rbind(df_sim_cat2,df)
}
df_sim_cat2 <- df_sim_cat2 %>% group_by(cat) %>%  mutate(mean = mean(power_nmin))

ggplot(df_sim_cat2, aes(x = power_nmin))+ 
  geom_histogram(fill = "grey", color = "black", position = "identity")+
  geom_vline(aes(xintercept = mean, group = cat), color = "red")+
  geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
  facet_wrap(cat ~.) +
  theme_bw()
# separate by method
ggplot(df_sim_cat2, aes(x = power_nmin, fill = method, colour = method))+ 
  facet_wrap(cat ~.) +
  geom_vline(aes(xintercept = mean, group = cat), color = "red")+
  geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
  geom_histogram(alpha = 0.3, position = "identity")+
  theme_bw()



### different n_pilot #####

### write function that does all this
vary_npilot <- function(p_C, p_E, niter = 10000, r = 1){
  npilot_vec <- c(50, 100, 200, 500, 1000, 2500, 5000, 10000, 15000, 20000)
  results <- vector(mode = "list", length = length(npilot_vec))
  for (i in seq_along(npilot_vec)) {
    sim <- simulation(p_C, p_E, r = r, niter = niter, n_pilot = npilot_vec[i])
    results[[i]] <- sim
  }
  return(results)
}

sim_npilot <- vary_npilot(p_C, p_E)

mean_sim_npilot <- c()
for (i in seq_along(1:length(sim_npilot))) {
  mean_sim_npilot[[i]] <- mean(sim_npilot[[i]]$actual_power_nmin)
}
mean_sim_npilot



# show all histograms
df_sim_npilot <- data.frame()
for (i in seq_along(1:length(sim_npilot))) {
  df <- data.frame("power_nmin" = sim_npilot[[i]]$actual_power_nmin, 
                   "method" = sim_npilot[[i]]$method,
                   "n_pilot"=c(50, 100, 200, 500, 1000, 2500, 5000, 10000, 15000, 20000)[i])
  df_sim_npilot <- rbind(df_sim_npilot,df)
}
df_sim_npilot <- df_sim_npilot %>% group_by(n_pilot) %>%  mutate(mean = mean(power_nmin))

ggplot(df_sim_npilot, aes(x = power_nmin))+ 
  facet_wrap(n_pilot ~.,) +
  geom_histogram(fill = "grey", color = "black", position = "identity")+
  geom_vline(aes(xintercept = mean, group = n_pilot), color = "red")+
  geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
  theme_bw()
# separate by method
ggplot(df_sim_npilot, aes(x = power_nmin, fill = method, colour = method))+ 
  facet_wrap(n_pilot ~., scales = "free") +
  geom_histogram(alpha = 0.4, position = "identity")+
  geom_vline(aes(xintercept = mean, group = n_pilot), color = "red")+
  geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
  theme_bw()

### different iterations #####

### write function that does all this
vary_niter <- function(p_C, p_E, n_pilot = 1000, r = 1){
  niter_vec <- c(50, 250, 500, 1000, 2500, 5000, 10000, 15000, 20000)
  results <- vector(mode = "list", length = length(niter_vec))
  for (i in seq_along(niter_vec)) {
    sim <- simulation(p_C, p_E, r = r, n_pilot = n_pilot, niter = niter_vec[i])
    results[[i]] <- sim
  }
  return(results)
}

sim_niter <- vary_niter(p_C, p_E, n_pilot = 10000)

mean_sim_niter <- c()
for (i in seq_along(1:length(sim_niter))) {
  mean_sim_niter[[i]] <- mean(sim_niter[[i]]$actual_power_nmin)
}
mean_sim_niter


# show all histograms
df_sim_niter <- data.frame()
for (i in seq_along(1:length(sim_niter))) {
  df <- data.frame("power_nmin" = sim_niter[[i]]$actual_power_nmin, 
                   "method" = sim_niter[[i]]$method,
                   "niter" = c(50, 250, 500, 1000, 2500, 5000, 10000, 15000, 20000)[i])
  df_sim_niter <- rbind(df_sim_niter,df)
}
df_sim_niter <- df_sim_niter %>% group_by(niter) %>%  mutate(mean = mean(power_nmin))

ggplot(df_sim_niter, aes(x = power_nmin))+ 
  facet_wrap(niter ~.) +
  geom_histogram(fill = "grey", color = "black", position = "identity")+
  geom_vline(aes(xintercept = mean, group = niter), color = "red")+
  geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
  
  theme_bw()
# separate by method
ggplot(df_sim_niter, aes(x = power_nmin, fill = method, colour = method))+ 
  facet_wrap(niter ~.,scales = "free") +
  geom_histogram(alpha = 0.3, position = "identity")+
  geom_vline(aes(xintercept = mean, group = niter), color = "red")+
  geom_vline(xintercept = 0.8, color = "red", linetype="dotted")+
  theme_bw()
