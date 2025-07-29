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
### for r=0.8
sim_0.8 <- simulation(p_C,p_E,n_pilot=1000,niter=10000, r=0.8)
mean(sim_0.8$actual_power_nmin) # Mean actual power of minimal sample
hist(sim_0.8$actual_power_nmin) # Histogram of power of minimal sample
df_sim_0.8 <- data.frame("power_nmin" = sim_0.8$actual_power_nmin, "method" = sim_0.8$method)
ggplot(df_sim_0.8, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(sim_0.8$method == "PO"))/10000
length(which(sim_0.8$method == "ttest"))/10000
length(which(sim_0.8$method == "AfS"))/10000
sim_0.8$actual_power[which(as.matrix(sim_0.8$n_needed) == min(sim_0.8$n_needed), arr.ind = TRUE)] # min actual power: 0.403

### for r=0.6
sim_0.6<-simulation(p_C,p_E,n_pilot=1000,niter=10000, r=0.6)
mean(sim_0.6$actual_power_nmin)
hist(sim_0.6$actual_power_nmin)
df_sim_0.6 <- data.frame("power_nmin" = sim_0.6$actual_power_nmin, "method" = sim_0.6$method)
ggplot(df_sim_0.6, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(sim_0.6$method == "PO"))/10000
length(which(sim_0.6$method == "ttest"))/10000
length(which(sim_0.6$method == "AfS"))/10000
sim_0.6$actual_power[which(as.matrix(sim_0.6$n_needed) == min(sim_0.6$n_needed), arr.ind = TRUE)] # min actual power: 0.3755

### for r=0.4
sim_0.4<-simulation(p_C,p_E,n_pilot=1000,niter=10000, r=0.4)
mean(sim_0.4$actual_power_nmin)
hist(sim_0.4$actual_power_nmin)
df_sim_0.4 <- data.frame("power_nmin" = sim_0.4$actual_power_nmin, "method" = sim_0.4$method)
ggplot(df_sim_0.4, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(sim_0.4$method == "PO"))/10000
length(which(sim_0.4$method == "ttest"))/10000
length(which(sim_0.4$method == "AfS"))/10000
sim_0.4$actual_power[which(as.matrix(sim_0.4$n_needed) == min(sim_0.4$n_needed), arr.ind = TRUE)] # min actual power: 0.358


### for r=0.2
sim_0.2<-simulation(p_C,p_E,n_pilot=1000,niter=10000, r=0.2)
mean(sim_0.2$actual_power_nmin)
hist(sim_0.2$actual_power_nmin)
df_sim_0.2 <- data.frame("power_nmin" = sim_0.2$actual_power_nmin, "method" = sim_0.2$method)
ggplot(df_sim_0.2, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(sim_0.2$method == "PO"))/10000
length(which(sim_0.2$method == "ttest"))/10000
length(which(sim_0.2$method == "AfS"))/10000
sim_0.2$actual_power[which(as.matrix(sim_0.2$n_needed) == min(sim_0.2$n_needed), arr.ind = TRUE)] # min actual power: 0.617


### vary p_C and p_E #####
##### normally distributed categories #####
# create vector that is normally distributed around a category (4)
set.seed(1)
p_C_norm <- normal_simplex(6,4)[[1]]
p_E_norm <- normal_simplex(6,4)[[2]]


sim_n <- simulation(p_C_norm, p_E_norm, n_pilot=1000, niter=10000)
mean(sim_n$actual_power_nmin)
hist(sim_n$actual_power_nmin)
# Plot of minimal power per Method
df_sim_n <- data.frame("power_nmin" = sim_n$actual_power_nmin, "method" = sim_n$method)
ggplot(df_sim_n, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")+
  theme_minimal()
length(which(sim_n$method == "PO"))/10000
length(which(sim_n$method == "ttest"))/10000
length(which(sim_n$method == "AfS"))/10000
sim_n$actual_power[which(as.matrix(sim_n$n_needed) == min(sim_n$n_needed), arr.ind = TRUE)] # min actual power: 0.403

# with more n_pilot
sim_n_2 <- simulation(p_C_norm, p_E_norm, n_pilot=5000, niter=10000)
mean(sim_n_2$actual_power_nmin)
hist(sim_n_2$actual_power_nmin)
# Plot of minimal power per Method
df_sim_n_2 <- data.frame("power_nmin" = sim_n_2$actual_power_nmin, "method" = sim_n_2$method)
ggplot(df_sim_n_2, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")+
  theme_minimal()
length(which(sim_n_2$method == "PO"))/10000
length(which(sim_n_2$method == "ttest"))/10000
length(which(sim_n_2$method == "AfS"))/10000
sim_n_2$actual_power[which(as.matrix(sim_n_2$n_needed) == min(sim_n_2$n_needed), arr.ind = TRUE)] # min actual power: 0.403

# now with normal distribution around 1
set.seed(1)
p_C_norm1 <- normal_simplex(6,1)[[1]]
p_E_norm1 <- normal_simplex(6,1)[[2]]

sim_n1 <- simulation(p_C_norm1, p_E_norm1, n_pilot=5000, niter=10000)
mean(sim_n1$actual_power_nmin)
hist(sim_n1$actual_power_nmin)
# Plot of minimal power per Method
df_sim_n1 <- data.frame("power_nmin" = sim_n1$actual_power_nmin, "method" = sim_n1$method)
ggplot(df_sim_n1, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")+
  theme_minimal()
length(which(sim_n1$method == "PO"))/10000
length(which(sim_n1$method == "ttest"))/10000
length(which(sim_n1$method == "AfS"))/10000
sim_n1$actual_power[which(as.matrix(sim_n1$n_needed) == min(sim_n1$n_needed), arr.ind = TRUE)] # min actual power: 0.403


### different number of categories #####
### 15 categories
set.seed(1)
p_C_15 <- generate_two_simplex_vectors(15)$a
p_E_15 <- generate_two_simplex_vectors(15)$b
sim_15 <- simulation(p_C_15,p_E_15,n_pilot=1000,niter=10000, r=1)
mean(sim_15$actual_power_nmin) # Mean actual power of minimal sample
hist(sim_15$actual_power_nmin) # Histogram of power of minimal sample
df_sim_15 <- data.frame("power_nmin" = sim_15$actual_power_nmin, "method" = sim_15$method)
ggplot(df_sim_15, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(sim_15$method == "PO"))/10000
length(which(sim_15$method == "ttest"))/10000
length(which(sim_15$method == "AfS"))/10000
sim_15$actual_power[which(as.matrix(sim_15$n_needed) == min(sim_15$n_needed), arr.ind = TRUE)] # min actual power

### 10 categories
set.seed(1)
p_C_10 <- generate_two_simplex_vectors(10)$a
p_E_10 <- generate_two_simplex_vectors(10)$b
sim_10 <- simulation(p_C_10,p_E_10,n_pilot=1000,niter=10000, r=1)
mean(sim_10$actual_power_nmin) # Mean actual power of minimal sample
hist(sim_10$actual_power_nmin) # Histogram of power of minimal sample
df_sim_10 <- data.frame("power_nmin" = sim_10$actual_power_nmin, "method" = sim_10$method)
ggplot(df_sim_10, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(sim_10$method == "PO"))/10000
length(which(sim_10$method == "ttest"))/10000
length(which(sim_10$method == "AfS"))/10000
sim_10$actual_power[which(as.matrix(sim_10$n_needed) == min(sim_10$n_needed), arr.ind = TRUE)] # min actual power

### 8 categories
set.seed(1)
p_C_8 <- generate_two_simplex_vectors(8)$a
p_E_8 <- generate_two_simplex_vectors(8)$b
sim_8 <- simulation(p_C_8,p_E_8,n_pilot=1000,niter=10000, r=1)
mean(sim_8$actual_power_nmin) # Mean actual power of minimal sample
hist(sim_8$actual_power_nmin) # Histogram of power of minimal sample
df_sim_8 <- data.frame("power_nmin" = sim_8$actual_power_nmin, "method" = sim_8$method)
ggplot(df_sim_8, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(sim_8$method == "PO"))/10000
length(which(sim_8$method == "ttest"))/10000
length(which(sim_8$method == "AfS"))/10000
sim_8$actual_power[which(as.matrix(sim_8$n_needed) == min(sim_8$n_needed), arr.ind = TRUE)] # min actual power

### 6 categories
set.seed(1)
p_C_6 <- generate_two_simplex_vectors(6)$a
p_E_6 <- generate_two_simplex_vectors(6)$b
sim_6 <- simulation(p_C_6,p_E_6,n_pilot=1000,niter=10000, r=1)
mean(sim_6$actual_power_nmin) # Mean actual power of minimal sample
hist(sim_6$actual_power_nmin) # Histogram of power of minimal sample
df_sim_6 <- data.frame("power_nmin" = sim_6$actual_power_nmin, "method" = sim_6$method)
ggplot(df_sim_6, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(sim_6$method == "PO"))/10000
length(which(sim_6$method == "ttest"))/10000
length(which(sim_6$method == "AfS"))/10000
sim_6$actual_power[which(as.matrix(sim_6$n_needed) == min(sim_6$n_needed), arr.ind = TRUE)] # min actual power

### 4 categories
set.seed(1)
p_C_4 <- generate_two_simplex_vectors(4)$a
p_E_4 <- generate_two_simplex_vectors(4)$b
sim_4 <- simulation(p_C_4,p_E_4,n_pilot=1000,niter=10000, r=1)
mean(sim_4$actual_power_nmin) # Mean actual power of minimal sample
hist(sim_4$actual_power_nmin) # Histogram of power of minimal sample
df_sim_4 <- data.frame("power_nmin" = sim_4$actual_power_nmin, "method" = sim_4$method)
ggplot(df_sim_4, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(sim_4$method == "PO"))/10000
length(which(sim_4$method == "ttest"))/10000
length(which(sim_4$method == "AfS"))/10000
sim_4$actual_power[which(as.matrix(sim_4$n_needed) == min(sim_4$n_needed), arr.ind = TRUE)] # min actual power

### 3 categories
set.seed(1)
p_C_3 <- generate_two_simplex_vectors(3)$a
p_E_3 <- generate_two_simplex_vectors(3)$b
sim_3 <- simulation(p_C_3,p_E_3,n_pilot=1000,niter=10000, r=1)
mean(sim_3$actual_power_nmin) # Mean actual power of minimal sample
hist(sim_3$actual_power_nmin) # Histogram of power of minimal sample
df_sim_3 <- data.frame("power_nmin" = sim_3$actual_power_nmin, "method" = sim_3$method)
ggplot(df_sim_3, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(sim_3$method == "PO"))/10000
length(which(sim_3$method == "ttest"))/10000
length(which(sim_3$method == "AfS"))/10000
sim_3$actual_power[which(as.matrix(sim_3$n_needed) == min(sim_3$n_needed), arr.ind = TRUE)] # min actual power

##### try to make vectors comparable
set.seed(1)
p_C_u <- uniform_simplex(6)[[1]]
p_E_u <- uniform_simplex(6)[[2]]

sim_u_6 <- simulation(p_C_u, p_E_u, n_pilot=1000, niter=10000, r=1)
mean(sim_u_6$actual_power_nmin) # Mean actual power of minimal sample
hist(sim_u_6$actual_power_nmin) # Histogram of power of minimal sample
df_sim_u_6 <- data.frame("power_nmin" = sim_u_6$actual_power_nmin, "method" = sim_u_6$method)
ggplot(df_sim_u_6, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(sim_u_6$method == "PO"))/10000
length(which(sim_u_6$method == "ttest"))/10000
length(which(sim_u_6$method == "AfS"))/10000
sim_u_6$actual_power[which(as.matrix(sim_u_6$n_needed) == min(sim_u_6$n_needed), arr.ind = TRUE)] # min actual power

## now different categories
# 3
set.seed(1)
p_C_u3 <- uniform_simplex(3)[[1]]
p_E_u3 <- uniform_simplex(3)[[2]]

sim_u_3 <- simulation(p_C_u3, p_E_u3, n_pilot=1000, niter=10000, r=1)
mean(sim_u_4$actual_power_nmin) # Mean actual power of minimal sample
hist(sim_u_4$actual_power_nmin) # Histogram of power of minimal sample
df_sim_u_4 <- data.frame("power_nmin" = sim_u_4$actual_power_nmin, "method" = sim_u_4$method)
ggplot(df_sim_u_6, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(sim_u_4$method == "PO"))/10000
length(which(sim_u_4$method == "ttest"))/10000
length(which(sim_u_4$method == "AfS"))/10000
sim_u_4$actual_power[which(as.matrix(sim_u_4$n_needed) == min(sim_u_4$n_needed), arr.ind = TRUE)] # min actual power

# 4
set.seed(1)
p_C_u4 <- uniform_simplex(4)[[1]]
p_E_u4 <- uniform_simplex(4)[[2]]

sim_u_4 <- simulation(p_C_u4, p_E_u4, n_pilot=1000, niter=10000, r=1)
mean(sim_u_4$actual_power_nmin) # Mean actual power of minimal sample
hist(sim_u_4$actual_power_nmin) # Histogram of power of minimal sample
df_sim_u_4 <- data.frame("power_nmin" = sim_u_4$actual_power_nmin, "method" = sim_u_4$method)
ggplot(df_sim_u_6, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(sim_u_4$method == "PO"))/10000
length(which(sim_u_4$method == "ttest"))/10000
length(which(sim_u_4$method == "AfS"))/10000
sim_u_4$actual_power[which(as.matrix(sim_u_4$n_needed) == min(sim_u_4$n_needed), arr.ind = TRUE)] # min actual power

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

vec_length_test <- vary_prob_length()
mean_sim_u <- c()
for (i in seq_along(1:13)) {
  mean_sim_u[[i]] <- mean(vec_length_test[[i]]$actual_power_nmin)
}
mean_sim_u
for (i in seq_along(1:13)) {
  hist(vec_length_test[[i]]$actual_power_nmin)
}

vec_length_test2 <- vary_prob_length(seed=2)
mean_sim_u <- c()
for (i in seq_along(1:13)) {
  mean_sim_u[[i]] <- mean(vec_length_test[[i]]$actual_power_nmin)
}
mean_sim_u
for (i in seq_along(1:13)) {
  hist(vec_length_test[[i]]$actual_power_nmin)
}


### different n_pilot #####
# n_pilot = 50
sim_50<-simulation(p_C,p_E,n_pilot=50,niter=10000)
mean(sim_50$actual_power_nmin)
median(sim_50$actual_power_nmin)
hist(sim_50$actual_power_nmin)
# plot of min power for the minimal sample size
df_sim_50 <- data.frame("power_nmin" = sim_50$actual_power_nmin, "method" = sim_50$method)
ggplot(df_sim_50, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(sim_50$method == "PO"))/10000
length(which(sim_50$method == "ttest"))/10000
length(which(sim_50$method == "AfS"))/10000
sim_50$actual_power[which(as.matrix(sim_50$n_needed) == min(sim_50$n_needed), arr.ind = TRUE)] # min actual power: 0.403
# plot of minimal sample size
df_sim_50_n <- data.frame("nmin" = sim_50$nmin, "method" = sim_50$method)
ggplot(df_sim_50_n, aes(x = nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")+
  xlim(0,4000)

# n_pilot = 250
sim_250<-simulation(p_C,p_E,n_pilot=250,niter=10000)
mean(sim_250$actual_power_nmin)
median(sim_250$actual_power_nmin)
hist(sim_250$actual_power_nmin)
# plot of min power for the minimal sample size
df_sim_250 <- data.frame("power_nmin" = sim_250$actual_power_nmin, "method" = sim_250$method)
ggplot(df_sim_250, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(sim_250$method == "PO"))/10000
length(which(sim_250$method == "ttest"))/10000
length(which(sim_250$method == "AfS"))/10000
sim_250$actual_power[which(as.matrix(sim_250$n_needed) == min(sim_250$n_needed), arr.ind = TRUE)] # min actual power: 0.403

# plot of minimal sample size
df_sim_250_n <- data.frame("nmin" = sim_250$nmin, "method" = sim_250$method)
ggplot(df_sim_250_n, aes(x = nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")+
  xlim(0,2500)

# n_pilot = 500
sim_500<-simulation(p_C,p_E,n_pilot=500,niter=10000)
mean(sim_500$actual_power_nmin)
median(sim_500$actual_power_nmin)
hist(sim_500$actual_power_nmin)
# plot of min power for the minimal sample size
df_sim_500 <- data.frame("power_nmin" = sim_500$actual_power_nmin, "method" = sim_500$method)
ggplot(df_sim_500, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(sim_500$method == "PO"))/10000
length(which(sim_500$method == "ttest"))/10000
length(which(sim_500$method == "AfS"))/10000
sim_500$actual_power[which(as.matrix(sim_500$n_needed) == min(sim_500$n_needed), arr.ind = TRUE)] # min actual power: 0.403

# plot of minimal sample size
df_sim_500_n <- data.frame("nmin" = sim_500$nmin, "method" = sim_500$method)
ggplot(df_sim_500_n, aes(x = nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")+
  xlim(0,1000)

# n_pilot = 1000
sim_1000<-simulation(p_C,p_E,n_pilot=1000,niter=10000)
mean(sim_1000$actual_power_nmin)
median(sim_1000$actual_power_nmin)
hist(sim_1000$actual_power_nmin)
# plot of min power for the minimal sample size
df_sim_1000 <- data.frame("power_nmin" = sim_1000$actual_power_nmin, "method" = sim_1000$method)
ggplot(df_sim_1000, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(sim_1000$method == "PO"))/10000
length(which(sim_1000$method == "ttest"))/10000
length(which(sim_1000$method == "AfS"))/10000
mean(sim_1000$actual_power[which(as.matrix(sim_1000$n_needed) == min(sim_1000$n_needed), arr.ind = TRUE)]) # min actual power: 0.403

# plot of minimal sample size
df_sim_1000_n <- data.frame("nmin" = sim_1000$nmin, "method" = sim_1000$method)
ggplot(df_sim_00_n, aes(x = nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")+
  xlim(0,1000)

# n_pilot = 5000
sim_5000<-simulation(p_C,p_E,n_pilot=5000,niter=10000)
mean(sim_5000$actual_power_nmin)
median(sim_5000$actual_power_nmin)
hist(sim_5000$actual_power_nmin)
# plot of min power for the minimal sample size
df_sim_5000 <- data.frame("power_nmin" = sim_5000$actual_power_nmin, "method" = sim_5000$method)
ggplot(df_sim_5000, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(sim_5000$method == "PO"))/10000
length(which(sim_5000$method == "ttest"))/10000
length(which(sim_5000$method == "AfS"))/10000
mean(sim_5000$actual_power[which(as.matrix(sim_5000$n_needed) == min(sim_5000$n_needed), arr.ind = TRUE)]) # min actual power: 0.403

# n_pilot = 10000
sim_10000<-simulation(p_C,p_E,n_pilot=10000,niter=10000)
mean(sim_10000$actual_power_nmin)
median(sim_10000$actual_power_nmin)
hist(sim_10000$actual_power_nmin)
# plot of min power for the minimal sample size
df_sim_10000 <- data.frame("power_nmin" = sim_10000$actual_power_nmin, "method" = sim_10000$method)
ggplot(df_sim_10000, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(sim_10000$method == "PO"))/10000
length(which(sim_10000$method == "ttest"))/10000
length(which(sim_10000$method == "AfS"))/10000
mean(sim_10000$actual_power[which(as.matrix(sim_10000$n_needed) == min(sim_10000$n_needed), arr.ind = TRUE)]) # min actual power: 0.403

### different iterations #####
# iter = 100
sim_i100<-simulation(p_C,p_E,n_pilot=1000,niter=100)
mean(sim_i100$actual_power_nmin)
median(sim_i100$actual_power_nmin)
hist(sim_i100$actual_power_nmin)
# plot of min power for the minimal sample size
df_sim_i100 <- data.frame("power_nmin" = sim_i100$actual_power_nmin, "method" = sim_i100$method)
ggplot(df_sim_i100, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(sim_i100$method == "PO"))/100
length(which(sim_i100$method == "ttest"))/100
length(which(sim_i100$method == "AfS"))/100
mean(sim_i100$actual_power[which(as.matrix(sim_i100$n_needed) == min(sim_i100$n_needed), arr.ind = TRUE)]) # min actual power: 0.403

# iter = 500
sim_i500<-simulation(p_C,p_E,n_pilot=1000,niter=500)
mean(sim_i500$actual_power_nmin)
median(sim_i500$actual_power_nmin)
hist(sim_i500$actual_power_nmin)
# plot of min power for the minimal sample size
df_sim_i500 <- data.frame("power_nmin" = sim_i500$actual_power_nmin, "method" = sim_i500$method)
ggplot(df_sim_i500, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(sim_i500$method == "PO"))/500
length(which(sim_i500$method == "ttest"))/500
length(which(sim_i500$method == "AfS"))/500
mean(sim_i500$actual_power[which(as.matrix(sim_i500$n_needed) == min(sim_i500$n_needed), arr.ind = TRUE)]) # min actual power: 0.403

# iter = 1000
sim_i1000<-simulation(p_C,p_E,n_pilot=1000,niter=1000)
mean(sim_i1000$actual_power_nmin)
median(sim_i1000$actual_power_nmin)
hist(sim_i1000$actual_power_nmin)
# plot of min power for the minimal sample size
df_sim_i1000 <- data.frame("power_nmin" = sim_i1000$actual_power_nmin, "method" = sim_i1000$method)
ggplot(df_sim_i1000, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(sim_i1000$method == "PO"))/1000
length(which(sim_i1000$method == "ttest"))/1000
length(which(sim_i1000$method == "AfS"))/1000
mean(sim_i1000$actual_power[which(as.matrix(sim_i1000$n_needed) == min(sim_i1000$n_needed), arr.ind = TRUE)]) # min actual power: 0.403

# iter = 5000
sim_i5000<-simulation(p_C,p_E,n_pilot=1000,niter=5000)
mean(sim_i5000$actual_power_nmin)
median(sim_i5000$actual_power_nmin)
hist(sim_i5000$actual_power_nmin)
# plot of min power for the minimal sample size
df_sim_i5000 <- data.frame("power_nmin" = sim_i5000$actual_power_nmin, "method" = sim_i5000$method)
ggplot(df_sim_i5000, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(sim_i5000$method == "PO"))/5000
length(which(sim_i5000$method == "ttest"))/5000
length(which(sim_i5000$method == "AfS"))/5000
mean(sim_i5000$actual_power[which(as.matrix(sim_i5000$n_needed) == min(sim_i5000$n_needed), arr.ind = TRUE)]) # min actual power: 0.403


# iter = 10000
sim_i10000<-simulation(p_C,p_E,n_pilot=1000,niter=10000)
mean(sim_i10000$actual_power_nmin)
median(sim_i10000$actual_power_nmin)
hist(sim_i10000$actual_power_nmin)
# plot of min power for the minimal sample size
df_sim_i10000 <- data.frame("power_nmin" = sim_i10000$actual_power_nmin, "method" = sim_i10000$method)
ggplot(df_sim_i10000, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(sim_i10000$method == "PO"))/10000
length(which(sim_i10000$method == "ttest"))/10000
length(which(sim_i10000$method == "AfS"))/10000
mean(sim_i10000$actual_power[which(as.matrix(sim_i10000$n_needed) == min(sim_i10000$n_needed), arr.ind = TRUE)]) # min actual power: 0.403

# iter = 15000
sim_i15000<-simulation(p_C,p_E,n_pilot=1000,niter=15000)
mean(sim_i15000$actual_power_nmin)
median(sim_i15000$actual_power_nmin)
hist(sim_i15000$actual_power_nmin)
# plot of min power for the minimal sample size
df_sim_i15000 <- data.frame("power_nmin" = sim_i15000$actual_power_nmin, "method" = sim_i15000$method)
ggplot(df_sim_i15000, aes(x = power_nmin, fill = method, colour = method)) + 
  geom_histogram(alpha = 0.3, position = "identity")
length(which(sim_i15000$method == "PO"))/15000
length(which(sim_i15000$method == "ttest"))/15000
length(which(sim_i15000$method == "AfS"))/15000
mean(sim_i15000$actual_power[which(as.matrix(sim_i15000$n_needed) == min(sim_i15000$n_needed), arr.ind = TRUE)]) # min actual power: 0.403
