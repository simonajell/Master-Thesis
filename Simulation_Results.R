
# First run the file "Functions.R" to load the packages and function needed in this File

#### Visualise the sample sizes for different allocations ####
  # Function to calculate the group sample sizes of a total sample size for different allocation ratios
  all_rat <- function(n=100, r_vec=c(seq(0.1, 0.9, 0.2), 1, seq(2,7,1))){
    # Create empty vectors
    n_C_vec <- c()
    n_E_vec <- c()
    
    # Loop over the different allocation ratios
    for (r in seq_along(1:length(r_vec))) {
      # Calculate the group sample sizes and round them up
      n_E_vec[r] <- ceiling(r_vec[r]/(1+r_vec[r]) * n)
      n_C_vec[r] <- ceiling(1/(1+r_vec[r]) * n)
    }
    return(data.frame("n" = c(n_E_vec, n_C_vec), 
                      "r" = rep(r_vec, 2),
                      "group"=c(rep("n_E", length(n_E_vec)),rep("n_C", length(n_C_vec))) ) )
  }
  
  # Apply function
  r_vec <- all_rat()
  
  # Plot the sample sizes and corresponding allocation ratios
  ggplot(r_vec, aes(x=r, y=n, color=group))+
    geom_line( )+
    geom_point( )+
    scale_color_manual(values = c("n_C" = "#669933", "n_E" = "#336699")) +
    ylab("Sample Size")+
    xlab("Allocation Ratio")+
    theme_bw()
  