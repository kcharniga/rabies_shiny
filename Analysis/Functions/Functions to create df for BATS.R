calc_prob_new <- function(size, prob, n_sims, 
                          prob_s_given_r, prob_s_given_not_r,
                          prob_h_given_r, prob_h_given_not_r,
                          prob_unprov_given_r, prob_unprov_given_not_r,
                          prob_prov_given_r, prob_prov_given_not_r) {
  
  # p1: positivity rate for state and animal
  if (prob[1] != 0){
    p1 <- rbinom(n = n_sims, 
                 size = size[1],
                 prob = prob[1])/size[1] # Just the first value
  } else {
    # If prob = 0 (i.e. pocket pets), sampling from the posterior distribution of binomial distribution
    # combined with conjugate prior distribution beta (uniform). Beta is the resulting distribution
    p1 <- rbeta(n = n_sims, 
                shape1 = 1 + 0, 
                shape2 = 1 + size[1] - 0) 
  }
  
  # Matrix to store results
    bat_mat <- matrix(NA, nrow = n_sims, ncol = 4)
  
  for (i in 1:length(p1)){
    
    # First, update the positivity rate for health status using
    # Bayes Rule
    p_sick    <- (prob_s_given_r*p1[i])/(prob_s_given_r*p1[i]+prob_s_given_not_r*(1-p1[i]))
    p_healthy <- (prob_h_given_r*p1[i])/(prob_h_given_r*p1[i]+prob_h_given_not_r*(1-p1[i]))
    
    # Next, adjust for whether the exposure was provoked
    
    # Sick
    p_sick_unprovoked <- (prob_unprov_given_r*p_sick)/(prob_unprov_given_r*p_sick+prob_unprov_given_not_r*(1-p_sick))
    p_sick_provoked   <- (prob_prov_given_r  *p_sick)/(prob_prov_given_r  *p_sick+prob_prov_given_not_r  *(1-p_sick))
    
    # Healthy
    p_healthy_unprovoked <- (prob_unprov_given_r*p_healthy)/(prob_unprov_given_r*p_healthy+prob_unprov_given_not_r*(1-p_healthy))
    p_healthy_provoked   <- (prob_prov_given_r  *p_healthy)/(prob_prov_given_r  *p_healthy+prob_prov_given_not_r  *(1-p_healthy))
    
      # Store probabilities in a matrix
      # NA for vaccination
      bat_mat[i,1] <- p_sick_unprovoked
      bat_mat[i,2] <- p_sick_provoked
      bat_mat[i,3] <- p_healthy_unprovoked
      bat_mat[i,4] <- p_healthy_provoked
    
  }
  
  # Calculate median and 95% CIs
  res <- apply(bat_mat, 2, FUN = quantile, probs = c(0.025, 0.5, 0.975))
  
  # Transpose the df
  final_df <- t(res)
  
  # Add column names and other variables
  colnames(final_df) <- c("lower_95CI",
                          "point_est",
                          "upper_95CI")
  
  return(final_df)
}

