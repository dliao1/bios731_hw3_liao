# Generates bootstrap data
#
# Inputs: Original data dataframe, true beta, sample size, alpha,
# number of bootstrap samples to take
# Returned value: nboot long list of estimated betas 
get_boot_data <- function (original_data, beta_true, sample_size, alpha, nboot) {
  all_boot_betas <- rep(NA, nboot)
  
  # Draws nboot number of bootstrap samples, fits model from current bootstrap
  # sample, and adds estimated beta to list.
  for (i in 1:nboot) {
    boot_sample <- slice_sample(original_data, n = sample_size, replace = TRUE)
    boot_sample_model <- fit_model(data = boot_sample)
    
    # Extract beta_treatment
    model_estims <- extract_estims(boot_sample_model, beta_true, alpha)
    all_boot_betas[i] <- ifelse(nrow(model_estims) == 1, model_estims$beta_hat, NA)
  }
  return (all_boot_betas)
}