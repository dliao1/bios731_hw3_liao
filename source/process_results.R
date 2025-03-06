# Processes results from each simulation scenario and summarizes in 
# a dataframe where each row corresponds to one combination of parameters
# Summary information includes biases, coverages for confidence intervals,
# power, type 1 error, and MCSEs for the previous informations.
# 
# Inputs: number of scenarios, number of simulations
# Return value: data frame with summary information
process_results <- function(num_scenarios, n_sim) {
  # Reads in finished scenario files
  all_wald_results <- vector("list", num_scenarios)
  all_boot_p_results <- vector("list", num_scenarios)
  
  # Loop through each scenario file
  for (i in 1:num_scenarios) {
    file_path <- here("results", "sim_wald", paste0("scenario_", i, ".RDA"))
    load(file_path)
    
    all_wald_results[[i]] <- all_wald_estim %>% 
      mutate(err_type = ifelse(err_type == 0, "Gamma", "Normal"))
    
    file_path <- here("results", "sim_boot_percentile", paste0("scenario_", i, ".RDA"))
    load(file_path)
    
    all_boot_p_results[[i]] <- all_boot_percent_estim %>%
      mutate(err_type = ifelse(err_type == 0, "Gamma", "Normal"))
  } 
  
  # Each list is of length 4 since 4 parameter combinations
  all_biases <- rep(NA, num_scenarios)
  var_hat <- rep(NA, num_scenarios)
  wald_coverage <- rep(NA, num_scenarios)
  boot_p_coverage <- rep(NA, num_scenarios)
  wald_time <- rep(NA, num_scenarios)
  boot_p_time <- rep(NA, num_scenarios)
  all_n <- rep(NA, num_scenarios)
  all_beta_true <- rep(NA, num_scenarios)
  all_err_type <- rep(NA, num_scenarios)
  scenario_num <- rep(NA, num_scenarios)
  mean_wald_se_beta <- rep(NA, num_scenarios)
  mean_boot_p_se_beta <- rep(NA, num_scenarios)
  se_hat <- rep(NA, num_scenarios)
  se_hat_se <- rep(NA, num_scenarios)
  wald_type1 <- rep(NA, num_scenarios)
  boot_p_type1 <- rep(NA, num_scenarios)
  wald_power <- rep(NA, num_scenarios)
  boot_p_power <- rep(NA, num_scenarios)
  
  bias_mcse <- rep(NA, num_scenarios)
  coverage_mcse <- rep(NA, num_scenarios)
  wald_power_mcse <- rep(NA, num_scenarios)
  wald_type1_mcse <- rep(NA, num_scenarios)
  
  
  for (i in 1:num_scenarios) {
    
    all_biases[i] <- (1/n_sim) * sum(all_wald_results[[i]]$beta_hat - all_wald_results[[i]]$beta_true)
    var_hat[i] <- sd(all_wald_results[[i]]$beta_hat)
    se_hat[i] <- mean(all_wald_results[[i]]$se_beta)
    se_hat_se[i] <- sd(all_wald_results[[i]]$se_beta)
    wald_coverage[i] <- mean(all_wald_results[[i]]$coverage == 1) 
    boot_p_coverage[i] <- mean(all_boot_p_results[[i]]$coverage == 1) 
    wald_time[i] <- mean(all_wald_results[[i]]$time)
    boot_p_time[i] <- mean(all_boot_p_results[[i]]$time)
    all_n[i] <- unique(all_wald_results[[i]]$n)[[1]]
    all_beta_true[i] <- unique(all_wald_results[[i]]$beta_true)
    all_err_type[i] <- unique(all_wald_results[[i]]$err_type)
    mean_wald_se_beta[i] <- mean(all_wald_results[[i]]$se_beta)
    mean_boot_p_se_beta[i] <- mean(all_boot_p_results[[i]]$se_beta)
    
    # Power = p(detecting an effect)
    # = p(rejecting null when null false)
    # = 1 - p(NOT rejecting null when null false)
    # = 1 - p(0 in CI when beta_true = 0.5)
    #  = p(0 not in CI when beta_true = 0.5)
    
    # type I error
    # p(rejecting null when null is true)
    # = p(0 not in CI when beta_true = 0)
    
    if (unique(all_wald_results[[i]]$beta_true)[[1]] == 0) {
      wald_type1[i] <- mean(all_wald_results[[i]]$ci_l > 0 | all_wald_results[[i]]$ci_u < 0)
      boot_p_type1[i] <- mean(all_boot_p_results[[i]]$ci_l > 0 | all_boot_p_results[[i]]$ci_u < 0)
      
      
    } else { # beta = 0.5
      wald_power[i] <- mean(all_wald_results[[i]]$ci_l > 0 | all_wald_results[[i]]$ci_u < 0)
      boot_p_power[i] <- mean(all_boot_p_results[[i]]$ci_l > 0 | all_boot_p_results[[i]]$ci_u < 0)
      
    }
    
    scenario_num[i] <- i
    bias_mcse[i] <- sqrt(sum((all_wald_results[[i]]$beta_hat - mean(all_wald_results[[i]]$beta_hat))^2) / (n_sim * (n_sim - 1)))
  }
  
  
  df <- bind_cols(
    scenario_num = scenario_num, 
    n = all_n, 
    beta_true = all_beta_true, 
    error_type = all_err_type, 
    bias = all_biases,
    bias_mcse = bias_mcse,
    var = var_hat,
    se_hat = se_hat,
    se_hat_se = se_hat_se,
    wald_coverage = wald_coverage,
    wald_coverage_mcse = sqrt(wald_coverage * (1-wald_coverage)/ n_sim) ,
    wald_time = wald_time,
    wald_power = wald_power,
    wald_power_mcse = sqrt(wald_power * (1-wald_power)/n_sim), 
    wald_type1 = wald_type1,
    wald_type1_mcse = sqrt(wald_type1 * (1-wald_type1)/n_sim),
    wald_se = mean_wald_se_beta,
    boot_p_coverage = boot_p_coverage,
    boot_p_coverage_mcse = sqrt(boot_p_coverage * (1-boot_p_coverage)/ n_sim) ,
    boot_p_time = boot_p_time,
    boot_p_se = mean_boot_p_se_beta,
    boot_p_power = boot_p_power,
    boot_p_power_mcse = sqrt(boot_p_power * (1-boot_p_power)/n_sim),
    boot_p_type1 = boot_p_type1,
    boot_p_type1_mcse = sqrt(boot_p_type1 * (1-boot_p_type1)/n_sim),
  )
  
  return(df)
  
}
