# Extracts estimates from linear model
#
# Input: linear regression model, true beta, alpha
# Returns: dataframe with columns: estimated beta, standard error,
# bias, upper/lower confidence interval limits, and coverage
extract_estims <- function(model, beta_true, alpha) {
  estims_df <- tidy(model, conf.int = TRUE) %>%
    filter(term == "x") %>%
    rename(beta_hat = estimate) %>%
    rename(se_beta = std.error) %>%
    mutate(beta_diff = beta_hat - beta_true) %>%
    mutate(ci_l = beta_hat - (qnorm(1 - alpha/2) *  se_beta)) %>%
    mutate(ci_u = beta_hat + (qnorm(1 - alpha/2) *  se_beta)) %>%
    mutate(coverage = ifelse(beta_true >= ci_l & beta_true <= ci_u, 1, 0)) %>%
    select(beta_hat, beta_diff, se_beta, ci_l, ci_u, coverage) %>%
  return (estims_df)
}