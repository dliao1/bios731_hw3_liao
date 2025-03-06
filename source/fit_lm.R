# Fits linear regression model
# 
# Inputs: dataframe with x,y columns
# Returns: linear regression model
fit_model <- function(data) {
  model <- lm(y ~ x, data = data)
  return (model)
}