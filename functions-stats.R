# Helper functions

#******************************************************************************************************************
# Statistical functions
#******************************************************************************************************************

# Perform (one-by-one, rudimentary; no interactions) variance decomposition on the bootstrap results
variance_decomp <- function(bootstrap_draws, output_var, var_output, calc_function) {
  vd <- tibble(Parameter = colnames(bootstrap_draws), variance = NA_real_)
  # For each parameter in turn, 'freeze' it (replace random draws with median,
  # removing all variability) and calculate variance in bootstrapped GPP
  for(i in seq_len(nrow(vd))) {
    x <- bootstrap_draws
    x[vd$Parameter[i]] <- median(unlist(x[vd$Parameter[i]]))
    vd$variance[i] <- var(unlist(calc_function(x)[output_var]))
  }
  vd$Contribution <- round((var_output - vd$variance) / var_output, 3)
  vd
}


# Given a vector `x` and the quantiles of the 'raw' distribution, i.e. the 
# standard top-down GPP or bottom-up Rs, compute the probability they agree
prob_agreement <- function(x, raw_quantiles) {
  stopifnot(length(raw_quantiles) == 3)
  
  # assuming x follows a normal distribution, we can calculate the probability
  # of it overlapping the distribution implied by the raw quantiles
  mu <- mean(x)
  s <- sd(x)
  p_high <- pnorm((raw_quantiles[3] - mu) / s)
  p_low <- pnorm((raw_quantiles[1] - mu) / s)
  (p_high - p_low) %>% round(4)
}


# Perform a Student's t-test and return results nicely formatted for a table
t_test <- function(x, y, alternative = "two.sided", ...) {
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(y))
  z <- t.test(x, y, alternative = alternative, ...)
  paste0("t = ", format(z$statistic, digits = 1, scientific = FALSE), 
         ", df = ", format(z$parameter, digits = 1, scientific = FALSE), 
         ", p-value = ", format(z$p.value, digits = 3, scientific = FALSE))
}

