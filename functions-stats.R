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

# Calculate overlap between two bootstrap samples
# This implements the method laid out in `toyproblem_overlap.Rmd`
# Returns the intersection point
calc_overlap <- function(left_sample, right_sample, intersection_interval) {
  # Test for normality
  left_normal_test <- shapiro.test(left_sample)
  left_normal <- left_normal_test$p.value >= 0.05
  print(left_normal_test)
  right_normal_test <- shapiro.test(right_sample)
  right_normal <- right_normal_test$p.value >= 0.05
  print(right_normal_test)
  
  # Calculate mu, sigma, and both normal and non-normal density functions
  mu_left <- mean(left_sample)
  sigma_left <- sd(left_sample)
  mu_right <- mean(right_sample)
  sigma_right <- sd(right_sample)
  
  if(mu_left > mu_right) {
    stop("It doesn't look like left_sample is actually on the left")  
  }
  
  density_left <- density(left_sample, n = 1024)
  normal_fn_left <- function(x) dnorm(x, mean = mu_left, sd = sigma_left)
  numerical_fn_left <- approxfun(x = density_left$x, y = density_left$y)
  density_right <- density(right_sample, n = 1024)
  normal_fn_right <- function(x) dnorm(x, mean = mu_right, sd = sigma_right)
  numerical_fn_right <- approxfun(x = density_right$x, y = density_right$y)
  
  # Decide which density functions we're using
  if(left_normal) {
    left_dist <- normal_fn_left
  } else {
    left_dist <- numerical_fn_left
  }
  if(right_normal) {
    right_dist <- normal_fn_right
  } else {
    right_dist <- numerical_fn_right
  }
  
  # Calculate the intersection point
  intersection <- uniroot(function(x) (left_dist(x) - right_dist(x)), 
                          interval = intersection_interval)$root
  cat("Intersection =", intersection, "Pg C/yr\n")
  
  # Plot the distributions and intersection
  tibble(x = seq(min(left_sample) - 50, max(right_sample) + 50, by = 5),
         left_normal = normal_fn_left(x),
         right_normal = normal_fn_right(x),
         left_non_normal = numerical_fn_left(x),
         right_non_normal = numerical_fn_right(x)) %>%
    gather(group, density, -x) %>% 
    ggplot(aes(x, density, color = group)) + geom_line() + 
    geom_vline(xintercept = intersection, linetype = 2) +
    annotate("label", x = mu_left, y = 0.01, label = left_normal) +
    annotate("label", x = mu_right, y = 0.01, label = right_normal) ->
    p
  print(p)
  
  # Calculate percentage of each distribution in overlap region
  
  # For vectors of x-axis points `x` and corresponding function evaluation points `y`,
  # the area under the curve according to the trapezoidal rule is:
  traprule <- function(x, y) {
    sum(diff(x) * (head(y, -1) + tail(y, -1)), na.rm = TRUE) / 2
  }
  
  # Left sample: define the x-axis points for integration
  dx <- 0.01   # integration step size
  left_xaxis_range <- seq(intersection, 4 * max(left_sample), by = dx ) # 4 is arbitrary
  # get the corresponding y-axis values
  left_yaxis_range <- left_dist(left_xaxis_range)
  # and get the area under the curve
  left_auc <- traprule(x = left_xaxis_range, y = left_yaxis_range)
  # Right sample calculated similarly except we can't go below zero, an easy cutoff
  right_xaxis_range <- seq(0, intersection, by = dx )
  right_yaxis_range <- right_dist(right_xaxis_range)
  right_auc <- traprule(x = right_xaxis_range, y = right_yaxis_range)
  
  cat("Left sample overlap =", round(left_auc * 100, 1), "%\n")
  cat("Right sample overlap =", round(right_auc * 100, 1), "%\n")
  
  cat("Left quantile of threshold:", ecdf(left_sample)(intersection), "\n")
  cat("Right quantile of threshold:", ecdf(right_sample)(intersection), "\n")
  intersection
}

