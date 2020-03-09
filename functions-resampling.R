# Helper functions

#******************************************************************************************************************
# Resampling functions
#******************************************************************************************************************


# * Method 1 - ignore SD and just resample with replacement
# * Method 2 - use SD; missing SDs calculated based on mean of coefficient of variability
# * Method 3 - like Method 2 but use the maximum CV if missing
# * Method 4 - using SD only when available, otherwise 0 SD

method1 <- function(n, x) { # ignore all errors
  mean(sample(x, n, replace = TRUE))
}
method2 <- function(n, x, s, f = median) { # replace NA errors with median error
  cv <- s / x
  cv <- replace_na(cv, f(cv, na.rm = TRUE))
  s <- if_else(is.na(s), x * cv, s)
  mean(rnorm(n, mean = x, sd = s))
}
method3 <- function(n, x, s) { # replace NA errors with max error
  method2(n, x, s, max)
}
method4 <- function(n, x, s) { # replace NA errors with zero
  method2(n, x, replace_na(s, 0))
}

# test sample size and repeat time
sample_size_tst <- function(test_size, N_SAMPLES) {

  # resample global soil respiration (RSG)
  resample_RSG <- tibble(
    Rs1_no_SD = sapply(rep(test_size, N_SAMPLES), method1, GlobalRs$Rs),
    Rs2_miss_SD_median = sapply(rep(test_size, N_SAMPLES), method2, GlobalRs$Rs, GlobalRs$SD),
    Rs3_miss_SD_max = sapply(rep(test_size, N_SAMPLES), method3, GlobalRs$Rs, GlobalRs$SD),
    Rs4_miss_SD_zero = sapply(rep(test_size, N_SAMPLES), method4, GlobalRs$Rs, GlobalRs$SD)
  )

  # comparing RSG resample results
  resample_RSG %>%
    gather(Rs_resample_type) %>%
    ggplot(aes(value, fill = Rs_resample_type)) +
    theme_cowplot() +
    geom_density(stat = "density", alpha = 0.5) +
    #  theme(legend.position = c(0.75, 0.75)) +
    # reverse the color scale so that teal (color 1) is always observations,
    # while pink is always implied fluxes
    scale_fill_discrete("", h.start = 180) +
    scale_color_discrete(h.start = 180) +
    ylab("Density") +
    xlab(expression(R[S] ~ (Pg ~ C ~ yr^{
      -1
    }))) -> RSG_sample_fig

  # resample gross primary production (GPP)
  resample_GPP <- tibble(
    GPP1_no_SD = sapply(rep(test_size, N_SAMPLES), method1, GPP$GPP),
    GPP2_miss_SD_median = sapply(rep(test_size, N_SAMPLES), method2, GPP$GPP, GPP$SD),
    GPP3_miss_SD_max = sapply(rep(test_size, N_SAMPLES), method3, GPP$GPP, GPP$SD),
    GPP4_miss_SD_zero = sapply(rep(test_size, N_SAMPLES), method4, GPP$GPP, GPP$SD)
  )

  # comparing GPP resample results
  resample_GPP %>%
    gather(GPP_resample_type) %>%
    ggplot(aes(value, fill = GPP_resample_type)) +
    theme_cowplot() +
    geom_density(stat = "density", alpha = 0.5) +
    #  theme(legend.position = c(0.75, 0.75)) +
    # reverse the color scale so that teal (color 1) is always observations,
    # while pink is always implied fluxes
    scale_fill_discrete("", h.start = 180) +
    scale_color_discrete(h.start = 180) +
    ylab("Density") +
    xlab(expression(GPP ~ (Pg ~ C ~ yr^{
      -1
    }))) -> GPP_sample_fig

  plot_grid(RSG_sample_fig, GPP_sample_fig, ncol = 1)
}
