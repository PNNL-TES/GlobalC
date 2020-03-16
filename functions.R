# Helper functions

#******************************************************************************************************************
# data processing function - data calculation, clean, and combine etc
#******************************************************************************************************************

# Function to prepare Rroot/Ra ratio and Rshoot/Ra ratio data
# Ra = Rroot + Rshoot (total autotrophic respiration)
# Rroot/Ra or Rshoot/Ra ratio were from two resources:
# 1) collected from published article, stored in FsFlFr.csv; 2) from srdb_v4.csv
# and results from above two data sources will be combined
cal_Froot <- function(Froot_data, srdb_data) {
  # Froot and Fshoot data from digitized papers
  Froot_data %>%
    select(Latitude, Longitude, Fshoot, Froot, Ecosystem, IGBP) %>%
    mutate(
      Fshoot = if_else(is.na(Fshoot), 100 - Froot, as.numeric(Fshoot)),
      Froot = if_else(is.na(Froot), 100 - Fshoot, as.numeric(Froot)),
      # data in FsFlFr were in percentage, scale to decimal (to match srdb_v4))
      Froot = Froot / 100,
      Fshoot = Fshoot / 100
    ) ->
  comb_data

  # use data from srdb_v4 to calculate Rroot/Ra ratio, Rroot/(RE-NPP)
  # Ra_annual in srdb_v4 is Rroot
  # Ra = Rroot + Rshoot = ER - Rh
  srdb_data %>%
    # select(Study_number, Biome, Ecosystem_type, Leaf_habit, Latitude, Longitude, Ra_annual, ER, GPP, NPP) %>%
    select(Latitude, Longitude, Ecosystem_type, Leaf_habit, Rs_annual, Ra_annual, Rh_annual, ER, GPP, NPP) %>%
    filter(!is.na(GPP) | !is.na(NPP) | !is.na(ER)) %>%
    mutate(
      ER = if_else(is.na(ER), GPP, ER), # If ER not available, use GPP (assume NEP very small)
      # If Rh_annual not available, use NPP (assume NEP very small)
      Rh_annual = if_else(is.na(Rh_annual), NPP, Rh_annual),
      Ra_annual = if_else(is.na(Ra_annual), Rs_annual - Rh_annual, Ra_annual)
    ) ->
  sub_srdb

  # If RE available but Ra_annual is not, we could also estimate Ra_annual based on the Bond_Lamberty (2004) model
  # We can get 22 more data, but don't know whether should keep those points (it does not change the results too much)
  # sub_srdb$Ra_annual <- ifelse(is.na(sub_srdb$Ra_annual), (-7.97+0.93*sqrt(sub_srdb$Rs_annual))^2, sub_srdb$Ra_annual)

  # In srdb Ra_annual = Rroot, so Froot = Rroot / Ra (fraction of Rroot to autotrophic respiration)
  sub_srdb$Froot <- sub_srdb$Ra_annual / (sub_srdb$ER - sub_srdb$Rh_annual)
  sub_srdb <- subset(sub_srdb, !is.na(Froot) & Froot < 0.9 & Froot > 0.1)

  # mean(sub_srdb$Froot)

  sub_srdb$Fshoot <- 1 - sub_srdb$Froot
  sub_srdb %>% select(Latitude, Longitude, Fshoot, Froot, Ecosystem_type, Leaf_habit) -> sub_srdb
  colnames(sub_srdb) <- colnames(comb_data)

  comb_data <- rbind(comb_data, sub_srdb)
  comb_data[comb_data$Ecosystem == "Cropland", ]$Ecosystem <- "Agriculture"
  return(comb_data)
}


# Prepare and calculate Ra-GPP ratio data
# data are from two sources:
# 1) from published papers, stored in RaGPP.csv; 2) from srdb_v4.csv
# data from above two sources will be combined
prep_RaGpp <- function(RaGPP_data, srdb_v4) {
  RaGPP_data %>%
    select(Ecosystem, Leaf_habit, RaGPP_ratio) %>%
    mutate(Source = "Meta-papers") -> comb_data

  srdb_v4 %>%
    select(Ecosystem_type, Leaf_habit, Rs_annual, Ra_annual, Rh_annual, ER, GPP, NPP) %>%
    filter(!is.na(GPP) | !is.na(NPP) | !is.na(ER)) %>%
    mutate(
      ER = if_else(is.na(ER), GPP, ER), # assume NEP very small
      RaGPP_ratio = (ER - NPP) / GPP,
      Source = "srdb_v4"
    ) %>%
    filter(RaGPP_ratio > 0, RaGPP_ratio < 1, !is.na(Ecosystem_type)) %>%
    select(Ecosystem_type, Leaf_habit, RaGPP_ratio, Source) %>%
    bind_rows(comb_data)
}



# Given a length 3 vector of quantiles (2.5%, 50%, 97.5%) return a nice label
make_label <- function(qt_data) {
  stopifnot(length(qt_data) == 3)
  paste0("Mean and 95% CI: ", qt_data[2], " (", qt_data[1], ", ", qt_data[3], ")")
}

# The primary top-down and bottom-up bootstrapping figures all have similar
# annotations. The following two functions handle this.
# Given a figure, the types (character vector) of information it has, the two
# quantile vectors (see above), and position information, return annotated figure
annotations <- function(figure, types, x_qt, x_qt_raw, xpos) {
  stopifnot(length(x_qt) == 3)
  stopifnot(length(x_qt_raw) == 3)

  ann_dat <- tibble(
    Type = types,
    GPP_type = NA, Rs_type = NA,
    x = xpos, y = c(0.035, 0.025),
    label = c(make_label(x_qt), make_label(x_qt_raw))
  )
  figure + geom_text(
    data = ann_dat,
    aes(x = x, y = y, label = label, color = Type),
    hjust = 0, show.legend = FALSE
  )
}



#******************************************************************************************************************
# functions to make figures
#******************************************************************************************************************
# plot RC (root respiration to soil respiration ratio), Ra-GPP ratio, and Rroot-Ra ratio sites spatial distribution
plot_sites <- function(srdb_v4, Froot_data) {
  # prepare data for word map
  worldMap <- map_data(map = "world")

  # RC sites from SRDB_v4
  srdb_v4 %>%
    select(Study_number, Biome, Leaf_habit, Latitude, Longitude, RC_annual) %>%
    filter(RC_annual > 0 & RC_annual < 1) %>%
    # do what the 'summarySE' function used to
    filter(!is.na(Latitude)) %>%
    group_by(Latitude, Longitude) %>%
    summarise(
      N = n(),
      median = median(RC_annual),
      sd = sd(RC_annual),
      RC_annual = mean(RC_annual)
    ) %>%
    mutate(
      se = sd / sqrt(N),
      ci = se * qt(0.95 / 2 + 0.5, N - 1)
    ) ->
  siteInfor

  # prepare Ra-GPP sites data
  srdb_v4 %>%
    select(Ecosystem_type, Leaf_habit, Rs_annual, Ra_annual, Rh_annual, ER, GPP, NPP, Latitude, Longitude) %>%
    filter(!is.na(GPP) | !is.na(NPP) | !is.na(ER)) %>%
    mutate(
      ER = if_else(is.na(ER), GPP, ER), # assume NEP very small
      RaGPP_ratio = (ER - NPP) / GPP,
      Source = "srdb_v4"
    ) %>%
    filter(RaGPP_ratio > 0, RaGPP_ratio < 1, !is.na(Ecosystem_type)) -> ra_gpp_srdb

  # prepare data for legend
  cc_legend <- tibble(
    x = rep(-170, 3),
    y = c(-25, -40, -55),
    size = c(1.25, 1.25, 1.25)
  )

  # Base map - word map
  sitemap <- ggplot(data = worldMap) +
    # geom_polygon(aes(x = long, y = lat , fill = region , group = group, alpha = 0.1), color = "white") +
    geom_polygon(aes(x = long, y = lat, group = group, alpha = 0.1), color = "white", fill = "gray") +
    coord_fixed(1.3) +
    theme(
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.title.y = element_text(size = 13, margin = margin(t = 0, r = 12, b = 0, l = 0)),
      axis.title.x = element_text(size = 13, margin = margin(t = 12, r = 0, b = 0, l = 0)),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1.25)
    ) +
    theme(legend.position = "none") +
    scale_x_continuous(
      name = "Longitude", breaks = seq(-180, 180, 30),
      labels = seq(-180, 180, 30)
    ) +
    scale_y_continuous(
      name = "Latitude", limits = c(-60, 90), breaks = seq(-90, 90, 15),
      labels = seq(-90, 90, 15)
    ) +
    # Rroot-Rs site
    geom_point(
      data = siteInfor, aes(x = siteInfor$Longitude, y = siteInfor$Latitude),
      color = "black", shape = 1, size = 2, alpha = 0.75
    ) +

    # Rroot-to-Ra ratio sites
    geom_point(
      data = Froot_data, aes(Longitude, Latitude), color = "red",
      shape = 2, size = 3.5, alpha = 1
    ) +

    # Ra-GPP ratio sited
    geom_point(
      data = ra_gpp_srdb, aes(Longitude, Latitude), color = "blue",
      shape = 3, size = 5, alpha = 1
    ) +

    # legend
    geom_point(
      data = cc_legend, aes(x, y, size = size), shape = c(1, 2, 3),
      color = c("black", "red", "blue"), alpha = c(1, 1, 1)
    ) +
    annotate("text", x = -155, y = -25, label = expression(R[root] ~ ":" ~ R[S] ~ ratio), size = 3.5, hjust = 0) +
    annotate("text", x = -155, y = -40, label = expression(R[root] ~ ":" ~ R[A] ~ ratio), size = 3.5, hjust = 0) +
    annotate("text", x = -155, y = -55, label = expression(R[A] ~ ":" ~ GPP ~ ratio), size = 3.5, hjust = 0) +
    guides(fill = FALSE) # do this to leave off the color legend
  print(sitemap)
}


# Plot Rroot-Ra ratio and group by ecosystem
plot_Rroot_Ra_ratio <- function(Froot_data) {
  # plot Froot (root respiration/autotrophic respiration ratio) by ecosystem type
  obs_DF <- Froot_data %>%
    filter(IGBP2 == "DF") %>%
    nrow()
  obs_EF <- Froot_data %>%
    filter(IGBP2 == "EF") %>%
    nrow()
  obs_MF <- Froot_data %>%
    filter(IGBP2 == "MF") %>%
    nrow()
  obs_Other <- Froot_data %>%
    filter(IGBP2 == "Other") %>%
    nrow()
  Fl_plot <- ggplot(Froot_data, aes(x = IGBP2, y = Froot)) +
    geom_violin() +
    # geom_jitter(shape = 16, position = position_jitter(0.2), col = "gray") +
    geom_quasirandom(col = "gray") +
    geom_boxplot(width = 0.1) +
    ylab("Froot") +
    # stat_summary(fun.y=median, geom="point", size=2, color="red") +
    scale_x_discrete(
      limits = c("DF", "EF", "MF", "Other"),
      labels = c(
        paste0("DF (n=", obs_DF, ")"), paste0("EF (n=", obs_EF, ")"),
        paste0("MF (n=", obs_MF, ")"), paste0("Other (n=", obs_Other, ")")
      )
    ) +
    theme(axis.title.x = element_blank()) +
    ylab(expression(R[root] ~ ":" ~ R[A] ~ ratio))

  print(Fl_plot)
}

# plot Ra-GPP ratio and grouped by exosystem
plot_RaGPP <- function(RaGPP_data) {
  var_obs <- nrow(RaGPP_data)
  obs_DF <- RaGPP_data %>%
    filter(IGBP2 == "Deciduous") %>%
    nrow()
  obs_EF <- RaGPP_data %>%
    filter(IGBP2 == "Evergreen") %>%
    nrow()
  obs_MF <- RaGPP_data %>%
    filter(IGBP2 == "Mixed") %>%
    nrow()
  obs_GRA <- RaGPP_data %>%
    filter(IGBP2 == "Grassland") %>%
    nrow()
  obs_Other <- RaGPP_data %>%
    filter(IGBP2 == "Other") %>%
    nrow()

  RaGPP_data$Global <- paste0("n = ", var_obs)

  Ra_GPP <- ggplot(RaGPP_data, aes(IGBP2, RaGPP_ratio)) +
    geom_violin() +
    # geom_jitter(shape = 16, position = position_jitter(0.2), col = 'gray') +
    geom_quasirandom(col = "gray") +
    geom_boxplot(width = 0.1) +
    # stat_summary(fun.y = median, geom = "point", size = 2, color = "red") +
    ylab(expression(R[A] ~ ":" ~ GPP ~ ratio)) +
    scale_x_discrete(
      limits = c("Deciduous", "Evergreen", "Mixed", "Grassland", "Other"),
      labels = c(
        paste0("DF (n = ", obs_DF, ")"), paste0("EF (n = ", obs_EF, ")"),
        paste0("MF (n = ", obs_MF, ")"), paste0("Grassland (n = ", obs_GRA, ")"),
        paste0("Other (n = ", obs_Other, ")")
      )
    ) +
    theme(axis.title.x = element_blank())

  print(mean(RaGPP_data$RaGPP_ratio))
  se <- sd(RaGPP_data$RaGPP_ratio) / sqrt(nrow(RaGPP_data))
  paste("Ra/GPP se =", se) %>% print()
  CI <- round(qt(0.975, df = nrow(RaGPP_data) - 1) * se, 3)
  paste("Ra/GPP 95% CI =", CI) %>% print()
  paste("Ra/GPP obs =", nrow(RaGPP_data)) %>% print()

  # print figure
  print(Ra_GPP)
}


# plot Rroot-to-Rs ratio and group by ecosystem
plot_Rroot_Rs_NPP <- function(sub_srdb, NPP_data) {

  # plot NPP panel
  NPP_data$Global <- "n = 251"
  obs_AG <- sub_srdb %>%
    filter(IGBP2 == "Agriculture") %>%
    nrow()
  obs_DF <- sub_srdb %>%
    filter(IGBP2 == "DF") %>%
    nrow()
  obs_EF <- sub_srdb %>%
    filter(IGBP2 == "EF") %>%
    nrow()
  obs_MF <- sub_srdb %>%
    filter(IGBP2 == "Mixed") %>%
    nrow()
  obs_GRA <- sub_srdb %>%
    filter(IGBP2 == "Grassland") %>%
    nrow()
  obs_SHR <- sub_srdb %>%
    filter(IGBP2 == "Shrubland") %>%
    nrow()
  obs_Other <- sub_srdb %>%
    filter(IGBP2 == "Other") %>%
    nrow()

  p_NPP <- ggplot(NPP_data, aes(Global, NPP)) +
    geom_violin() +
    # geom_jitter(shape = 16, position = position_jitter(0.2), col = 'gray') +
    geom_quasirandom(col = "gray", varwidth = TRUE) +
    geom_boxplot(width = 0.1) +
    ylab(expression(NPP ~ (Pg ~ yr^{-1}))) +
    # stat_summary(fun.y = median, geom = "point", size = 2, color = "red") +
    theme(axis.title.x = element_blank())
  print(median(NPP_data$NPP))
  se <- sd(NPP_data$NPP) / sqrt(nrow(NPP_data))
  paste("NPP se =", se) %>% print()
  CI <- round(qt(0.975, df = nrow(NPP_data) - 1) * se, 3)
  paste("NPP 95% CI =", CI) %>% print()

  # plot Fire panel
  p_fire <- tibble(Fire = c(2, 3.5, 7.3, 4, 5.1, 2.02, 2.71, 3.02, 2.08)) %>%
    ggplot(aes(x = "n = 9", y = Fire)) +
    geom_violin() +
    # geom_jitter(shape = 16, position = position_jitter(0.2), col = 'gray') +
    geom_quasirandom(col = "gray", varwidth = TRUE) +
    geom_boxplot(width = 0.1) +
    ylab(expression(C[fire] ~ (Pg ~ yr^{-1}))) +
    # stat_summary(fun.y = median, geom = "point", size = 2, color = "red") +
    theme(axis.title.x = element_blank())

  p_mul <- plot_grid(p_NPP, p_fire,
    nrow = 1, labels = c("( a )", "( b )"),
    vjust = c(3), hjust = c(-2.8, -2.8, -2.5)
  )

  # plot Rroot_Rs ratio panel
  sub_srdb$Global <- paste0("n = ", nrow(sub_srdb))

  RC_plot <- ggplot(sub_srdb, aes(x = IGBP2, y = RC_annual)) +
    geom_violin() +
    # geom_jitter(shape = 16, position = position_jitter(0.2), col = 'gray') +
    geom_quasirandom(col = "gray", varwidth = TRUE) +
    geom_boxplot(width = 0.1) +
    scale_x_discrete(
      limits = c("Agriculture", "DF", "EF", "Mixed", "Grassland", "Shrubland", "Other"),
      labels = c(
        paste0("Agriculture (n = ", obs_AG, ")"), paste0("DF (n = ", obs_DF, ")"),
        paste0("EF (n = ", obs_EF, ")"), paste0("MF (n = ", obs_MF, ")"),
        paste0("GRA (n = ", obs_GRA, ")"), paste0("SHR (n = ", obs_SHR, ")"),
        paste0("Other (n = ", obs_Other, ")")
      )
    ) +
    # stat_summary(fun.y = median, geom = "point", size = 2, color = "red") +
    ylab(expression(R[root] ~ ":" ~ R[S] ~ ratio)) +
    theme(axis.title.x = element_blank())

  print(mean(sub_srdb$RC_annual))
  se <- sd(sub_srdb$RC_annual) / sqrt(nrow(sub_srdb))
  paste0("RC se = ", se) %>% print()
  CI <- round(qt(0.975, df = nrow(sub_srdb) - 1) * se, 3)
  paste("RC 95% CI =", CI) %>% print()
  paste("RC obs =", nrow(sub_srdb)) %>% print()

  # output figure
  plot_grid(p_mul, RC_plot,
    nrow = 2, labels = c("", "( c )"),
    vjust = c(3), hjust = c(-2.8)
  )
}



#*****************************************************************************************************************
# functions not used
#*****************************************************************************************************************
# boots function
# bootsting <- function (sdata) {
#   k = 10000
#   mysamples = replicate (k, sample(sdata, 2, replace = T))
#   mymeans = apply(mysamples, 2, mean)
#
#   boots_plot <- tibble(Boots = mymeans) %>%
#     ggplot(., aes(Boots)) + geom_histogram(color = "black", fill = "gray", bins = 30)
#
#   print(boots_plot)
#   return(mymeans)
# }


# bottom_up <- function () {
#   var_Rs <- GlobalRs %>% filter(!is.na(Rs)) %>% select(Rs)
#   samp_Rs <- sample(var_Rs$Rs, 1, replace = T)
#
#   # RC
#   var_RC <- srdb_v4 %>% select(RC_annual) %>% filter(RC_annual > 0 & RC_annual < 1 & !is.na(RC_annual))
#   # mean(var_Rroot_Rs_ratio$RC_annual)
#   samp_RC <- sample(var_RC$RC_annual, 1, replace = T)
#
#   # Froot and Froot/Fshoot ratio
#   var_Froot <- Froot %>% select (Froot)
#   samp_Froot <- sample(var_froot$Froot, 1, replace = T)
#   samp_FrFs_ratio <- 1/(1-sam_Froot) # Froot/Fshoot ratio
#
#   # sampling NPP
#   samp_NPP <- sample(NPP$NPP, 1, replace = T)
#
#   Rroot <- samp_Rs * samp_RC
#   Rshoot <- samp_Rs * samp_FrFs_ratio
#   samp_GPP <- samp_NPP + Rroot + Rshoot
# }
#
# result <- replicate(1000, bottom_up())
# summary(result)
# hist(result)
# abline(v = quantile(result, c(0.025, 0.975)), col = "red", lty = "dashed")

# plot GPP and RSG
# plot GPP, GPP vs published year relationship, and reported GPP trend
# plot RSG, RSG vs published year relationship, and reported RSG trend
# plot_GPP_RSG <- function (GPP_data, RSG_data) {
#
#   obs_gpp <- nrow(GPP_data)
#   GPP_data$Global <- paste0("n = ", obs_gpp)
#   p_GPP <- ggplot(GPP_data, aes(x = Global, y = GPP)) + geom_violin() +
#     geom_quasirandom(col = 'gray') +
#     # geom_jitter(shape = 16, position = position_jitter(0.2), col = 'gray') +
#     geom_boxplot(width = 0.1) +
#     # stat_summary(fun.y = median, geom = "point", size = 2, color = "red") +
#     ylab(expression(GPP~(Pg~yr^{-1})))
#
#   # GPP_His <- ggplot(data = GPP_data, aes(GPP)) +
#   # geom_histogram(bins = 30, col = "black", fill = "gray",alpha = 0.6) + #labs(title = "Histogram for GPP") +
#   # labs(x = "GPP", y = "Number of estimates") # xlim(c(70,185)) + # ylim(c(0,4))
#   # GPP_His <- GPP_His + annotate("text", x = 85, y = 6.5, label = "(a)" , size = 5)
#
#   # Trend
#   subGPP <- subset(GPP_data, GPP > 90 & GPP < 180 )
#   Trend_GPP <- ggplot(subGPP, aes(Pub_year, GPP)) +
#     geom_point(col = "gray") +
#     # geom_smooth(aes(Pub_year, GPP), method = "lm", col = "black", se = FALSE) +
#     geom_smooth(method = "lm", col = "red", data = GPP_data, se = TRUE) +
#     theme(legend.title = element_blank()) +
#     # one outlier
#     # geom_point(data = subset(GPP_data, GPP_data$GPP < 90 | GPP_data$GPP > 180) ,aes(x = Pub_year, y = GPP), color = c("red"), shape = 16 ) +
#     theme(legend.position = "none") +
#     theme(axis.title.y = element_blank()) +
#     labs(x = expression(Data~published~time~"(year)"))
#   # annotate("text", x = 1958, y = 175, size = 5)
#
#   lm_model1 <- lm(GPP ~ Pub_year, data = subGPP)
#   lm_model2 <- lm(GPP ~ Pub_year, data = GPP)
#
#   print(summary(lm_model1))
#   print(summary(lm_model2))
#
#   slope1 <- coefficients(summary(lm_model1))[2,1]
#   slope2 <- coefficients(summary(lm_model2))[2,1]
#
#   GPP_data$Global <- paste("n =", GPP %>% filter(!is.na(Trend)) %>% nrow())
#   p_trend <- ggplot(GPP_data, aes(Global, Trend)) + geom_violin() +
#     geom_quasirandom(col = 'gray') +
#     # geom_jitter(shape = 16, position = position_jitter(0.2), col = 'gray') +
#     geom_boxplot(width = 0.1) +
#     # geom_point(aes(y = slope1), col = "black", shape = 13) +
#     # geom_point(aes(y = slope2), col = "red", shape = 13) +
#     # stat_summary(fun.y = median, geom = "point", size = 2, color = "red") +
#     ylab(expression(GPP~trend~(Pg~yr^{-2})))
#
#   print(paste0("median = ", median(GPP_data$GPP)))
#   se <- sd(GPP_data$GPP) / sqrt(nrow(GPP_data))
#   paste0("se = ", se) %>% print()
#   CI <- round(qt(0.975, df = nrow(GPP_data)-1) * se, 3)
#   print(paste0("95% CI = ", CI))
#   print(paste0("Trend = ", nrow(GPP[!is.na(GPP$Trend),])))
#
#   # plot global Rs
#   p_Rs <- ggplot(RSG_data, aes(x = paste("n =", nrow(RSG_data)), y = Rs)) +
#     geom_violin() +
#     geom_quasirandom(col = 'gray') +
#     # geom_jitter(shape = 16, position = position_jitter(0.2), col = 'gray') +
#     geom_boxplot(width = 0.1) +
#     # stat_summary(fun.y = median, geom = "point", size = 2, color = "red") +
#     ylab(expression(R[SG]~(Pg~yr^{-1}))) +
#     xlab("Global")
#
#   # plot Rs~Pub_year relationship
#   sum_mod <- lm(Rs~Pub_year, data = RSG_data)
#   slope3 <- coefficients(summary(sum_mod))[2,1]
#
#   Trend_Rs <- ggplot(RSG_data, aes(Pub_year, Rs)) + geom_point(col = "gray") +
#     geom_smooth(method = "lm", col = "red") +
#     xlim(min(GPP_data$Pub_year), max(GPP_data$Pub_year)) +
#     theme(legend.title = element_blank(), axis.title.y = element_blank()) +
#     annotate("text", x = 1992, y = 80,
#              label = paste("Trend =", round(slope3, 3)), na.rm = TRUE) +
#     labs(x = expression(Data~published~time~"(year)"))
#   print(sum_mod)
#
#   # plot trend
#   Rs_IRate <- RSG_data %>%
#     filter(!is.na(IncreaseRate)) %>%
#     ggplot(aes(x = "n = 7", y = IncreaseRate)) +
#     geom_violin() +
#     geom_quasirandom(col = 'gray') +
#     # geom_jitter(shape = 16, position = position_jitter(0.2), col = 'gray') +
#     geom_boxplot(width = 0.1) +
#     ylab(expression(R[S]~trend~(Pg~yr^{-2}))) +
#     xlab("Global")
#   # geom_point(aes(y = slope3), col = "red", shape = 13)
#
#   plot_grid(p_Rs, Trend_Rs, Rs_IRate, p_GPP, Trend_GPP, p_trend,
#             ncol = 3,
#             labels = c('( a )', '( b )', '( c )', "( d )", "( e )", "( f )"),
#             vjust = c(3),
#             hjust = c(-2.35, -1.5, -2.25, -2.25, -1.5, -2.65))
# }


# This is the old joint figure, with sloping lines that everyone found unclear
# ```{r old-joint, results='show'}
# n_points <- 20
# # w is the weight of top-down GPP dominance, 0 to 1
# w <- seq(0, 1, length.out = 20)
# w <- c(0.25, 0.5, 0.75, w) %>%
#   unique() %>%
#   sort() # ensure 0.25, 0.5, and 0.75 are included
# 
# gpp_wm <- tibble(
#   q025 = rep(NA_real_, n_points),
#   q500 = rep(NA_real_, n_points),
#   q975 = rep(NA_real_, n_points)
# )
# rs_wm <- gpp_wm
# 
# for (i in seq_along(w)) {
#   weighted_gpp <- gpp_results$GPP_raw * w[i] + gpp_results$GPP2 * (1 - w[i])
#   gpp_wm[i, ] <- quantile(weighted_gpp, probs = QUANTILES)
#   weighted_rs <- rs_results$Rs_raw * (1 - w[i]) + rs_results$Rs_topdown2 * w[i]
#   rs_wm[i, ] <- quantile(weighted_rs, probs = QUANTILES)
# }
# 
# gpp_wm$w <- w
# rs_wm$w <- w
# results <- bind_rows("GPP" = gpp_wm, "Rs" = rs_wm, .id = "which")
# 
# # add col_scale column, which enable color gradual change in Figure 1 panel 3
# results %>% mutate(col_scale = 1:46) -> results
# 
# results_w75 <- filter(results, w == 0.75)
# results_w50 <- filter(results, w == 0.5)
# results_w25 <- filter(results, w == 0.25)
# myarrow <- arrow(length = unit(0.05, "inches"))
# 
# p_joint <- ggplot(results, aes(w, q500, group = which)) +
#   geom_vline(xintercept = 0.5, size = 0.25, linetype = 2) +
#   # scale_fill_manual(values=c("#999999", "#E69F00")) +
#   
#   geom_line(size = 1.5) +
#   geom_ribbon(aes(fill = which, ymin = q025, ymax = q975),
#               col = NA, alpha = 0.25, show.legend = FALSE
#   ) +
#   coord_cartesian(ylim = c(50, 170), expand = FALSE) +
#   scale_color_discrete("", labels = c("GPP", expression(R[S]))) +
#   theme(legend.text.align = 0) +
#   # scale_fill_manual(values=c("brown", "blueviolet")) +
#   
#   # 75%
#   geom_segment(
#     data = results_w75, aes(x = 0, xend = w, yend = q500),
#     linetype = 2, size = 0.5
#   ) +
#   geom_label(
#     data = results_w75, aes(x = 0.125, y = q500, label = round(q500, 0)),
#     size = 2, show.legend = FALSE
#   ) +
#   # 50%
#   geom_segment(
#     data = results_w50, aes(x = 0, xend = w, yend = q500),
#     linetype = 2, size = 1
#   ) +
#   geom_label(
#     data = results_w50, aes(x = 0.075, y = q500, label = round(q500, 0)),
#     size = 3.5, fontface = "bold", show.legend = FALSE
#   ) +
#   
#   # 25%
#   geom_segment(
#     data = results_w25, aes(x = 0, xend = w, yend = q500),
#     linetype = 2, size = 0.5
#   ) +
#   geom_label(
#     data = results_w25, aes(x = 0.025, y = q500, label = round(q500, 0)),
#     size = 2, show.legend = FALSE
#   ) +
#   
#   # Arrows and labels
#   annotate("segment", x = 0.1, xend = 0.025, y = 60, yend = 60, arrow = myarrow) +
#   annotate("text",
#            x = 0.11, y = 60, size = 3, hjust = "inward",
#            label = expression(bold(Determined ~ more ~ by ~ R[S]))
#   ) +
#   annotate("segment", x = 0.9, xend = 0.975, y = 60, yend = 60, arrow = myarrow) +
#   annotate("text",
#            x = 0.89, y = 60, size = 3, hjust = "inward", fontface = "bold",
#            label = "Determined more by GPP"
#   ) +
#   
#   # Other
#   xlab("Weight of GPP") +
#   ylab(expression(Flux ~ (Pg ~ C ~ yr^-1)))
# 
# # change color
# # scale_color_manual(values=c("#999999", "#E69F00"))
# 
# print(p_joint)
# knitr::kable(results_w50, digits = 1, format = "markdown")
# ```
