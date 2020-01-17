# Helper functions

#******************************************************************************************************************
# data processing function - data calculation, clean, and combine etc
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

# Function to calculate Rroot/Ra ratio and Rshoot/Ra ratio
# Ra = Rroot + Rshoot (total autotrophic respiration)
# Rroot/Ra or Rshoot/Ra ratio were from two resources: 1) collected from published article, stored in FsFlFr.csv; 2) from srdb_v4.csv
# and results from above two data sources will be combined
cal_Froot <- function (sdata, sdata2) {
  # Froot and Fshoot data from digitized papers
  sdata %>%
    select(Latitude, Longitude, Fshoot, Froot, Ecosystem, IGBP) %>% 
    mutate(Fshoot = if_else(is.na(Fshoot), 100 - Froot, as.numeric(Fshoot)),
           Froot = if_else(is.na(Froot), 100 - Fshoot, as.numeric(Froot)),
           # data in FsFlFr were in percentage, scale to decimal (to match srdb_v4))
           Froot = Froot / 100,
           Fshoot = Fshoot / 100) ->
    comb_data
  
  # use data from srdb_v4 to calculate Rroot/Ra ratio, Rroot/(RE-NPP)
  # Ra_annual in srdb_v4 is Rroot
  # Ra = Rroot + Rshoot = ER - Rh
  sdata2 %>% 
    # select(Study_number, Biome, Ecosystem_type, Leaf_habit, Latitude, Longitude, Ra_annual, ER, GPP, NPP) %>% 
    select(Latitude, Longitude, Ecosystem_type, Leaf_habit, Rs_annual, Ra_annual, Rh_annual, ER, GPP, NPP) %>% 
    filter(!is.na(GPP) | !is.na(NPP) | !is.na(ER)) %>% 
    mutate(ER = if_else(is.na(ER), GPP, ER), # If ER not available, use GPP (assume NEP very small)
           # If Rh_annual not available, use NPP (assume NEP very small)
           Rh_annual = if_else(is.na(Rh_annual), NPP, Rh_annual),
           Ra_annual = if_else(is.na(Ra_annual), Rs_annual - Rh_annual, Ra_annual)) -> 
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
  comb_data[comb_data$Ecosystem  ==  'Cropland', ]$Ecosystem <- "Agriculture"
  return (comb_data)
}


# Prepare and calculate Ra-GPP ratio 
# data are from two sources: 1) from published papers, stored in RaGPP.csv; 2) from srdb_v4.csv
# data from above two sources will be combined
prep_RaGpp <- function(RaGPP, srdb_v4) {
  RaGPP %>% 
    select(Ecosystem, Leaf_habit, RaGPP_ratio) %>% 
    mutate(Source = "Meta-papers") -> sdata2
  
  srdb_v4 %>% 
    select(Ecosystem_type, Leaf_habit, Rs_annual, Ra_annual, Rh_annual, ER, GPP, NPP) %>% 
    filter(!is.na(GPP) | !is.na(NPP) | !is.na(ER)) %>% 
    mutate(ER = if_else(is.na(ER), GPP, ER), # assume NEP very small
           RaGPP_ratio = (ER - NPP) / GPP,
           Source = "srdb_v4") %>% 
    filter(RaGPP_ratio > 0, RaGPP_ratio < 1, !is.na(Ecosystem_type)) %>%
    select(Ecosystem_type, Leaf_habit, RaGPP_ratio, Source) %>% 
    bind_rows(sdata2)
}


#******************************************************************************************************************
# functions to make figures
#******************************************************************************************************************
# plot RC (root respiration to soil respiration ratio), Ra-GPP ratio, and Rroot-Ra ratio sites spatial distribution
plot_sites <- function (sdata, sdata2, srdb_v4) {
  # prepare data for word map
  worldMap <- map_data(map = "world")
  
  # RC sites from SRDB_v4
  sdata %>% 
    select(Study_number, Biome, Leaf_habit, Latitude, Longitude, RC_annual) %>% 
    filter(RC_annual > 0 & RC_annual < 1) %>% 
    # do what the 'summarySE' function used to
    filter(!is.na(Latitude)) %>% 
    group_by(Latitude, Longitude) %>% 
    summarise(N = n(),
              median = median(RC_annual),
              sd = sd(RC_annual),
              RC_annual = mean(RC_annual)) %>% 
    mutate(se = sd / sqrt(N),
           ci = se * qt(0.95 / 2 + 0.5, N - 1)) ->
    siteInfor
  
  # prepare Ra-GPP sites data
  srdb_v4 %>% 
    select(Ecosystem_type, Leaf_habit, Rs_annual, Ra_annual, Rh_annual, ER, GPP, NPP, Latitude, Longitude) %>% 
    filter(!is.na(GPP) | !is.na(NPP) | !is.na(ER)) %>% 
    mutate(ER = if_else(is.na(ER), GPP, ER), # assume NEP very small
           RaGPP_ratio = (ER - NPP) / GPP,
           Source = "srdb_v4") %>% 
    filter(RaGPP_ratio > 0, RaGPP_ratio < 1, !is.na(Ecosystem_type)) -> ra_gpp_srdb
    
  # prepare data for legend
  cc_legend <- tibble(x = rep(-170, 3), 
                      y = c(-25, -40, -55),
                      size = c(1.25, 1.25, 1.25))
  
  # Base map - word map
  sitemap <- ggplot(data = worldMap) + 
    # geom_polygon(aes(x = long, y = lat , fill = region , group = group, alpha = 0.1), color = "white") + 
    geom_polygon(aes(x = long, y = lat, group = group, alpha = 0.1), color = "white", fill = "gray") + 
    coord_fixed(1.3) +
    theme(axis.text.y   = element_text(size = 12),
          axis.text.x   = element_text(size = 12),
          axis.title.y   = element_text(size = 13, margin = margin(t = 0, r = 12, b = 0, l = 0)),
          axis.title.x   = element_text(size = 13, margin = margin(t = 12, r = 0, b = 0, l = 0)),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, size = 1.25))+
    theme(legend.position = "none")+
    scale_x_continuous(name = "Longitude", breaks = seq(-180, 180, 30),
                       labels = seq(-180, 180, 30)) +
    scale_y_continuous(name = "Latitude", limits = c(-60, 90), breaks = seq(-90, 90, 15),
                       labels = seq(-90,90,15))+
    geom_point(data = siteInfor, aes(x = siteInfor$Longitude, y = siteInfor$Latitude), 
               color = "black", shape = 1, size = 2, alpha = 0.75) + 
    
    # Rroot-to-Ra ratio sites
    geom_point(data = sdata2, aes(Longitude, Latitude), color = "red",
               shape = 2, size = 3.5, alpha = 1) +
    
    # Ra-GPP ratio sited
    geom_point(data = ra_gpp_srdb, aes(Longitude, Latitude), color = "blue",
               shape = 3, size = 5, alpha = 1) +
    
    # legend
    geom_point(data = cc_legend, aes(x, y, size = size), shape = c(1, 2, 3), 
               color = c("black", "red", "blue"), alpha = c(1, 1, 1)) +
    annotate("text", x = -150, y = -25, label = expression("Rroot-"~R[S]~ratio), size = 4, hjust = 0) +
    annotate("text", x = -150, y = -40, label = "Rroot-Ra ratio", size = 4, hjust = 0) +
    annotate("text", x = -150, y = -55, label = "Ra-GPP ratio", size = 4, hjust = 0) +
    guides(fill = FALSE)  # do this to leave off the color legend
  print(sitemap)
}


# plot GPP and RSG 
# plot GPP, GPP vs published year relationship, and reported GPP trend
# plot RSG, RSG vs published year relationship, and reported RSG trend
plot_GPP_RSG <- function (sdata, sdata2) {
  
  obs_gpp <- nrow(sdata)
  sdata$Global <- paste0("n = ", obs_gpp)
  p_GPP <- ggplot(sdata, aes(x = Global, y = GPP)) + geom_violin() +
    geom_quasirandom(col = 'gray') +
    # geom_jitter(shape = 16, position = position_jitter(0.2), col = 'gray') +
    geom_boxplot(width = 0.1) +
    # stat_summary(fun.y = median, geom = "point", size = 2, color = "red") +
    ylab(expression(GPP~(Pg~yr^{-1})))
  
  # GPP_His <- ggplot(data = sdata, aes(GPP)) + 
  # geom_histogram(bins = 30, col = "black", fill = "gray",alpha = 0.6) + #labs(title = "Histogram for GPP") +
  # labs(x = "GPP", y = "Number of estimates") # xlim(c(70,185)) + # ylim(c(0,4)) 
  # GPP_His <- GPP_His + annotate("text", x = 85, y = 6.5, label = "(a)" , size = 5)
  
  # Trend
  subGPP <- subset(sdata, GPP > 90 & GPP < 180 )
  Trend_GPP <- ggplot(subGPP, aes(Pub_year, GPP)) +
    geom_point(col = "gray") + 
    # geom_smooth(aes(Pub_year, GPP), method = "lm", col = "black", se = FALSE) +
    geom_smooth(method = "lm", col = "red", data = sdata, se = TRUE) +
    theme(legend.title = element_blank()) +
    # one outlier
    # geom_point(data = subset(sdata, sdata$GPP < 90 | sdata$GPP > 180) ,aes(x = Pub_year, y = GPP), color = c("red"), shape = 16 ) +
    theme(legend.position = "none") +
    theme(axis.title.y = element_blank()) +
    labs(x = expression(Data~published~time~"(year)"))
  # annotate("text", x = 1958, y = 175, size = 5)
  
  lm_model1 <- lm(GPP ~ Pub_year, data = subGPP)
  lm_model2 <- lm(GPP ~ Pub_year, data = GPP)
  
  print(summary(lm_model1))
  print(summary(lm_model2))
  
  slope1 <- coefficients(summary(lm_model1))[2,1]
  slope2 <- coefficients(summary(lm_model2))[2,1]
  
  sdata$Global <- paste("n =", GPP %>% filter(!is.na(Trend)) %>% nrow())
  p_trend <- ggplot(sdata, aes(Global, Trend)) + geom_violin() +
    geom_quasirandom(col = 'gray') +
    # geom_jitter(shape = 16, position = position_jitter(0.2), col = 'gray') +
    geom_boxplot(width = 0.1) +
    # geom_point(aes(y = slope1), col = "black", shape = 13) + 
    # geom_point(aes(y = slope2), col = "red", shape = 13) +
    # stat_summary(fun.y = median, geom = "point", size = 2, color = "red") +
    ylab(expression(GPP~trend~(Pg~yr^{-2})))
  
  print(paste0("median = ", median(sdata$GPP)))
  se <- sd(sdata$GPP) / sqrt(nrow(sdata))
  paste0("se = ", se) %>% print()
  CI <- round(qt(0.975, df = nrow(sdata)-1) * se, 3)
  print(paste0("95% CI = ", CI))
  print(paste0("Trend = ", nrow(GPP[!is.na(GPP$Trend),])))
  
  # plot global Rs
  p_Rs <- ggplot(sdata2, aes(x = paste("n =", nrow(sdata2)), y = Rs)) + 
    geom_violin() +
    geom_quasirandom(col = 'gray') +
    # geom_jitter(shape = 16, position = position_jitter(0.2), col = 'gray') +
    geom_boxplot(width = 0.1) +
    # stat_summary(fun.y = median, geom = "point", size = 2, color = "red") +
    ylab(expression(R[SG]~(Pg~yr^{-1}))) +
    xlab("Global")
  
  # plot Rs~Pub_year relationship
  sum_mod <- lm(Rs~Pub_year, data = sdata2)
  slope3 <- coefficients(summary(sum_mod))[2,1]

  Trend_Rs <- ggplot(sdata2, aes(Pub_year, Rs)) + geom_point(col = "gray") +
    geom_smooth(method = "lm", col = "red") +
    xlim(min(sdata$Pub_year), max(sdata$Pub_year)) +
    theme(legend.title = element_blank(), axis.title.y = element_blank()) +
    annotate("text", x = 1992, y = 80, 
             label = paste("Trend =", round(slope3, 3)), na.rm = TRUE) +
    labs(x = expression(Data~published~time~"(year)"))
  print(sum_mod)
  
  # plot trend
  Rs_IRate <- sdata2 %>% 
    filter(!is.na(IncreaseRate)) %>% 
    ggplot(aes(x = "n = 7", y = IncreaseRate)) + 
    geom_violin() +
    geom_quasirandom(col = 'gray') +
    # geom_jitter(shape = 16, position = position_jitter(0.2), col = 'gray') +
    geom_boxplot(width = 0.1) +
    ylab(expression(R[S]~trend~(Pg~yr^{-2}))) +
    xlab("Global") 
    # geom_point(aes(y = slope3), col = "red", shape = 13)
  
  plot_grid(p_Rs, Trend_Rs, Rs_IRate, p_GPP, Trend_GPP, p_trend, 
            ncol = 3, 
            labels = c('( a )', '( b )', '( c )', "( d )", "( e )", "( f )"), 
            vjust = c(3), 
            hjust = c(-2.35, -1.5, -2.25, -2.25, -1.5, -2.65))
}


# Plot Rroot-Ra ratio and group by ecosystem
plot_Rroot_Ra_ratio <- function (sdata) {
  # plot Froot (root respiration/autotrophic respiration ratio) by ecosystem type
  Fl_plot <- ggplot(sdata, aes(x = IGBP2, y=Froot)) + geom_violin() +
    # geom_jitter(shape = 16, position = position_jitter(0.2), col = "gray") +
    geom_quasirandom(col = 'gray') +
    geom_boxplot(width = 0.1) +
    ylab("Froot") +
    # stat_summary(fun.y=median, geom="point", size=2, color="red") +
    scale_x_discrete(limits = c("DF", "EF", "MF", "Other"),
                     labels = c("DF (n=13)","EF (n=91)", "MF (n=7)", "Other (n=11)") ) +
    theme(axis.title.x=element_blank()) +
    ylab("Rroot-Ra ratio")
  
  print(Fl_plot)
}

# plot Ra-GPP ratio and grouped by exosystem
plot_RaGPP <- function (sdata2) {
  
  var_obs <- nrow(sdata2)
  
  sdata2$Global <- paste0("n = ", var_obs)
  
  Ra_GPP <- ggplot(sdata2, aes(IGBP2, RaGPP_ratio)) + geom_violin() +
    # geom_jitter(shape = 16, position = position_jitter(0.2), col = 'gray') +
    geom_quasirandom(col = 'gray') +
    geom_boxplot(width = 0.1) +
    # stat_summary(fun.y = median, geom = "point", size = 2, color = "red") +
    ylab("Ra-GPP ratio") +
    scale_x_discrete(limits = c("Deciduous", "Evergreen", "Forest", "Grassland", "Other"),
                     labels = c("Deciduous (n = 20)", "Evergreen (n = 127)", "Forest (n = 68)", "Grassland (n = 11)", "Other (n = 14)") ) +
    theme(axis.title.x = element_blank())
  
  print(mean(sdata2$RaGPP_ratio))
  se <- sd(sdata2$RaGPP_ratio) / sqrt(nrow(sdata2))
  paste("Ra/GPP se =", se) %>% print()
  CI <- round(qt(0.975, df = nrow(sdata2) - 1) * se, 3)
  paste("Ra/GPP 95% CI =", CI) %>% print()
  paste("Ra/GPP obs =", nrow(sdata2)) %>% print()
  
  # print figure
  print(Ra_GPP)
}



# plot Rroot-to-Rs ratio and group by ecosystem
plot_Rroot_Rs_NPP <- function (sdata, sdata2, sdata3) {
  sub_Rs <- subset(sdata, !is.na(sdata$Rs))
  Rs_obs <- nrow(sub_Rs)
  sub_Rs$label <- paste0("n = ", Rs_obs)
  p_Rs <- ggplot(sub_Rs, aes(x = label, y = Rs)) + geom_violin() +
    # geom_jitter(shape = 16, position = position_jitter(0.2), col = 'gray') +
    geom_quasirandom(col = 'gray', varwidth = TRUE) +
    geom_boxplot(width = 0.1) +
    # geom_point(aes(y = 71.67), col = "red") + geom_point(aes(y = 76.68), col = "red") +
    xlab("Global Rs") + ylab(expression(Rs~(Pg~yr^{-1}))) +
    theme(axis.title.x = element_blank())
  
  sub_Rh <- subset(sdata, !is.na(sdata$Rh))
  sub_Rh$label <- "n = 2"
  p_Rh <- ggplot(sub_Rh, aes(x = label, y = Rh)) + geom_boxplot(width = .1) +
    geom_point(colour = "gray", size = 2) +
    geom_point(aes(y = 46.47), col = "red") +
    xlab("Global Rh") + ylab(expression(Rh~(Pg~yr^{-1}))) +
    theme(axis.title.x = element_blank())
  
  sub_Ra <- subset(sdata, !is.na(sdata$Ra))
  sub_Ra$label <- "n = 1"
  p_Ra <- ggplot(sub_Ra, aes(x = label, y = Ra)) + geom_boxplot(width = 0.1) +
    geom_point(colour = "gray", size = 2) +
    geom_point(y = 25.20, col = "red") + geom_point(y = 30.21, col = "red") +
    xlab("Global Ra") + ylab(expression(Ra~(Pg~yr^{-1}))) +
    theme(axis.title.x = element_blank())
  
  # RC and RH distribution
  # RC
  sub_srdb %>%
    count(IGBP2) %>%
    print()
  
  #sub_srdb %>% group_by(Ecosystem_type) %>% summarise(count = n()) %>% print ()
  
  RC_obs <- nrow(sub_srdb)
  # sub_srdb %>% filter(Ra_Rh_Ratio < 5) -> sub_srdb
  sub_srdb$Global <- paste0("n = ", RC_obs)
  
  RC_plot <- ggplot(sub_srdb, aes(x = Global, y = RC_annual)) + geom_violin() +
    # geom_jitter(shape = 16, position = position_jitter(0.2), col = 'gray') +
    geom_quasirandom(col = 'gray', varwidth = TRUE) +
    geom_boxplot(width = 0.1) +
    # stat_summary(fun.y = median, geom = "point", size = 2, color = "red") +
    ylab("RC") +
    theme(axis.title.x = element_blank())
  
  # RC by biome
  RC_plot2 <- ggplot(sub_srdb, aes(x = IGBP2, y = RC_annual)) + geom_violin() +
    # geom_jitter(shape = 16, position = position_jitter(0.2), col = 'gray') +
    geom_quasirandom(col = 'gray', varwidth = TRUE) +
    geom_boxplot(width = 0.1) +
    scale_x_discrete(limits = c("Agriculture", "DF", "EF", "Mixed", "Grassland", "Shrubland", "Other"),
                     labels = c("Agriculture (n = 40)", "DF (n = 189)", "EF (n = 256)", "Mixed (n = 28)", "Grassland (n = 74)", "Shrubland (n = 14)", "Other (n = 16)")) +
    # stat_summary(fun.y = median, geom = "point", size = 2, color = "red") +
    ylab("Rroot-R"[S]~"ratio") +
    theme(axis.title.x = element_blank())
  
  
  print(mean(sub_srdb$RC_annual))
  se <- sd(sub_srdb$RC_annual) / sqrt(nrow(sub_srdb))
  paste0("RC se = ", se) %>% print()
  CI <- round(qt(0.975, df = nrow(sub_srdb)-1) * se, 3)
  paste("RC 95% CI =", CI) %>% print()
  paste("RC obs =", nrow(sub_srdb)) %>% print()
  
  # plot NPP
  sdata3$Global <- "n = 251"
  p_NPP <- ggplot(sdata3, aes(Global, NPP)) + geom_violin() +
    # geom_jitter(shape = 16, position = position_jitter(0.2), col = 'gray') +
    geom_quasirandom(col = 'gray', varwidth = TRUE) +
    geom_boxplot(width = 0.1) +
    ylab(expression(NPP~(Pg~yr^{-1}))) +
    # stat_summary(fun.y = median, geom = "point", size = 2, color = "red") +
    theme(axis.title.x = element_blank())
  print(median(sdata3$NPP))
  se <- sd(sdata3$NPP) / sqrt(nrow(sdata3))
  paste("NPP se =", se) %>% print()
  CI <- round(qt(0.975, df = nrow(sdata3) - 1) * se, 3)
  paste("NPP 95% CI =", CI) %>% print()
  
  # plot Fire
  p_fire <- tibble(Fire = c(2, 3.5, 7.3, 4, 5.1, 2.02, 2.71, 3.02, 2.08)) %>% 
    ggplot(aes(x = "Fire", y = Fire)) + geom_violin() +
    # geom_jitter(shape = 16, position = position_jitter(0.2), col = 'gray') +
    geom_quasirandom(col = 'gray', varwidth = TRUE) +
    geom_boxplot(width = 0.1) +
    ylab(expression(Fire~burned~(Pg~yr^{-1}))) +
    # stat_summary(fun.y = median, geom = "point", size = 2, color = "red") +
    theme(axis.title.x = element_blank())
  
  p_mul <- plot_grid(p_NPP, p_fire, 
                     nrow = 1, labels = c('( a )', '( b )'), 
                     vjust = c(3), hjust = c(-2.8, -2.8, -2.5)) 
  
  plot_grid(p_mul, RC_plot2, 
            nrow = 2, labels = c('', '( c )'),
            vjust = c(3), hjust = c(-2.8))
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
  
  ann_dat <- tibble(Type = types,
                    GPP_type = NA, Rs_type = NA,
                    x = xpos, y = c(0.035, 0.025),
                    label = c(make_label(x_qt), make_label(x_qt_raw)))
  figure + geom_text(data = ann_dat,
                     aes(x = x, y = y, label = label, color = Type),
                     hjust = 0, show.legend = FALSE)
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


