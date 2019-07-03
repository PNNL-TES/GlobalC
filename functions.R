
#*****************************************************************************************************************
# Baasic functions
#*****************************************************************************************************************
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=F,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     median = median   (xx[[col]], na.rm=na.rm),
                     #do.call("rbind", tapply(xx[[col]], measurevar, quantile, c(0.25, 0.5, 0.75)))
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#************************************************************
# plot Figure 1, sites distribution

plot_sites <- function (sdata, sdata2) {
  worldMap <- map_data(map="world")
  
  # RC sites from SRDB_v4
  sdata %>% 
    select(Study_number, Biome, Leaf_habit, Latitude, Longitude, RC_annual) %>% 
    filter(RC_annual > 0 & RC_annual < 1) -> 
    sub_srdb 
  # sub_srdb$RC_RH_Ratio <- sub_srdb$RC_annual/(1-sub_srdb$RC_annual)
  
  x <- rep(-170, 2)
  size <- c(2, 1.5)
  y <- c(-25, -40)
  
  cc_legend <- tibble (x, y, size)
  
  siteInfor <- summarySE (data=sub_srdb, measurevar='RC_annual', groupvars=c("Latitude","Longitude"))
  siteInfor <- siteInfor[, c(1:3)]
  siteInfor <- siteInfor[which(!is.na(siteInfor$Latitude)),]
 
  sitemap <- ggplot(data = worldMap) + 
    geom_polygon(aes(x = long, y = lat , fill = region , group = group, alpha = 0.1), color = "white") + 
    coord_fixed(1.3) +
    theme(axis.text.y   = element_text(size=12),
          axis.text.x   = element_text(size=12),
          axis.title.y   = element_text(size=13, margin = margin(t = 0, r = 12, b = 0, l = 0)),
          axis.title.x   = element_text(size=13, margin = margin(t = 12, r = 0, b = 0, l = 0)),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1.25))+
    theme(legend.position="none")+
    scale_x_continuous(name="Longitude", breaks=seq(-180,180, 30),labels = seq(-180,180, 30))+
    scale_y_continuous(name="Latitude", limits = c(-60,90),breaks=seq(-90,90,15),labels = seq(-90,90,15))+
    geom_point(data = siteInfor, aes(x=siteInfor$Longitude, y=siteInfor$Latitude), color = "black"
               ,shape = 18, size = 2, alpha = 1) + 
    
    
    # Rshoot/Rroot sites
    geom_point(data = sdata2, aes(x=sdata2$Longitude, y=sdata2$Latitude), color = "red"
               ,shape = 3, size = 1.5, stroke = 1.5, alpha = 1, fill = "red") +
    # legend
    geom_point(data = cc_legend ,aes(x=cc_legend$x, y=cc_legend$y), shape = c(18, 3), stroke = c(1, 1.5) 
               , color = c("black", "red"), size = cc_legend$size, alpha = c(1, 1)) +
    annotate("text", x = -150, y = -10, label = "Legend", size = 4, hjust = 0) +
    annotate("text", x = -150, y = -25, label = "RC sites", size = 4, hjust = 0) +
    annotate("text", x = -150, y = -40, label = "Rroot/Rshoot sites", size = 4, hjust = 0) +
    guides(fill=FALSE)  # do this to leave off the color legend
  print(sitemap)
}


# plot GPP
plot_GPP <- function (sdata) {
  
  sdata$Global <- "n=50"
  p_GPP <- ggplot(sdata, aes(x = Global, y=GPP)) + geom_violin() +
    geom_jitter(shape=16, position=position_jitter(0.2), col = 'gray') +
    geom_boxplot(width=.1) +
    # stat_summary(fun.y=median, geom="point", size=2, color="red") +
    ylab(expression(GPP~"("~Pg~yr^{-1}~")"))
  
  # GPP_His <- ggplot(data=sdata, aes(GPP)) + 
    # geom_histogram(bins = 30, col="black", fill="gray",alpha = 0.6) + #labs(title="Histogram for GPP") +
    # labs(x="GPP", y="Number of estimates") # xlim(c(70,185)) + # ylim(c(0,4)) 
  # GPP_His <- GPP_His + annotate("text", x = 85, y = 6.5, label = "(a)" , size = 5)

  # Trend
  subGPP <- subset( sdata, sdata$GPP > 90 & sdata$GPP < 180 )
  Trend_GPP <- ggplot (subGPP, aes(x = Year, y = GPP)) + geom_point() 
  Trend_GPP <- Trend_GPP + geom_smooth(aes(Year, GPP), method = "lm") +
    theme(legend.title=element_blank()) +
    # one outlier
    geom_point(data = subset(sdata, sdata$GPP < 90 | sdata$GPP > 180) ,aes(x=Year, y=GPP), color = c("red"), shape = 4 ) +
    theme(legend.position = "none") +
    theme(axis.title.y = element_blank())
    # annotate("text", x = 1958, y = 175, size = 5)
  
  lm_model1 <- lm(GPP~Year, data = subGPP)
  lm_model2 <- lm(GPP~Year, data = GPP)
  
  sdata$Global <- "n=16"
  p_trend <- ggplot(sdata, aes(x = Global, y=Trend)) + geom_violin() +
    geom_jitter(shape=16, position=position_jitter(0.2), col = 'gray') +
    geom_boxplot(width=.1) +
    geom_point(aes(y=0.2), col="red") + geom_point(aes(y=0.43), col="red") +
    # stat_summary(fun.y=median, geom="point", size=2, color="red") +
    ylab(expression(GPP~trend~"("~Pg~yr^{-2}~")"))
  
  print(paste0("median=", median(sdata$GPP)))
  se <- sd(sdata$GPP)/sqrt(nrow(sdata))
  paste0("se=", se) %>% print()
  CI <- round(qt(0.975,df=nrow(sdata)-1)*se, 3)
  print(paste0("95% CI=", CI))
  print(paste0("number of trend=", nrow(GPP[!is.na(GPP$Trend),])))
  print(summary(lm_model1))
  print(summary(lm_model2))
  
  plot_grid(p_GPP, Trend_GPP, p_trend, ncol = 3
            , labels = c('( a )', '( b )', '( c )')
            , vjust = c(3,3, 3), hjust = c(-2.35, -1.5, -2.25))

}



# plot NPP and Rab/Rh
plot_RaGPP <- function (sdata2, sdata3) {
  # Ra/GPP
  sdata2$Source <- "Meta-papers"
  sdata2 %>% select(Ecosystem, RaGPP_ratio, Source) -> sdata2
  
  
  sdata3 %>% 
    # select(Study_number, Biome, Ecosystem_type, Leaf_habit, Latitude, Longitude, Ra_annual, ER, GPP, NPP) %>% 
    select(Ecosystem_type, Rs_annual, Ra_annual, Rh_annual, ER, GPP, NPP) %>% 
    filter(!is.na(GPP) | !is.na(NPP) | !is.na(ER) ) -> sub_srdb
  sub_srdb$ER <- ifelse(is.na(sub_srdb$ER), sub_srdb$GPP, sub_srdb$ER) # assume NEP very small
  sub_srdb$RaGPP_ratio <- (sub_srdb$ER - sub_srdb$NPP)/sub_srdb$GPP
  sub_srdb %>% filter (RaGPP_ratio > 0 & RaGPP_ratio <1 & !is.na(Ecosystem_type)) -> sub_srdb
  sub_srdb$Source <- "srdb_v4"
  sub_srdb %>% select(Ecosystem_type, RaGPP_ratio, Source) -> sub_srdb
  colnames(sub_srdb) <- colnames(sdata2)
  
  sdata2 <- rbind(sdata2, sub_srdb)
  sdata2[sdata2$Ecosystem == "Crop", ]$Ecosystem <- "Agriculture"
  sdata2[sdata2$Ecosystem == "Tundra", ]$Ecosystem <- "Grassland" # only one tundra site, and IGBP has no tundra class, put it to grassland
  var_obs <- nrow(sdata2)
 
  sdata2$Global <- paste0("n=", var_obs)
  
  Ra_GPP <- ggplot(sdata2, aes(x = Ecosystem, y=RaGPP_ratio)) + geom_violin() +
    geom_jitter(shape=16, position=position_jitter(0.2), col = 'gray') +
    geom_boxplot(width=.1) +
    # stat_summary(fun.y=median, geom="point", size=2, color="red") +
    ylab("Ra/GPP ratio") +
    scale_x_discrete(limits = c("Agriculture", "Forest", "Grassland", "Savanna", "Wetland"),
                     labels = c("Agriculture (n=6)", "Forest (n=207)", "Grassland (n=7)", "Savanna (n=5)", "Wetland (n=7)") ) +
    theme(axis.title.x=element_blank())
  
  print(mean(sdata2$RaGPP_ratio))
  se <- sd(sdata2$RaGPP_ratio)/sqrt(nrow(sdata2))
  paste0("Ra/GPP se=", se) %>% print()
  CI <- round(qt(0.975,df=nrow(sdata2)-1)*se, 3)
  paste0("Ra/GPP 95% CI=", CI) %>% print()
  paste0("Ra/GPP obs=", nrow(sdata2)) %>% print()
  
  # print figure
  print(Ra_GPP)
  # plot_grid(p_NPP, Ra_GPP, labels = c('( a )', '( b )')
  #           , vjust = c(3.5,3.5), hjust = c(-2.5, -2.25) )
  
  return (sdata2)
}



# FlFsFr %>%
#   select(Fshoot, Froot, Ecosystem) ->
#   comb_data
# comb_data$Fshoot <- ifelse(is.na(comb_data$Fshoot), 100-comb_data$Froot, comb_data$Fshoot)
# comb_data$Froot <- ifelse(is.na(comb_data$Froot), 100-comb_data$Fshoot, comb_data$Froot)
# gather(comb_data, key = "Fraction", value = "Percentage", -Ecosystem) ->
#   comb_data


# Function to calculate Rroot/Ra ratio and Rshoot/Ra ratio
# Ra = Rroot + Rshoot (total autotrophic respiration)
# Rrrot/Ra or Rshoot/Ra ratio were from two resources: 1) collected from published article, stored in FsFlFr (sdata); 2) from srdb_v4 (sdata2)
#install.packages("reshapes") # Old package replaced by gather
# library(reshape2)
# FsFfFrReshape <- melt(FsFfFr,id.vars='IGBP', measure.vars=c('Fl','Fs','Fr'))
cal_Froot <- function (sdata, sdata2) {
  # Froot and Fshoot data from digitized papers
  sdata %>%
    select(Fshoot, Froot, Ecosystem) ->
    comb_data
  comb_data$Fshoot <- ifelse(is.na(comb_data$Fshoot), 100-comb_data$Froot, comb_data$Fshoot) 
  comb_data$Froot <- ifelse(is.na(comb_data$Froot), 100-comb_data$Fshoot, comb_data$Froot)
  comb_data$Froot <- comb_data$Froot/100 # data in FsFlFr were in percentage, scale to decimal (to match srdb_v4)
  comb_data$Fshoot <- comb_data$Fshoot/100 # data in FsFlFr were in percentage, scale to decimal (to match srdb_v4)
  
  # use data from srdb_v4 to calculate Rroot/Ra ratio, Rroot/(RE-NPP)
  # Ra_annual in srdb_v4 is Rroot
  # Ra = Rroot + Rshoot = ER - Rh
  sdata2 %>% 
    # select(Study_number, Biome, Ecosystem_type, Leaf_habit, Latitude, Longitude, Ra_annual, ER, GPP, NPP) %>% 
    select(Ecosystem_type, Rs_annual, Ra_annual, Rh_annual, ER, GPP, NPP) %>% 
    filter(!is.na(GPP) | !is.na(NPP) | !is.na(ER) ) -> 
    sub_srdb
  sub_srdb$ER <- ifelse(is.na(sub_srdb$ER), sub_srdb$GPP, sub_srdb$ER) # IF ER not available, use GPP (assume NEP very small).
  sub_srdb$Rh_annual <- ifelse(is.na(sub_srdb$Rh_annual), sub_srdb$NPP, sub_srdb$Rh_annual) # IF Rh_annual not available, use NPP (assume NEP very small)
  sub_srdb$Ra_annual <- ifelse(is.na(sub_srdb$Ra_annual), sub_srdb$Rs_annual - sub_srdb$Rh_annual, sub_srdb$Ra_annual) # When Ra_annual = Rs_annual - Rh_annual
  # If RE available but Ra_annual not availabel, we could also estimate Ra_annual based on the Bond_Lamberty (2004) model
  # We can get 22 more data, but don't know whether should keep those points (it does not change the results too much)
  # sub_srdb$Ra_annual <- ifelse(is.na(sub_srdb$Ra_annual), (-7.97+0.93*sqrt(sub_srdb$Rs_annual))^2, sub_srdb$Ra_annual )
  
  # In srdb Ra_annual = Rroot, so Froot = Rroot / Ra (fraction of Rroot to autotrophic respiration)
  sub_srdb$Froot <- sub_srdb$Ra_annual/(sub_srdb$ER - sub_srdb$Rh_annual)
  sub_srdb <- subset(sub_srdb, !is.na(Froot) & Froot < 0.9 & Froot > 0.1)
  
  # mean(sub_srdb$Froot)
  
  sub_srdb$Fshoot <- 1 - sub_srdb$Froot
  sub_srdb %>% select(Fshoot, Froot, Ecosystem_type) -> sub_srdb
  colnames(sub_srdb) <- colnames(comb_data)
  
  comb_data <- rbind(comb_data, sub_srdb)
  comb_data[comb_data$Ecosystem == 'Cropland', ]$Ecosystem <- "Agriculture"
  return (comb_data)
}



# plot global Rs
plot_Rs <- function (sdata, sdata2, sdata3) {
  sub_Rs <- subset(sdata, !is.na(sdata$Rs))
  Rs_obs <- nrow(sub_Rs)
  sub_Rs$label <- paste0("n=", Rs_obs)
  p_Rs <- ggplot(sub_Rs, aes(x = label, y=Rs)) + geom_violin() +
    geom_jitter(shape=16, position=position_jitter(0.2), col = 'gray') +
    geom_boxplot(width=.1) +
    # geom_point(aes(y=71.67), col="red") + geom_point(aes(y=76.68), col="red") +
    xlab("Global Rs") + ylab(expression(Rs~"("~Pg~yr^{-1}~")")) +
    theme(axis.title.x=element_blank())
  
  sub_Rh <- subset(sdata, !is.na(sdata$Rh))
  sub_Rh$label <- "n=2"
  p_Rh <- ggplot(sub_Rh, aes(x = label, y=Rh)) + geom_boxplot(width=.1) +
    geom_point(colour = "gray", size = 2) +
    geom_point(aes(y=46.47), col="red") +
    xlab("Global Rh") + ylab(expression(Rh~"("~Pg~yr^{-1}~")")) +
    theme(axis.title.x=element_blank())
    
  
  sub_Ra <- subset(sdata, !is.na(sdata$Ra))
  sub_Ra$label <- "n=1"
  p_Ra <- ggplot(sub_Ra, aes(x = label, y=Ra)) + geom_boxplot(width=.1) +
    geom_point(colour = "gray", size = 2) +
    geom_point(aes(y=25.20), col="red") + geom_point(aes(y=30.21), col="red") +
    xlab("Global Ra") + ylab(expression(Ra~"("~Pg~yr^{-1}~")")) +
    theme(axis.title.x=element_blank())
  
  # RC and RH distribution
  # RC
  sdata2 %>% 
    select(Study_number, Biome, Ecosystem_type, Leaf_habit, Latitude, Longitude, RC_annual, Ra_annual, Rh_annual) %>% 
    filter(RC_annual > 0 & RC_annual < 1) -> 
    sub_srdb 
  
  sub_srdb %>% group_by(Ecosystem_type) %>% summarise(count = n()) %>% print ()
  
  RC_obs <- nrow(sub_srdb)
  # sub_srdb %>% filter(Ra_Rh_Ratio < 5) -> sub_srdb
  sub_srdb$Global <- paste0("n=", RC_obs)
  
  RC_plot <- ggplot(sub_srdb, aes(x = Global, y=RC_annual)) + geom_violin() +
    geom_jitter(shape=16, position=position_jitter(0.2), col = 'gray') +
    geom_boxplot(width=.1) +
    # stat_summary(fun.y=median, geom="point", size=2, color="red") +
    ylab("RC") +
    theme(axis.title.x=element_blank())
  
  # RC by biome
  RC_plot2 <- ggplot(sub_srdb, aes(x = Ecosystem_type, y=RC_annual)) + geom_violin() +
    geom_jitter(shape=16, position=position_jitter(0.2), col = 'gray') +
    geom_boxplot(width=.1) +
    scale_x_discrete(limits = c("Agriculture", "Desert", "Forest", "Grassland", "Savanna", "Shrubland", "Wetland"),
                     labels = c("Agriculture (n=40)", "Desert (n=3)", "Forest (n=480)", "Grassland (n=74)", "Savanna (n=3)", "Shrubland(n=11)", "Wetland (n=6)") ) +
    # stat_summary(fun.y=median, geom="point", size=2, color="red") +
    ylab("RC") +
    theme(axis.title.x=element_blank())
  
  
  print(mean(sub_srdb$RC_annual))
  se <- sd(sub_srdb$RC_annual)/sqrt(nrow(sub_srdb))
  paste0("RC se=", se) %>% print()
  CI <- round(qt(0.975,df=nrow(sub_srdb)-1)*se, 3)
  paste0("RC 95% CI=", CI) %>% print()
  paste0("RC obs=", nrow(sub_srdb)) %>% print()
  
  # plot NPP
  sdata3$Global <- "n=251"
  p_NPP <- ggplot(sdata3, aes(x = Global, y=NPP)) + geom_violin() +
    geom_jitter(shape=16, position=position_jitter(0.2), col = 'gray') +
    geom_boxplot(width=.1) +
    ylab(expression(NPP~"("~Pg~yr^{-1}~")")) +
    # stat_summary(fun.y=median, geom="point", size=2, color="red") +
    theme(axis.title.x = element_blank())
  print(median(sdata3$NPP))
  se <- sd(sdata3$NPP)/sqrt(nrow(sdata3))
  paste0("NPP se=", se) %>% print()
  CI <- round(qt(0.975,df=nrow(sdata3)-1)*se, 3)
  paste0("NPP 95% CI=", CI) %>% print()
  
  # plot Fire
  p_fire <- tibble(Fire = c(2, 3.5, 7.3, 4, 5.1, 2.02, 2.71, 3.02, 2.08)) %>% ggplot(., aes(x="Fire", y = Fire)) + geom_violin() +
    geom_jitter(shape=16, position=position_jitter(0.2), col = 'gray') +
    geom_boxplot(width=.1) +
    ylab(expression(Fire~burned~"("~Pg~yr^{-1}~")")) +
    # stat_summary(fun.y=median, geom="point", size=2, color="red") +
    theme(axis.title.x = element_blank())
    
  
  p_mul <- plot_grid(p_Rs, p_NPP, p_fire, nrow = 1
            , labels = c('( a )', '( b )', '( c )')
            , vjust = c(3), hjust = c(-2.8, -2.8, -2.5)) 
  
  plot_grid(p_mul, RC_plot2, nrow = 2
            , labels = c('', '( d )')
            , vjust = c(3), hjust = c(-2.8) )

}


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
#   var_Rs <- GlobalRs %>% filter(!is.na(Rs) ) %>% select(Rs)
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

