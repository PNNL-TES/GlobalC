
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
  
  sdata %>% 
    select(Study_number, Biome, Leaf_habit, Latitude, Longitude, RC_annual) %>% 
    filter(RC_annual > 0 & RC_annual < 1) -> 
    sub_srdb 
  # sub_srdb$RC_RH_Ratio <- sub_srdb$RC_annual/(1-sub_srdb$RC_annual)
  
  x <- rep(-170, 2)
  size <- c(2, 1.5)
  y <- c(-25, -40)
  
  cc_legend <- cbind (x, y, size)
  cc_legend <- as.data.frame(cc_legend)
  
  
  siteInfor <- summarySE (data=sub_srdb, measurevar='RC_annual', groupvars=c("Latitude","Longitude"))
  siteInfor <- siteInfor[, c(1:3)]
  siteInfor <- siteInfor[which(!is.na(siteInfor$Latitude)),]
    
  # siteInfor$N <- siteInfor$N * 0.25
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
    # FlFsFr sites
    geom_point(data = sdata2, aes(x=sdata2$Longitude, y=sdata2$Latitude), color = "red"
               ,shape = 3, size = 1.5, stroke = 1.5, alpha = 1, fill = "red") +
    # legend
    geom_point(data = cc_legend ,aes(x=cc_legend$x, y=cc_legend$y), shape = c(18, 3), stroke = c(1, 1.5) 
               , color = c("black", "red"), size = cc_legend$size, alpha = c(1, 1)) +
    annotate("text", x = -150, y = -10, label = "Legend", size = 4, hjust = 0) +
    annotate("text", x = -150, y = -25, label = "Rab/Rh sites", size = 4, hjust = 0) +
    annotate("text", x = -150, y = -40, label = "FlFsFr sites", size = 4, hjust = 0) +
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
plot_NPP_Rab <- function (sdata, sdata2) {
  sdata$Global <- "n=251"
  p_NPP <- ggplot(sdata, aes(x = Global, y=NPP)) + geom_violin() +
    geom_jitter(shape=16, position=position_jitter(0.2), col = 'gray') +
    geom_boxplot(width=.1) +
    ylab(expression(NPP~"("~Pg~yr^{-1}~")")) +
    # stat_summary(fun.y=median, geom="point", size=2, color="red") +
    theme(axis.title.x = element_blank())
  print(median(sdata$NPP))
  se <- sd(sdata$NPP)/sqrt(nrow(sdata))
  paste0("NPP se=", se) %>% print()
  CI <- round(qt(0.975,df=nrow(sdata)-1)*se, 3)
  paste0("NPP 95% CI=", CI) %>% print()
  
  # Rab/Rh
  sdata2 %>% 
    select(Study_number, Biome, Leaf_habit, Latitude, Longitude, Ra_annual, Rh_annual) %>% 
    filter(Ra_annual > 0 & Rh_annual > 0) -> 
    sub_srdb 
  sub_srdb$Rab_Rh <- sub_srdb$Ra_annual/sub_srdb$Rh_annual
  # sub_srdb %>% filter(Ra_Rh_Ratio < 5) -> sub_srdb
  sub_srdb$Global <- "n=602"
  
  Rab_plot <- ggplot(sub_srdb, aes(x = Global, y=Rab_Rh)) + geom_violin() +
    geom_jitter(shape=16, position=position_jitter(0.2), col = 'gray') +
    geom_boxplot(width=.1) +
    # stat_summary(fun.y=median, geom="point", size=2, color="red") +
    ylab("Rab / Rh") +
    theme(axis.title.x=element_blank())
  
  print(median(sub_srdb$Rab_Rh))
  se <- sd(sub_srdb$Rab_Rh)/sqrt(nrow(sub_srdb))
  paste0("Rab se=", se) %>% print()
  CI <- round(qt(0.975,df=nrow(sub_srdb)-1)*se, 3)
  paste0("Rab 95% CI=", CI) %>% print()
  paste0("Rab obs=", nrow(sub_srdb)) %>% print()
  
  # print figure
  plot_grid(p_NPP, Rab_plot, labels = c('( a )', '( b )')
            , vjust = c(3.5,3.5), hjust = c(-2.5, -2.25) )
}


# plot Ff
#install.packages("reshapes") # Old package to do the same thing as gather
# library(reshape2)
# FsFfFrReshape <- melt(FsFfFr,id.vars='IGBP', measure.vars=c('Fl','Fs','Fr'))
plot_Ff <- function (sdata) {
  sdata %>%
    select(Fl, Fs, Fr) %>% 
    gather(Fraction, Percentage, Fl:Fr, na.rm = TRUE) ->
    comb_data
  
  Fl_plot <- ggplot(comb_data, aes(x = Fraction, y=Percentage)) + geom_violin() +
    geom_jitter(shape=16, position=position_jitter(0.2), col = 'gray') +
    geom_boxplot(width=.1) +
    ylab("Percentage (%)") +
    # stat_summary(fun.y=median, geom="point", size=2, color="red") +
    scale_x_discrete(labels=c("Fl" = "Fl (n=20)", "Fr" = "Fr (n=27)", "Fs" = "Fs (n=18)")) +
    theme(axis.title.x=element_blank())
  
  print(median(sdata$Fr, na.rm = T) / (median(sdata$Fl, na.rm = T)+median(sdata$Fs, na.rm = T)+median(sdata$Fr, na.rm = T))*100)
  se <- sd(sdata$Fr, na.rm = T)/nrow(subset(sdata, !is.na(sdata$Fr)))
  paste0("se=", se) %>% print()
  CI <- round(qt(0.975,df=nrow(subset(sdata, !is.na(sdata$Fr)))-1)*se, 3)
  paste0("Fr 95% CI=", CI," ,obs=", nrow(subset(sdata, !is.na(sdata$Fr)))) %>% print()
  
  print(median(sdata$Fs, na.rm = T) / (median(sdata$Fl, na.rm = T)+median(sdata$Fs, na.rm = T)+median(sdata$Fr, na.rm = T))*100)
  se <- sd(sdata$Fs, na.rm = T)/nrow(subset(sdata, !is.na(sdata$Fs)))
  paste0("se=", se) %>% print()
  CI <- round(qt(0.975,df=nrow(subset(sdata, !is.na(sdata$Fs)))-1)*se, 3)
  paste0("Fs 95% CI=", CI, " ,obs=", nrow(subset(sdata, !is.na(sdata$Fs)))) %>% print()
  
  print(median(sdata$Fl, na.rm = T) / (median(sdata$Fl, na.rm = T)+median(sdata$Fs, na.rm = T)+median(sdata$Fr, na.rm = T))*100)
  se <- sd(sdata$Fl, na.rm = T)/nrow(subset(sdata, !is.na(sdata$Fl)))
  paste0("se=", se) %>% print()
  CI <- round(qt(0.975,df=nrow(subset(sdata, !is.na(sdata$Fl)))-1)*se,3) 
  paste0("Fl 95% CI=", CI ," ,obs=", nrow(subset(sdata, !is.na(sdata$Fl)))) %>% print()
  
  print(Fl_plot)
}

# plot global Rs

plot_Rs <- function (sdata) {
  sub_Rs <- subset(sdata, !is.na(sdata$Rs))
  sub_Rs$label <- "n=17"
  p_Rs <- ggplot(sub_Rs, aes(x = label, y=Rs)) + geom_violin() +
    geom_jitter(shape=16, position=position_jitter(0.2), col = 'gray') +
    geom_boxplot(width=.1) +
    geom_point(aes(y=71.67), col="red") + geom_point(aes(y=76.68), col="red") +
    xlab("Global Rs") + ylab(expression(Rs~"("~Pg~yr^{-1}~")")) +
    theme(axis.title.x=element_blank())
  
  sub_Rh <- subset(sdata, !is.na(sdata$Rh))
  sub_Rh$label <- "n=3"
  p_Rh <- ggplot(sub_Rh, aes(x = label, y=Rh)) + geom_boxplot(width=.1) +
    geom_point(colour = "gray", size = 2) +
    geom_point(aes(y=46.47), col="red") +
    xlab("Global Rh") + ylab(expression(Rh~"("~Pg~yr^{-1}~")")) +
    theme(axis.title.x=element_blank())
    
  
  sub_Ra <- subset(sdata, !is.na(sdata$Ra))
  sub_Ra$label <- "n=3"
  p_Ra <- ggplot(sub_Ra, aes(x = label, y=Ra)) + geom_boxplot(width=.1) +
    geom_point(colour = "gray", size = 2) +
    geom_point(aes(y=25.20), col="red") + geom_point(aes(y=30.21), col="red") +
    xlab("Global Ra") + ylab(expression(Ra~"("~Pg~yr^{-1}~")")) +
    theme(axis.title.x=element_blank())
    
  
  plot_grid(p_Rs, p_Rh, p_Ra, ncol = 3
            , labels = c('( a )', '( b )', '( c )')
            , vjust = c(3, 3, 3), hjust = c(-2, -2.25, -2.25)) 

}

 



