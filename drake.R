# load packages
library(tibble)
library(lubridate)
library(kableExtra)
library(raster)
library(ncdf4)
library(tidyr)
library(dplyr)
library(magrittr)
library(data.table)

## functions
# create a function to get precipitation and temperature for mgrsd
get_climate <- function (sdata) {
  drake::readd(precip_global) -> wordprec
  drake::readd(tm_global) -> wordtm
  drake::readd(land_evi) -> wordevi
  
  for(i in 1:nrow(sdata)){
  # for(i in 1:100){
    target_lat <- sdata$Latitude[i]
    target_lon <- sdata$Longitude[i]
    target_mon <- sdata$Meas_Month2[i]
    
    if (is.na(target_lat) | is.na(target_lon) | is.na(target_mon)) {
      sdata$Pmonth[i] = NA
      sdata$Tmonth[i] = NA
      sdata$EVI[i] = NA
    }
    else{
      ilon <- wordprec$x[which.min(abs(wordprec$x - target_lon))]
      ilat <- wordprec$y[which.min(abs(wordprec$y - target_lat))]
      
      iperc <- wordprec %>% filter(x == ilon & y == ilat)
      itm <- wordtm %>% filter(x == ilon & y == ilat)
      ievi <- wordevi %>% filter(Longitude == ilon & Latitude == ilat)
      
      if (nrow(iperc) == 0 | nrow(itm) == 0 | nrow(ievi) == 0) {
        sdata$Pmonth[i] = NA
        sdata$Tmonth[i] = NA
        sdata$EVI[i] = NA
      }
      
      else {
        if (target_mon == 1) {
          sdata$Pmonth[i] = iperc$prec1
          sdata$Tmonth[i] = itm$tmean1/10
          sdata$EVI[i] = ievi$EVI_1 }
        else if (target_mon == 2) {
          sdata$Pmonth[i] = iperc$prec2
          sdata$Tmonth[i] = itm$tmean2/10
          sdata$EVI[i] = ievi$EVI_2 }
        else if (target_mon == 3) {
          sdata$Pmonth[i] = iperc$prec3
          sdata$Tmonth[i] = itm$tmean3/10
          sdata$EVI[i] = ievi$EVI_3 }
        else if (target_mon == 4) {
          sdata$Pmonth[i] = iperc$prec4
          sdata$Tmonth[i] = itm$tmean4/10
          sdata$EVI[i] = ievi$EVI_4 }
        else if (target_mon == 5) {
          sdata$Pmonth[i] = iperc$prec5
          sdata$Tmonth[i] = itm$tmean5/10
          sdata$EVI[i] = ievi$EVI_5 }
        else if (target_mon == 6) {
          sdata$Pmonth[i] = iperc$prec6
          sdata$Tmonth[i] = itm$tmean6/10
          sdata$EVI[i] = ievi$EVI_6 }
        else if (target_mon == 7) {
          sdata$Pmonth[i] = iperc$prec7
          sdata$Tmonth[i] = itm$tmean7/10
          sdata$EVI[i] = ievi$EVI_7 }
        else if (target_mon == 8) {
          sdata$Pmonth[i] = iperc$prec8
          sdata$Tmonth[i] = itm$tmean8/10
          sdata$EVI[i] = ievi$EVI_8 }
        else if (target_mon == 9) {
          sdata$Pmonth[i] = iperc$prec9
          sdata$Tmonth[i] = itm$tmean9/10
          sdata$EVI[i] = ievi$EVI_9 }
        else if (target_mon == 10) {
          sdata$Pmonth[i] = iperc$prec10
          sdata$Tmonth[i] = itm$tmean10/10
          sdata$EVI[i] = ievi$EVI_10 }
        else if (target_mon == 11) {
          sdata$Pmonth[i] = iperc$prec11
          sdata$Tmonth[i] = itm$tmean11/10
          sdata$EVI[i] = ievi$EVI_11 }
        else {
          sdata$Pmonth[i] = iperc$prec12
          sdata$Tmonth[i] = itm$tmean12/10
          sdata$EVI[i] = ievi$EVI_12 }
      }
    }
    print(paste0(target_lon, "-->", ilon, "***", target_lat, "-->", ilat, "***", i))
  }
  return(sdata)
}


#*****************************************************************************************************************
# plan 
#*****************************************************************************************************************
plan <- drake::drake_plan(
  mgrsd = read.csv(here::here("data", "MGRhD_SRDBV5_20200819.csv")), 
  
  # Extract above/belowground biomass
  BM_aboveground = raster(here::here("data/extdata", "aboveground_biomass_carbon_2010.tif")) %>%
    raster::extract(mgrsd %>% dplyr::select(Longitude, Latitude)),
  BM_belowground = raster(here::here("data/extdata", "belowground_biomass_carbon_2010.tif")) %>%
    raster::extract(mgrsd %>% dplyr::select(Longitude, Latitude)),
  
  BM_aboveground_global = raster(here::here("data/extdata", "aboveground_biomass_carbon_2010.tif")) %>%
    raster::extract(landgrids %>% dplyr::select(Longitude, Latitude)),
  BM_belowground_global = raster(here::here("data/extdata", "belowground_biomass_carbon_2010.tif")) %>%
    raster::extract(landgrids %>% dplyr::select(Longitude, Latitude)),
  
  ## Extracts N deposition data for each mgrsd point
  #  https://www.isimip.org/gettingstarted/availability-input-data-isimip2b/
  N_dep_half_deg = raster(here::here("data/extdata", "global_mean_Nitrogen_depostion_1980_2017_half_degree.tif")) %>% 
    raster::extract(mgrsd %>% dplyr::select(Longitude, Latitude)),
  
  N_dep_global = raster(here::here("data/extdata", "global_mean_Nitrogen_depostion_1980_2017_half_degree.tif")) %>% 
    raster::extract(landgrids %>% dplyr::select(Longitude, Latitude)),
  
  # extract BD, clay percentage, and SOC from SoilGrids
  bd_soilgrids = raster("data/extdata/BLDFIE_M_sl2_1km_ll.tif") %>%
    raster::extract(mgrsd %>% dplyr::select(Longitude, Latitude)),
  clay_soilgrids = raster("data/extdata/CLYPPT_M_sl2_1km_ll.tif") %>%
    raster::extract(mgrsd %>% dplyr::select(Longitude, Latitude)),
  soc_soilgrids = raster("data/extdata/OCSTHA_M_sd2_1km_ll.tif") %>%
    raster::extract(mgrsd %>% dplyr::select(Longitude, Latitude)),
  
  bd_soilgrids_global = raster("data/extdata/BLDFIE_M_sl2_1km_ll.tif") %>%
    raster::extract(landgrids %>% dplyr::select(Longitude, Latitude)),
  clay_soilgrids_global = raster("data/extdata/CLYPPT_M_sl2_1km_ll.tif") %>%
    raster::extract(landgrids %>% dplyr::select(Longitude, Latitude)),
  soc_soilgrids_global = raster("data/extdata/OCSTHA_M_sd2_1km_ll.tif") %>%
    raster::extract(landgrids %>% dplyr::select(Longitude, Latitude)),
  
  # get EVI 
  landgrids = read.table(here::here("data/extdata", "landgrids.dat"), sep = ",") %>% 
    rename(Longitude = V2, Latitude = V1),
  EVI_1 = raster("data/EVIs/EVI_mean_month_1.tif") %>%
    raster::extract(landgrids %>% dplyr::select(Longitude, Latitude)),
  EVI_2 = raster("data/EVIs/EVI_mean_month_2.tif") %>%
    raster::extract(landgrids %>% dplyr::select(Longitude, Latitude)),
  EVI_3 = raster("data/EVIs/EVI_mean_month_3.tif") %>%
    raster::extract(landgrids %>% dplyr::select(Longitude, Latitude)),
  EVI_4 = raster("data/EVIs/EVI_mean_month_4.tif") %>%
    raster::extract(landgrids %>% dplyr::select(Longitude, Latitude)),
  EVI_5 = raster("data/EVIs/EVI_mean_month_5.tif") %>%
    raster::extract(landgrids %>% dplyr::select(Longitude, Latitude)),
  EVI_6 = raster("data/EVIs/EVI_mean_month_6.tif") %>%
    raster::extract(landgrids %>% dplyr::select(Longitude, Latitude)),
  EVI_7 = raster("data/EVIs/EVI_mean_month_7.tif") %>%
    raster::extract(landgrids %>% dplyr::select(Longitude, Latitude)),
  EVI_8 = raster("data/EVIs/EVI_mean_month_8.tif") %>%
    raster::extract(landgrids %>% dplyr::select(Longitude, Latitude)),
  EVI_9 = raster("data/EVIs/EVI_mean_month_9.tif") %>%
    raster::extract(landgrids %>% dplyr::select(Longitude, Latitude)),
  EVI_10 = raster("data/EVIs/EVI_mean_month_10.tif") %>%
    raster::extract(landgrids %>% dplyr::select(Longitude, Latitude)),
  EVI_11 = raster("data/EVIs/EVI_mean_month_11.tif") %>%
    raster::extract(landgrids %>% dplyr::select(Longitude, Latitude)),
  EVI_12 = raster("data/EVIs/EVI_mean_month_12.tif") %>%
    raster::extract(landgrids %>% dplyr::select(Longitude, Latitude)),
  
  land_evi = cbind(landgrids, EVI_1/10000, EVI_2/10000, EVI_3/10000, EVI_4/10000, EVI_5/10000, EVI_6/10000,
                   EVI_7/10000, EVI_8/10000, EVI_9/10000, EVI_10/10000, EVI_11/10000, EVI_12/10000) %>% 
    rename("EVI_1" = "EVI_1/10000", "EVI_2" = "EVI_2/10000", "EVI_3" = "EVI_3/10000", "EVI_4" = "EVI_4/10000", "EVI_5" = "EVI_5/10000", "EVI_6" = "EVI_6/10000",
           "EVI_7" = "EVI_7/10000", "EVI_8" = "EVI_8/10000", "EVI_9" = "EVI_9/10000", "EVI_10" = "EVI_10/10000", "EVI_11" = "EVI_11/10000", "EVI_12" = "EVI_12/10000"),
  
  # temperature and precipitation data
  # Temp data is stored in degC * 10, so we need to divide to get back to degC
  precip = getData("worldclim", path = here::here(), var = "prec", res = 10, download = !file.exists("/Users/Documents/PNNL/SRPartitioning/wc10/prec1.hdr")),
  tmean = getData("worldclim", path = here::here(), var = "tmean", res = 10, download = !file.exists("/Users/Documents/PNNL/SRPartitioning/wc10/tmean1.hdr")),
  
  land_prec = precip %>% raster::extract(landgrids %>% dplyr::select(Longitude, Latitude)),
  land_tm = tmean %>% raster::extract(landgrids %>% dplyr::select(Longitude, Latitude)),
  
  # Calculate annual sum for precipitation...
  precip_global = raster::as.data.frame(precip, xy = TRUE) %>% drop_na(),
  
  # Calculate annual sum for temperature...
  tm_global = raster::as.data.frame(tmean, xy = TRUE) %>% drop_na(), 
  
  # mgrsd with all factors extracted
  mgrsd_all = cbind(mgrsd, BM_aboveground, BM_belowground, N_dep_half_deg, bd_soilgrids, clay_soilgrids, soc_soilgrids) %>% 
    # unit from kg/m3 to g/cm3
    mutate(bd_soilgrids = bd_soilgrids / 1000),
  
  # get monthly world climate and evi for mgrsd
  mgrsd_climate = get_climate(mgrsd_all),
  
  # prepare global 0.1*0.1 resolution data for global RC_prediction map
  global_land = cbind(land_evi, land_prec, land_tm, BM_aboveground_global, BM_belowground_global, N_dep_global,
                      bd_soilgrids_global, clay_soilgrids_global, soc_soilgrids_global) %>% 
    mutate(tmean1 = tmean1/10, tmean2 = tmean2/10, tmean3 = tmean3/10, tmean4 = tmean4/10, tmean5 = tmean5/10, tmean6 = tmean6/10,
           tmean7 = tmean7/10, tmean8 = tmean8/10, tmean9 = tmean9/10, tmean10 = tmean10/10, tmean11 = tmean11/10, tmean12 = tmean12/10)
)

drake::make(plan)

