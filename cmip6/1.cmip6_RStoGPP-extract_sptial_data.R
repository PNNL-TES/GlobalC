## 1.cmip6_RStoGPP-extract_spatial_data.R
## Calculate annual ratio of RS to GPP for each grid cell. 
## It determines the gridcell soil respiration (RS) to gross primary production (GPP)
## ratio where soil respiration is equal to the the sum of root and soil respiraiton 
## (raRoot and rhSoil). 

# 0. Set Up --------------------------------------------------------------------------------------
# Load the required R packages. 
library(dplyr)
library(tidyr)
library(tibble)
library(ncdf4)

# Define directories.  
PROJECT_DIR    <- '/pic/projects/GCAM/Dorheim/Dorheim/GlobalC' # Where the GlobalC project lives on pic 
INPUT_DIR      <- file.path(PROJECT_DIR,'cmip6', 'output')             # Where the csv files containing the information about the files to process live.
OUTPUT_DIR     <- file.path(INPUT_DIR) # Where to save the final csv file. 
INTERMED_DIR   <- file.path(OUTPUT_DIR, 'intermediate-girdded-csv'); dir.create(INTERMED_DIR) # Where to save all the intermediate ncdf files 
INTERMEDNC_DIR <- file.path(OUTPUT_DIR, 'intermediate-ncs'); dir.create(INTERMEDNC_DIR) 
CMIP6_DIR      <- '/pic/projects/GCAM/CMIP6'

# Define the path to the CDO
CDO_EXE <- "/share/apps/netcdf/4.3.2/gcc/4.4.7/bin/cdo" # Define the cdo directory.

# Option to delete all of the intermediate nc and csv files. 
CELAN_UP <- FALSE

# 1. Find Files to Process  ------------------------------------------------------------

# Import the cmip6 archive index. 
read.csv(file.path(CMIP6_DIR, 'cmip6_archive_index.csv'), stringsAsFactors = FALSE) %>% 
  # We are looking for the monthly GPP, raRoot, rhLitter, rh, and rhSoil for the 
  # historical experiments.
  filter(variable %in% c('gpp', 'raRoot', 'rhSoil')) %>% 
  filter(experiment %in% c('historical')) -> 
  data_files

# Figure out what files have sufficent data to process (need rhSoil, raRoot, and gpp files). 
data_files %>% 
  select(-domain) %>%  
  spread(variable, file) %>%
  na.omit -> 
  data_to_process

read.csv(file.path(CMIP6_DIR, 'cmip6_archive_index.csv'), stringsAsFactors = FALSE) %>% 
  filter(variable %in% c('areacella', 'sftlf')) %>% 
  select(file, variable, model, experiment, ensemble, grid) %>% 
  spread(variable, file) %>% 
  na.omit -> 
  fx_files 

data_to_process %>% 
  left_join(fx_files) %>% 
  na.omit() %>% 
  # Make sure that we are processing the netcdf files that actually have 
  # data for the years that we want to extract. 
  filter(grepl('201412', time)) -> 
  to_process

# The years of data to extract, if we end up being interested in ratios that span  mulitple 
# netcdf files the code will have to be adjusted so that the netcdfs for the same model are 
# concatenated together. 
years <- 2010
  
# 2. Calculate the Rs to GPP Ratio  ------------------------------------------------------------

apply(to_process, 1, function(input){
  
  # Define the base name to use for the intermediate files that are saved during this process. 
  base_name <- paste(input[["model"]], input[["experiment"]], input[["ensemble"]], input[["grid"]], sep  = '_')

  # Messages
  print('------------------------------------------------') 
  print(base_name)  
  
  # Create the intermediate netcdf files. 
  annual_gppNC    <- file.path(INTERMEDNC_DIR, paste0(base_name, 'annual_gpp.nc'))
  annual_raRootNC <- file.path(INTERMEDNC_DIR, paste0(base_name, 'annual_raRoot.nc'))
  annual_rhSoilNC <- file.path(INTERMEDNC_DIR, paste0(base_name, 'annual_rhSoil.nc'))
  time_nc         <- file.path(INTERMED_DIR, paste0(base_name, '_time.nc'))
  
  # Take the annual mean from monthly data. 
  system2(CDO_EXE, args = c('yearmonmean', input[['gpp']], annual_gppNC), stdout = TRUE, stderr = TRUE )
  system2(CDO_EXE, args = c('yearmonmean', input[['raRoot']], annual_raRootNC), stdout = TRUE, stderr = TRUE )
  system2(CDO_EXE, args = c('yearmonmean', input[['rhSoil']], annual_rhSoilNC), stdout = TRUE, stderr = TRUE )
  
  # Convert from realtive time to absolute time. 
  system2(CDO_EXE, args = c('-a', '-copy', annual_raRootNC, time_nc), stdout = TRUE, stderr = TRUE )
  
  # Import a nc to use as a refernce. 
  time_nc %>% 
    nc_open -> 
    nc
  
  # Import all of the annual netcdf files. 
  annual_raRootNC %>% 
    nc_open %>% 
    ncvar_get('raRoot') -> 
    raRoot 
  
  annual_rhSoilNC %>% 
    nc_open %>% 
    ncvar_get('rhSoil') -> 
    rhSoil
  
  annual_gppNC %>% 
    nc_open %>% 
    ncvar_get('gpp') -> 
    gpp
  
  input[['sftlf']] %>%  
    nc_open %>% 
    ncvar_get('sftlf') -> 
    sftlf
  
  input[['areacella']] %>% 
    nc_open %>% 
    ncvar_get('areacella') -> 
    areacella
  
  # Define what we expect 0 to be. 
  practically_zero <- 3e-11
  
  # Identify all of the 0s in the netcdfs, these 
  # index values will be used to replace value from the 
  # ratio with 0s.
  raRoot_zero <- which(abs(raRoot) <= practically_zero)
  rhSoil_zero <- which(abs(rhSoil) <= practically_zero)
  gpp_zero    <- which(abs(gpp) <= practically_zero)
  zeros       <- c(raRoot_zero, rhSoil_zero, gpp_zero)
  
  # Calcualte the total below ground respiration and use 
  # this value to get the Rs to GPP ratio. Finally 
  # replace values with 0s, (this will remove the values that 
  # were Inf or large because of the machine percision).
  rs_total     <- raRoot + rhSoil
  ratio        <- rs_total / gpp
  ratio[zeros] <- 0 
  
  # Convert the time column into years. 
  year      <-   as.numeric(substr(ncvar_get(nc, 'time'), 1, 4))
  
  # For each of the years to save extract the data from the ratio 
  # 3d array and format it as a data frame. 
  lapply(years, function(y){
    
    # The data is stored as lon, lat, year
    to_format <- ratio[ , ,which(year == y)]
    lat <- ncvar_get(nc, 'lat')
    lon <- ncvar_get(nc, 'lon')
    
    
    expand.grid(longitude = lon, 
                latitude = lat) %>%  
      mutate(value = as.vector(ratio[,,which(year == 2010)])) %>%  
      mutate(era = 'cmip6', 
             variable = 'ratio', 
             model = input[["model"]], 
             experiment = input[["experiment"]],
             ensemble = input[["ensemble"]], 
             grid = input[["grid"]]) %>%  
      mutate(year = y)  -> 
      rslt 
    
    out_file <- file.path(OUTPUT_DIR, paste0(base_name, '-', y, '-gridded-ratios.csv'))
    write.csv(rslt, file = out_file, row.names = FALSE)
  out_file
  }) %>% 
    unlist 
}) %>% 
  unlist -> 
  files 










