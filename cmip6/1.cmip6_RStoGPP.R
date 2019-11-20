## 1.cmip6_RStoGPP.R
## This script processes the CMIP6 files for the GlobalC project. 
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
PROJECT_DIR   <- '/pic/projects/GCAM/Dorheim/Dorheim/GlobalC' # Where the GlobalC project lives on pic 
INPUT_DIR     <- file.path(PROJECT_DIR,'cmip6', 'output')             # Where the csv files containing the information about the files to process live.
OUTPUT_DIR    <- file.path(INPUT_DIR) # Where to save the final csv file. 
INTERMED_DIR  <- file.path(OUTPUT_DIR, 'intermediate-ratio-csv'); dir.create(INTERMED_DIR) # Where to save all the intermediate ncdf files 
CMIP6_DIR     <- '/pic/projects/GCAM/CMIP6'

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
  filter(experiment %in% c('historical', 'esm-hist')) -> 
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
  na.omit() -> 
  to_process



# 2. Calculate the Weighted Mean Ratio  ------------------------------------------------------------

apply(to_process, 1, function(input){
  
  # Define the base name to use for the intermediate files that are saved during this process. 
  base_name <- paste(input[["model"]], input[["experiment"]], input[["ensemble"]], input[["grid"]],input[["time"]], sep  = '_')

  # Messages
  print('------------------------------------------------') 
  print(base_name)  

  # Import a nc to use as a refernce. 
  input[['raRoot']] %>% 
    nc_open -> 
    nc
  
  # Import all of the data from the netcdfs. 
  input[['raRoot']] %>% 
    nc_open %>% 
    ncvar_get('raRoot') -> 
    raRoot 
  
  input[['rhSoil']] %>% 
    nc_open %>% 
    ncvar_get('rhSoil') -> 
    rhSoil
  
  input[['gpp']] %>% 
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
  
  rs_total     <- raRoot + rhSoil
  ratio        <- rs_total / gpp
  ratio[zeros] <- 0 
  
  
  # Calculate the land area with in each grid cell, this will 
  # be used as weights for the global weighted mean. 
  land_area <- areacella * (sftlf / 100)
  
  # Calculate the weighted global average. 
  value <- apply(ratio, MARGIN = 3, weighted.mean, w = land_area) 
  
  # Format the output. 
  data.frame(value = value, 
             units = NA, 
             time = ncvar_get(nc, 'time')) %>%  
    mutate(era = 'cmip6', 
           variable = 'ratio', 
           model = input[["model"]], 
           experiment = input[["experiment"]],
           ensemble = input[["ensemble"]], 
           grid = input[["grid"]]) -> 
    rslt 
  
  write.csv(rslt, file = file.path(INTERMED_DIR, paste0(base_name, '.csv')), row.names = FALSE)
  
  rslt
  
}) %>%  
  bind_rows %>% 
  mutate(year = substr(time, 1, 4), 
         motnh = substr(time, 5, 6)) %>%
  write.csv(., file = file.path(OUTPUT_DIR, 'cmip6-ratio-RS2GPP-fromgridcell.csv'), row.names = FALSE)









