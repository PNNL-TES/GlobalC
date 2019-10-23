## 1.process_cmip5.R
## This script processes the CMIP5 files for the GlobalC project. 
## Running time depends on the resolution of the regions, the smaller the 
## regions are the longer the script will take to run. However, if you don't 
## clean up the intermediate files they can be used to prevent rerunning 
## portions of the script if the job terminates on pic early.

# 0. Set Up --------------------------------------------------------------------------------------
# Load the required R packages. 
library(dplyr)
library(tidyr)
library(tibble)
library(ncdf4)

# Define directories.  
PROJECT_DIR   <- '/pic/projects/GCAM/Dorheim/Dorheim/GlobalC' # Where the GlobalC project lives on pic 
INPUT_DIR     <- file.path(PROJECT_DIR,'cmip5', 'output')             # Where the csv files containing the information about the files to process live.
OUTPUT_DIR    <- file.path(INPUT_DIR) # Where to save the final csv file. 
INTERMED_DIR  <- file.path(OUTPUT_DIR, 'intermediate-ncs'); dir.create(INTERMED_DIR) # Where to save all the intermediate ncdf files 

# Define the path to the CDO
CDO_EXE <- "/share/apps/netcdf/4.3.2/gcc/4.4.7/bin/cdo" # Define the cdo directory.

# Option to delete all of the intermediate nc and csv files. 
CELAN_UP <- FALSE

# 1. Import Infomration To Process ------------------------------------------------------------
fx_files   <- read.csv(list.files(INPUT_DIR, pattern = '0.fx_files_CMIP5.csv', full.names = TRUE), stringsAsFactors = FALSE)
data_files <- read.csv(list.files(INPUT_DIR, pattern = '0.to_process_CMIP5.csv', full.names = TRUE), stringsAsFactors = FALSE)

# 2. Calculate the Land Area Weight .ncs -------------------------------------------------------
# Ensure that we are only processing the land area weights for the model, experiment, and ensemble
# members we have access to.
fx_files %>%  
  left_join(data_files %>%  
              select(model, experiment, ensemble) %>% 
              distinct, by = c("model", "experiment", "ensemble")) -> 
  meta_files_to_process

# Use cdo to create a nc of the land area weight, multiply the cell area by the fraction of the cell that is land.
# In the next section these files are going to be used to weight the averages and calculate total land per lattitude band.
apply(meta_files_to_process, MARGIN = 1, FUN = function(input){
  
  # Create the path to the output nc file
  out1_nc <- file.path(INTERMED_DIR, paste0(basename(input[['sftlf']]), '_LandArea1.nc'))
  out_nc <- file.path(INTERMED_DIR, paste0(basename(input[['sftlf']]), '_LandArea.nc'))
  if(file.exists(out_nc)){file.remove(out_nc)} # will probably want to remove this line
  assertthat::assert_that(!file.exists(out_nc), msg = 'land area weight exists')
  
  # Execute the cdo to create a nc of the land area weight, multiply the cell area by the fraction of the cell that is land.
  system2(CDO_EXE, args = c("-divc,100", input[['sftlf']], out1_nc), stdout = TRUE, stderr = TRUE)
  system2(CDO_EXE, args = c("-mul", input[['areacella']], out1_nc, out_nc), stdout = TRUE, stderr = TRUE)
  
  # Add the LandArea netcdf to the input information
  as_tibble(input) %>% 
    mutate(LandArea = out_nc)
  
  input[['LandArea']] <- out_nc
  input
  
}) %>% 
  # Format into a nice data frame, that can easily be joined with the data to process tibble. 
  unlist %>%  
  matrix(ncol = ncol(meta_files_to_process) + 1, 
         byrow = TRUE, 
         dimnames = list(NULL, c(names(meta_files_to_process), 'LandArea'))) %>%  
  as.data.frame(stringsAsFactors = FALSE) -> 
  LandAreaNc_tib

# Pair each of the data files with the land area file, using the land area files to weight the averages prevents 
# values from the ocean from being incorporated into the averages. This also ensures that the global and 
# large regional averages will not be skewed by the high latitude regions.
data_files %>% 
  full_join(LandAreaNc_tib %>% select(model, experiment, LandArea), 
            by = c("model", "experiment")) %>% 
  # This will remove files that do not have corresponding land area netcdf (usually this is because of the ensemble members)
  # TODO figure out some way to handle this better. 
  na.omit ->
  data_LandArea

# 3A. Regional Average Functions  -------------------------------------------------------
# Generate a tibble of the region to process. 
# Args
#   input: the netcdf file to process, it will provide the min and max lat and lon coordinates. 
#   region: if set to NULL will define the region as global, otherwise it may be a list of the resolution 
#           to split up the lat and lon regions into. 
# Returns - a data frame of the region and min/max lat/lon bounds to select when using the cdo to prcoess 
# the regional average for the netcdf in the cdo_regional_means function. 
create_region_df <- function(input, region){
  
  if(is.data.frame(region)){
    
    # TODO add code here to ensure that the region data frame is 
    # formatted correctly
  } else {
    
    nc        <- nc_open(input[['files']])
    lat_range <- signif(range(ncvar_get(nc, 'lat')), digits = 2)
    lon_range <- signif(range(ncvar_get(nc, 'lon')), digits = 2)
  }
  
  
  if(is.list(region) && all(names(region) %in% c('lat', 'lon'))){
    
    if(is.null(region$lat)){
      lat_values <- sort(lat_range)
    } else if(is.numeric(region$lat)){
      lat_values <- seq(from = min(lat_range), to = max(lat_range), by = region$lat )
    }
    
    if(is.null(region$lon)){
      lon_values <- sort(lon_range)
    } else if(is.numeric(region$lon)){
      lon_values <- seq(from = min(lon_range), to = max(lon_range), by = region$lon)
    }
    
    # Generate the regaion_dataframe that will be used to process 
    # the data netcdf files. 
    full_join(tibble(min_lat = lat_values[1:length(lat_values)-1], 
                     max_lat = lat_values[2:length(lat_values)], 
                     join = 1),
              tibble(min_lon = lon_values[1:length(lon_values) -1], 
                     max_lon = lon_values[2:length(lon_values)], 
                     join = 1), by = 'join') %>%  
      select(-join) %>% 
      mutate(region = paste0(min_lat, '_', max_lat, '_', min_lon, '_', max_lon)) -> 
      region_df
    
  } else if(is.null(region)){
    
    # Assume that the region to process is global. 
    tibble(min_lat = min(lat_range), 
           max_lat = max(lat_range), 
           min_lon = min(lon_range), 
           max_lon = max(lon_range)) %>% 
      mutate(region = 'global') -> 
      region_df
    
  } else {
    
    # The function is not set up to handle this sort of 
    stop('Format of region is unrecognized')
  }
  
  region_df
  
}

# Calucalte the weighted regional average of a netcdf. 
# Args 
#   data_LandArea: a dataframe containing the netcdf to process, CMIP6 meta data, and the land area netcdf. 
#   region_input: default is set to NULL, which will process the global weighted average. Otherwise it can be set to 
#   a list of the size of the lat and lon boxes to process the weighted regional average of. 
# Returns
#   messages - prints out messages about the netcdf and region that is being processed, this can make the slurm files messy but is helpful for debugging.  
#   intermediate netcdfs - Because CMIP6 files cannot be processed with pipelined CDO, this function saves lots of intermediate netcdf files. They eventually need to be cleaned up. 
#   intermediate csv files - In order to help with the debugging process and prevent process from being lost of the pic script terminates early this function saves all of the processed regional data as csv files. 
cdo_regional_means <- function(data_LandArea, region_input = NULL){
  
  # For every data netcdf file listed in the data_LandArea input calculate 
  # the regional average. 
  apply(data_LandArea, 1, function(input, region = region_input){
    
    # Messages
    print('------------------------------------------------') 
    print(input[['files']])  
    
    # Define the base name to use for the intermediate files that are saved during this process. 
    base_name <- paste(input[['era']], input[["variable"]], input[["domain"]],  input[["model"]], 
                       input[["experiment"]], input[["ensemble"]], input[["time"]], sep  = '_')
    
    # Generate the data frame that contains the lat and lon 
    # boundaries for the cdo calls. 
    region_df <- create_region_df(input, region)
    
    # For each of the defined regions calculate the regional average.
    apply(region_df, 1, function(region_input){
      
      # Extract regional infromation. 
      min_lat = region_input[['min_lat']]
      max_lat = region_input[['max_lat']]
      min_lon = region_input[['min_lon']]
      max_lon = region_input[['max_lon']] 
      region = region_input[['region']]
      
      print(region)
      
      # Make the intermediate netcdfs
      region_weights       <- file.path(INTERMED_DIR, paste0(base_name, '_', region, '.nc'))
      region_data          <- file.path(INTERMED_DIR, paste0(base_name, '_', region, '.nc'))
      region_data_weighted <- file.path(INTERMED_DIR, paste0(base_name, '_Weighted_', region, '.nc'))
      region_mean          <- file.path(INTERMED_DIR, paste0(base_name, '_Mean_', region, '.nc'))
      final_data           <- file.path(INTERMED_DIR, paste0(base_name, '_', region, '_final.nc'))
      
      # Format the sellonlatbox cdo call, this is the cdo command that will select a sepecfic lat / lon box.
      # The gubsub call was added when it turned out that switching from a negative value to a positivie 
      # coordinate value added a space that caused the cdo pipeline to fail.
      lat_cdo <- gsub(pattern = ' ', replacement = '', x = paste0( "sellonlatbox", ",", min_lon,",", max_lon,",", min_lat,",", max_lat, sep ="," ))
      
      # Select the Land Area weights for the specific region and calculate the total land area in the region.
      system2(CDO_EXE, args = c(lat_cdo, input[['LandArea']], region_weights), stdout = TRUE, stderr = TRUE)
      
      # Save the total land area and the units. 
      nc <- nc_open(region_weights)
      region_LandArea       <- sum(ncvar_get(nc, 'areacella'))
      region_LandArea_units <- ncatt_get(nc, 'areacella')$units
      
      # Extract the data for a specfic region. 
      system2(CDO_EXE, args = c(lat_cdo, input[['files']], region_data), stdout = TRUE, stderr = TRUE)
      # Define the area grid for the data netcdf to process. 
      system2(CDO_EXE, args = c(paste0("setgridarea,", region_weights), region_data, region_data_weighted), stdout = TRUE, stderr = TRUE)
      # Calculate the weighted average for the specific reigon. 
      system2(CDO_EXE, args = c('fldmean', region_data_weighted, region_mean), stdout = TRUE, stderr = TRUE )
      # Convert from relative time to absolute time. 
      system2(CDO_EXE, args = c('-a', '-copy', region_mean, final_data), stdout = TRUE, stderr = TRUE )
      
      # Extract data and format output.
      nc <- nc_open(final_data)
      data.frame(value = ncvar_get(nc, input[['variable']]), 
                 units = ncatt_get(nc, input[['variable']])$units, 
                 time = ncvar_get(nc, 'time'), 
                 region = region) %>%  
        mutate(LandArea = region_LandArea, 
               LandArea_units = region_LandArea_units, 
               era = input[['era']], 
               variable = input[["variable"]], 
               domain = input[["domain"]], 
               model = input[["model"]], 
               experiment = input[["experiment"]],
               ensemble = input[["ensemble"]]) 
      
    }) %>% 
      # Concatenate all of the regional data for a single ESM into a single data frame.
      bind_rows()  -> 
      regional_data
    
    # Save the regional averages for the ESM.
    write.csv(regional_data, file = file.path(INTER_OUTPUT_CSV, paste0(base_name, '.csv')), row.names = TRUE)
    
    # Return the regional weighted average data. 
    regional_data
    
    
  }) %>%  
    # Format all of the regional data for all of the ESMs into a single data frame.
    bind_rows()
  
  
}


# 3B. Weighted Global Average -------------------------------------------------------
INTER_OUTPUT_CSV <- file.path(OUTPUT_DIR, 'intermediate-csv-global')
dir.create(INTER_OUTPUT_CSV)

cdo_regional_means(data_LandArea, NULL) %>% 
  write.csv(file = file.path(OUTPUT_DIR, 'cmip6_global_means.csv'), row.names = FALSE)

# 34. Weighted Regional Average -------------------------------------------------------
# Split up into regions every 5 degrees of lattitude.
cdo_regional_means(data_LandArea, list('lat' = 5, 'lon' = NULL)) %>% 
  write.csv(file = file.path(OUTPUT_DIR, 'cmip6_region_means.csv'), row.names = FALSE)


# 4. Intermediate Files -------------------------------------------------------------- 

## If for some reason there is a need to use the intermediatei files un comment this section. 
if(CELAN_UP){
  
  file.remove(list.files(INTERMED_DIR, full.names = TRUE, recursive = TRUE))
  
} else{
  
  list.files(INTER_OUTPUT_CSV, pattern = '.csv', full.names = TRUE) %>% 
    lapply(FUN = read.csv, stringsAsFactors = FALSE) %>% 
    bind_rows %>% 
    select(-X) -> 
    data2 
  
  data2 %>% 
    write.csv(file = file.path(OUTPUT_DIR, 'cmip6_region_meansNEW.csv'), row.names = FALSE)
  
  
}




