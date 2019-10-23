## 0.Find_CMIP5_files.R
## Find the CMIP6 gppLut, rhLitter, raRoot, and rhSoil netcdf files. After the CMIP6 archive index is 
## in place this script will be obsolete. 

## Creates
##  0.to_process_CMIP5.csv - a table of all the data CMIP6 files that need to be processed.
##  0.fx_files_CMIP5.csv - a table of all the meta files (cell area and land fraction) that are needed to generate the land area weights. 

# 0. Set Up --------------------------------------------------------------------------------------
# Load the required R packages. 
library(dplyr)
library(tidyr)
library(tibble)

# Define directories.  
CMIP_DIR    <- '/pic/projects/GCAM/CMIP5-KDorheim' # Where the CMIP6 data lives on PIC.
PROJECT_DIR <- '/pic/projects/GCAM/Dorheim/Dorheim/GlobalC' # Where the GlobalC project lives on pic 
OUTPUT_DIR  <- file.path(PROJECT_DIR, 'cmip5', 'output'); dir.create(OUTPUT_DIR, recursive = TRUE)

# Define the variable that we are interested in looking for. 
vars <- c("gpp", "rh")

# 1. Find the Files To Process -------------------------------------------------------------------
# Search the CMIP6 directory for the CMIP6 variable files. 
files <- list.files(path = CMIP_DIR, pattern = paste(vars, collapse = '|'), full.names = TRUE, recursive = TRUE)

# Format the files to process into a data frame that contains information about the file path, variable, domain, 
# model, experiment, ensemble, grid, and time. This information will be saved as a csv file that 
# will be used in the 1.Prcoess_GlobalC_Netcdf files. 
tibble(files = files) %>% 
  mutate(basename = gsub(pattern = '.nc', replacement = '', x = basename(files))) %>% 
  separate(col = basename, into = c('variable', 'domain', 'model', 'experiment', 'ensemble', 'time'), sep = '_') %>% 
  mutate(era = 'cmip5') -> 
  to_process

write.csv(to_process, file = file.path(OUTPUT_DIR, '0.to_process_CMIP5.csv'), row.names = FALSE)


# 2. Find fx (meta) files -------------------------------------------------------------------

# Define the variable that we are interested in looking for. 
fx_vars <- c("sftlf", "areacella")

#Search the CMIP6 directory for the CMIP6 variable files. 
files <- list.files(path = file.path(CMIP_DIR, 'fx'), pattern = paste(fx_vars, collapse = '|'), 
                    full.names = TRUE, recursive = TRUE)

# Format the files to process into a data frame that contains information about the file path, variable, domain, 
# model, experiment, ensemble, grid, and time. This information will be saved as a csv file that 
# will be used in the 1.Prcoess_GlobalC_Netcdf files. 
tibble(files = files) %>% 
  mutate(basename = gsub(pattern = '.nc', replacement = '', x = basename(files))) %>% 
  separate(col = basename, into = c('variable', 'meta', 'model', 'experiment', 'ensemble'), sep = '_') %>% 
  spread(variable, files) %>% 
  na.omit %>% 
  filter(model %in% to_process$model & experiment %in% to_process$experiment) %>% 
  mutate(era = 'cmip5') ->
  fx_table

write.csv(fx_table, file = file.path(OUTPUT_DIR, '0.fx_files_CMIP5.csv'), row.names = FALSE)






