

## Average monthly minimum temperature
## By: H. Achicanoy / J. Ramirez-Villegas
## July, 2023
## Edited by Fabio Castro-Llanos
## March, 2025

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))

# Root folder 
root <- '/home/jovyan/common_data'

# Extent CHIRPS
msk <- terra::rast('/home/jovyan/common_data/chirps_wrld/chirps-v2.0.1981.01.01.tif')
ext <- ext(msk)

# Function
calc_tmax <- function(yr, mn){
  
  #yr <- '2021'; mn <- '01'
  
  outfile <- paste0(out_dir, 'TMAX-', yr, '-', mn, '.tif')
  cat(basename(outfile), '\n')
  
  if(!file.exists(outfile)){
    
    ## Create output directory
    dir.create(dirname(outfile), F, T)
    
    ## Tidy the dates
    cat('Processing -------> ', yr, mn, '\n')
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    fls <- paste0(ts_pth, '/', 'tasmax_', dts, '.tif')
    fls <- fls[file.exists(fls)]
    
    ## Read the precipitation data
    tas <- rast(fls)
    
    ## Read the precipitation data
    tas <- rast(fls)
    tas <- terra::crop(tas, ext)
    
    ## To calculate the maximum temperature (average)
    tas <- mean(tas, na.rm = T)
    terra::writeRaster(x = tas, filename = outfile, overwrite = T)
    rm(tas); gc(reset = T)
    
    
  }
  
}

# Future setup
for(ssp in c('ssp126', 'ssp245', 'ssp370', 'ssp585')){
  for(gcm in c('ACCESS-ESM1-5', 'MPI-ESM1-2-HR', 'EC-Earth3', 'INM-CM5-0', 'MRI-ESM2-0')){
    
    ## Parameters
    cmb <- paste0(ssp, '_', gcm)
    yrs <- 2021:2100
    mnt <- c(paste0('0', 1:9), 10:12)
    stp <- base::expand.grid(yrs, mnt) %>% base::as.data.frame() %>% setNames(c('yrs', 'mnt')) %>% arrange(yrs, mnt) %>% as.data.frame(); rm(yrs, mnt)
    
    ## Setup in/out files
    ts_pth  <- paste0(root, '/nex-gddp-cmip6/tasmax/', ssp, '/', gcm)
    out_dir <- paste0(root, '/nex-gddp-cmip6_indices/', ssp, '_', gcm, '/TMAX/') 
    
    1:nrow(stp) %>% 
      purrr::map(.f = function(i){calc_tmax(yr = stp$yrs[i], mn = stp$mnt[i])})
    
    ##
    cat('----Finish----\n')
    
  }
}

# Check the results 
# fles <- dir_ls('/home/jovyan/common_data/nex-gddp-cmip6_indices/ssp126_ACCESS-ESM1-5/PTOT') %>% as.character()
# length(fles) == (length(2021:2100) * 12) #
