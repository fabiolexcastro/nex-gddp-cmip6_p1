## Number of soil water stress days (NDWS)
## By: H. Achicanoy
## December, 2022

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))
source('https://raw.githubusercontent.com/CIAT-DAPA/agro-clim-indices/main/_main_functions.R')

root <- '/home/jovyan/common_data'
ref <- terra::rast('/home/jovyan/common_data/chirps_wrld/chirps-v2.0.1981.01.01.tif')

# Soil variables
scp <- terra::rast(paste0(root,'/atlas_hazards/soils/sscp_world.tif'))
scp <- scp %>% terra::resample(ref) %>% terra::mask(ref)
sst <- terra::rast(paste0(root,'/atlas_hazards/soils/ssat_world.tif'))
sst <- sst %>% terra::resample(ref) %>% terra::mask(ref)

# Calculate NDWS function
calc_ndwl0 <- function(yr, mn){
  outfile <- paste0(out_dir,'/NDWL0-',yr,'-',mn,'.tif')
  cat(outfile, "\n")
  if(!file.exists(outfile)){
    dir.create(dirname(outfile),F,T)
    # Last day of the month
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    # Sequence of dates
    if(as.numeric(yr) > 2020 & mn == '02'){
      dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-28')), by = 'day')
    } else {
      dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    }
    
    # Files
    pr_fls <- paste0(pr_pth,'/pr_', dts,'.tif')
    pr_fls <- pr_fls[file.exists(pr_fls)]
    
    tx_fls <- paste0(tx_pth, '/tasmax_', dts, '.tif')
    tx_fls <- tx_fls[file.exists(tx_fls)]
    
    tm_fls <- paste0(tm_pth, '/tasmin_', dts, '.tif')
    tm_fls <- tm_fls[file.exists(tm_fls)]
    
    sr_fls <- paste0(sr_pth, '/rsds_', dts, '.tif')
    sr_fls <- sr_fls[file.exists(sr_fls)]
    
    # Read variables
    prc <- terra::rast(pr_fls)
    prc <- prc %>% terra::crop(terra::ext(ref))
    
    tmx <- terra::rast(tx_fls)
    tmx <- terra::crop(tmx, terra::ext(ref))
    
    tmn <- terra::rast(tm_fls)
    tmn <- terra::crop(tmn, terra::ext(ref))
    
    tav <- (tmx + tmn)/2
    
    srd <- terra::rast(sr_fls)
    srd <- srd %>% terra::crop(terra::ext(ref))
    
    # Maximum evapotranspiration
    ETMAX <- terra::lapp(x = terra::sds(srd,tmn,tav,tmx), fun = peest)
    ETMAX
    rm(list=c("tmn", "tmx", "tav", "srd"))
    gc(verbose=F, full=T, reset=T)
    
    ## 
    scp <- terra::resample(scp, ETMAX, method = 'bilinear')
    sst <- terra::resample(sst, ETMAX, method = 'bilinear')
    
    # Compute water balance model
    date <- paste0(yr,'-',mn)
    if(date %in% c('2021-01')){
      AVAIL <<- ETMAX * 0
      AVAIL[!is.na(AVAIL)] <- 0
    } else {
      avail_fl <- list.files(path=dirname(outfile), pattern="AVAIL-")
      avail_fl <- avail_fl[grep(pattern="\\.tif", avail_fl)]
      avail_fl <- avail_fl[length(avail_fl)]
      AVAIL <<- terra::rast(paste0(dirname(outfile),"/", avail_fl))
      AVAIL <<- AVAIL[[terra::nlyr(AVAIL)]]
    }
    
    eabyep_calc <- function(soilcp = scp, soilsat = ssat, avail = AVAIL,rain = prc[[1]], evap = ETMAX[[1]]){
      
      avail   <- min(avail, soilcp)
      
      # ERATIO
      percwt <- min(avail/soilcp*100, 100)
      percwt <- max(percwt, 1)
      eratio <- min(percwt/(97-3.868*sqrt(soilcp)), 1)
      
      demand  <- eratio * evap
      result  <- avail + rain - demand
      logging <- result - soilcp
      logging <- max(logging, 0)
      logging <- min(logging, soilsat)
      # runoff  <- result - logging + soilcp
      avail   <- min(soilcp, result)
      avail   <- max(avail, 0)
      # runoff  <- max(runoff, 0)
      
      out     <- list(Availability = c(AVAIL, avail),
                      # Demand       = demand,
                      # Eratio       = eratio
                      Logging      = logging
      )
      return(out)
    }
    
    watbal <- 1:terra::nlyr(ETMAX) %>%
      purrr::map(.f = function(i){
        water_balance <- eabyep_calc(soilcp  = scp,
                                     soilsat = sst,
                                     avail   = AVAIL[[terra::nlyr(AVAIL)]],
                                     rain    = prc[[i]],
                                     evap    = ETMAX[[i]])
        AVAIL <<- water_balance$Availability
        return(water_balance)
      })
    
    LOGGING <- watbal %>% purrr::map('Logging') %>% terra::rast()
    
    # Calculate number of soil waterlogging days (if logging is above 0)
    # Note NDWL50 uses ssat * 0.5 (so this means soil is at 50% toward saturation)
    # here it suffices if the soil is above field capacity
    NDWL0  <- sum(LOGGING > 0)
    
    ## To write the raster
    terra::writeRaster(NDWL0, outfile)
    terra::writeRaster(AVAIL, paste0(dirname(outfile),'/AVAIL-',yr,'-',mn,'.tif'))
    
    ## Clean up
    rm(list=c("prc", "ETMAX", "AVAIL", "watbal", "ERATIO", "LOGGING", "NDWL0"))
    gc(verbose=F, full=T, reset=T)
    
  }
  
}

# Future setup
# gcm <- 'ACCESS-ESM1-5' #'ACCESS-ESM1-5' 'MPI-ESM1-2-HR' 'EC-Earth3' 'INM-CM5-0' 'MRI-ESM2-0'
# ssp <- 'ssp126'

# gcm <-  'ACCESS-ESM1-5'
# ssp <- 'ssp126'

for (gcm in c("ACCESS-ESM1-5", "MPI-ESM1-2-HR", "EC-Earth3", "INM-CM5-0", "MRI-ESM2-0")) {
  for (ssp in c('ssp126', 'ssp245', 'ssp370', 'ssp585')) {
    
    cmb <- paste0(ssp,'_',gcm)
    cat('To process -----> ', cmb, '\n')
    yrs <- 2021:2100
    mns <- c(paste0('0',1:9),10:12)
    stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
    names(stp) <- c('yrs','mns')
    stp <- stp %>%
      dplyr::arrange(yrs, mns) %>%
      base::as.data.frame()
    
    pr_pth <- paste0(root,'/nex-gddp-cmip6/pr/', ssp, '/', gcm) # Precipitation
    tm_pth <- paste0(root,'/nex-gddp-cmip6/tasmin/', ssp, '/', gcm) # Minimum temperature
    tx_pth <- paste0(root,'/nex-gddp-cmip6/tasmax/', ssp, '/', gcm) # Maximum temperature
    sr_pth <- paste0(root,'/nex-gddp-cmip6/rsds/',ssp,'/',gcm) # Solar radiation
    out_dir <- paste0(root,'/nex-gddp-cmip6_indices/',cmb,'/NDWL0') # Output directory
    
    yrs_mpg <- data.frame(Baseline = as.character(rep(1995:2014, 4)),
                          Future = as.character(c(2021:2040,2041:2060,
                                                  2061:2080,2081:2100)))
    
    1:nrow(stp) %>%
      purrr::map(.f = function(i){
        calc_ndwl0(yr = stp$yrs[i], mn = stp$mns[i])
        gc(verbose = T, full = T, reset = T)
        if (i%%5 == 0) {
          tmpfls <- list.files(tempdir(), full.names=TRUE)
          1:length(tmpfls) %>% purrr::map(.f = function(k) {system(paste0("rm -f ", tmpfls[k]))})
        }
      })
    
  }
}
