
## Thornthwaite's Aridity Index (TAI)
## By: H. Achicanoy
## December, 2022
## Edited by: F. Castro-Llanos
## March, 2025

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate,envirem))

#source('https://raw.githubusercontent.com/CIAT-DAPA/agro-clim-indices/main/_main_functions.R')

root <- '/home/jovyan/common_data'
ref <- terra::rast('/home/jovyan/common_data/chirps_wrld/chirps-v2.0.1981.01.01.tif')
# sce_climate <- "historical" #"future"terra::rast('/home/jovyan/common_data/chirps_wrld/chirps-v2.0.1981.01.01.tif'

# ET SRAD, load and process only once
srf <- list.dirs(paste0(root,'/ET_SolRad'), full.names = T, recursive = F)
srf <- srf[-length(srf)]
srf <- srf |> gtools::mixedsort()
srd <- srf |> terra::rast()
srd <- srd |> terra::resample(ref)
srd <- srd |> terra::crop(ref)
srd <- srd |> terra::mask(ref)
names(srd) <- c(paste0('SRAD_0',1:9),paste0('SRAD_', 10:12))

cat('Srad read\n')

# Calculate TAI function
calc_tai <- function(yr){
  
  outfile <- paste0(out_dir,'/TAI-',yr,'.tif')
  cat(outfile, "\n")
  
  if(!file.exists(outfile)){
    
    # To create the output raster
    dir.create(dirname(outfile),F,T)
    
    # Sequence of dates
    dts <- seq(from = as.Date(paste0(yr,'-01-01')), to = as.Date(paste0(yr,'-12-31')), by = 'day')
    
    # Files
    pr_fls <- paste0(pr_pth,'/pr_', dts,'.tif')
    pr_fls <- pr_fls[file.exists(pr_fls)]
    
    tx_fls <- paste0(tx_pth, '/tasmax_', dts, '.tif')
    tx_fls <- tx_fls[file.exists(tx_fls)]
    
    tm_fls <- paste0(tm_pth, '/tasmin_', dts, '.tif')
    tm_fls <- tm_fls[file.exists(tm_fls)]
    
    # Read variables, and calculate monthly summaries, do by month to reduce memory consumption
    prc_ls <- tav_ls <- rng_ls <- c()
    
    for (j in 1:12) {
      
      cat("month=", j, "\n")
      
      # Load precip
      mth_fls <- paste0(pr_pth, '/pr_', yr, '-', sprintf("%02.0f", j), "-", sprintf("%02.0f", 1:31), '.tif')
      prc <- terra::rast(pr_fls[pr_fls %in% mth_fls])
      prc <- prc %>% terra::crop(terra::ext(ref))
      prc_month <- sum(prc, na.rm=TRUE)
      rm(prc)
      gc(verbose=F, full=T, reset=T)
      
      # Load tmax
      mth_fls <- paste0(tx_pth, '/tasmax_', yr, '-', sprintf("%02.0f", j), "-", sprintf("%02.0f", 1:31), '.tif')
      tmx <- terra::rast(tx_fls[tx_fls %in% mth_fls])
      tmx <- terra::crop(tmx, ext(ref))
      
      ## Load tmin
      mth_fls <- paste0(tm_pth, '/tasmin_', yr, '-', sprintf("%02.0f", j), "-", sprintf("%02.0f", 1:31), '.tif')
      tmn <- terra::rast(tm_fls[tm_fls %in% mth_fls])
      tmn <- terra::crop(tmn, ext(ref))
      
      ## Calculate average temperature 
      tav <- (tmx + tmn) / 2
      tav_month <- mean(tav, na.rm = T)
      rm(tav)
      gc(verbose = F, full = T, reset = T)
      
      ## Calculate temperature range
      rnge <- abs(tmx - tmn)
      rng_month <- mean(rnge, na.rm=TRUE)
      rm(rnge, tmx, tmn)
      gc(verbose=F, full=T, reset=T)
      
      # Append
      prc_ls <- c(prc_ls, prc_month)
      tav_ls <- c(tav_ls, tav_month)
      rng_ls <- c(rng_ls, rng_month)
      
    }
    
    # Assign precipitation names in envirem environment
    envirem::assignNames(solrad='SRAD_##', tmean = 'TMEAN_##', precip = 'PREC_##')
    
    TMEAN <- tav_ls %>% terra::rast() %>% raster::stack()
    PREC  <- prc_ls %>% terra::rast() %>% raster::stack()
    TRNG  <- rng_ls %>% terra::rast() %>% raster::stack()
    
    names(TMEAN) <- c(paste0('TMEAN_0',1:9),paste0('TMEAN_', 10:12))
    names(PREC)  <- c(paste0('PREC_0',1:9),paste0('PREC_', 10:12))
    names(TRNG)  <- c(paste0('TRNG_0',1:9),paste0('TRNG_', 10:12))
    
    #clean up
    rm(tav_ls, prc_ls, rng_ls, tav_month, prc_month, rng_month)
    gc(verbose=F, full=T, reset=T)
    
    # Thornthwaite's Aridity Index
    srd <- terra::resample(srd, rast(TMEAN))
    srd <- terra::crop(srd, ext(ref))
    
    PET <- envirem::monthlyPET(rast(TMEAN), srd, rast(TRNG))
    names(PET)  <- c(paste0('PET_0',1:9), paste0('PET_', 10:12))
    TAI <- envirem::aridityIndexThornthwaite(rast(PREC), PET)
    names(TAI) <- yr
    
    #write output
    terra::writeRaster(TAI, outfile)
    
    #clean up
    rm(PET, TMEAN, PREC, TRNG, TAI)
    gc(verbose=F, full=T, reset=T)
    
  }
}


# Future setup
#gcm <- 'EC-Earth3' #'ACCESS-ESM1-5' 'MPI-ESM1-2-HR' 'EC-Earth3' 'INM-CM5-0' 'MRI-ESM2-0'

for (gcm in c("ACCESS-ESM1-5", "MPI-ESM1-2-HR", "EC-Earth3", "INM-CM5-0", "MRI-ESM2-0")) {
  
  for (ssp in c('ssp126', 'ssp245', 'ssp370', 'ssp585')) {
    
    # ssp <- 'ssp585'
    # prd <- '2021_2040'
    # gcm <- 'ACCESS-ESM1-5'
    
    cmb <- paste0(ssp, '_', gcm)
    cat(cmb, '\n')
    
    pr_pth <- paste0(root,'/nex-gddp-cmip6/pr/', ssp, '/', gcm) # Precipitation
    tm_pth <- paste0(root,'/nex-gddp-cmip6/tasmin/', ssp, '/', gcm) # Minimum temperature
    tx_pth <- paste0(root,'/nex-gddp-cmip6/tasmax/', ssp, '/', gcm) # Maximum temperature
    out_dir <- paste0(root,'/nex-gddp-cmip6_indices/', cmb, '/TAI')
    
    2021:2100 %>%
      purrr::map(.f = function(i){calc_tai(yr = i)})
    
  }
}

