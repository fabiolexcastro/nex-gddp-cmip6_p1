## Number of soil water stress days (NDWS)
## By: H. Achicanoy & A. Mendez
## April, 2025

# R options
# args <- commandArgs(trailingOnly = T)
options(warn = -1, scipen = 999) # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate,compiler,raster,ncdf4))

peest2 <- function(srad, tmin, tmean, tmax){
  
  # Convert rasters to matrices for faster processing
  srad_m  <- terra::values(srad) # this variable comes in .nc format (original raster::values)
  tmin_m  <- terra::values(tmin)
  tmean_m <- terra::values(tmean)
  tmax_m  <- terra::values(tmax)
  
  # Constants
  albedo  <- 0.2
  vpd_cte <- 0.7
  a_eslope <- 611.2
  b_eslope <- 17.67
  c_eslope <- 243.5
  
  # Pre-allocate matrices
  n_cells <- ncol(srad_m)
  n_layers <- nrow(srad_m)
  rn <- matrix(0, nrow = n_layers, ncol = n_cells)
  eslope <- matrix(0, nrow = n_layers, ncol = n_cells)
  et_max <- matrix(0, nrow = n_layers, ncol = n_cells)
  
  # Optimized calculations
  rn <- (1-albedo) * srad_m
  rm(srad_m)
  
  temp_denom <- tmean_m + c_eslope
  eslope <- (a_eslope * b_eslope * c_eslope * exp(b_eslope * tmean_m/temp_denom)) / (temp_denom^2)
  rm(tmean_m, temp_denom)
  
  vpd <- vpd_cte * (0.61120 * (exp((17.67*tmax_m)/(tmax_m+243.5)) - 
                                 exp((17.67*tmin_m)/(tmin_m+243.5))))
  rm(tmin_m, tmax_m)
  
  # Constants
  pt_const <- 1.26
  pt_fact  <- 1
  vpd_ref  <- 1
  psycho   <- 62
  rho_w    <- 997
  rlat_ht  <- 2.26E6
  
  # Final calculation
  pt_coef <- 1 + (pt_fact*pt_const-1) * vpd / vpd_ref
  conversion_factor <- 1E6 * 100 / (rlat_ht * rho_w) * 10
  et_max <- pt_coef * rn * eslope/(eslope+psycho) * conversion_factor
  
  # Convert back to raster
  et_max_rast <- tmin
  terra::values(et_max_rast) <- et_max
  
  return(et_max_rast)
  
}

root <- '/home/jovyan/common_data'

# Soil variables
ref <- terra::rast('/home/jovyan/common_data/chirps_wrld/chirps-v2.0.1981.01.01.tif')

# Soil variables
scp <- terra::rast(paste0(root,'/atlas_hazards/soils/sscp_world.tif'))
scp <- scp |> terra::resample(ref) |> terra::mask(ref)
sst <- terra::rast(paste0(root,'/atlas_hazards/soils/ssat_world.tif'))
sst <- sst |> terra::resample(ref) |> terra::mask(ref)

# Calculate NDWL50 function
# yr <- '2021'
# mn <- '01'

calc_ndws <- function(yr, mn){
  outfile <- paste0(out_dir,'/NDWS-',yr,'-',mn,'.tif')
  cat(outfile,'\n')
  if (!file.exists(outfile)) {
    dir.create(dirname(outfile),F,T)
    # Last day of the month
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    # Sequence of dates
    if(as.numeric(yr) > 2020 & mn == '02'){
      dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-28')), by = 'day')
    } else {
      dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    }
    
    cat(">>> Iniciando proceso:", yr, "-", mn, " \n")
    
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
    prc <- prc |> terra::crop(terra::ext(ref))
    
    tmx <- terra::rast(tx_fls)
    tmx <- terra::crop(tmx, terra::ext(ref))
    
    tmn <- terra::rast(tm_fls)
    tmn <- terra::crop(tmn, terra::ext(ref))
    
    tav <- (tmx + tmn)/2
    
    srd <- terra::rast(sr_fls)
    srd <- srd |> terra::crop(terra::ext(ref))
    
    # Maximum evapotranspiration
    ETMAX <<- peest2(srad = srd, tmin = tmn, tmean = tav, tmax = tmx)
    rm(list = c('tmn','tmx','tav','srd'))
    
    scp <- terra::resample(scp, ETMAX[[1]], method = 'bilinear')
    sst <- terra::resample(sst, ETMAX[[1]], method = 'bilinear')
    
    # Compute water balance model
    date <- paste0(yr,'-',mn)
    if(date %in% c('2021-01')){
      AVAIL <<- ETMAX[[1]] * 0
      AVAIL[!is.na(AVAIL)] <- 0
    } else {
      avail_fl <- list.files(path=dirname(outfile), pattern="AVAIL-")
      avail_fl <- avail_fl[grep(pattern="\\.tif", avail_fl)]
      avail_fl <- avail_fl[length(avail_fl)]
      AVAIL <<- terra::rast(paste0(dirname(outfile),"/", avail_fl))
      AVAIL <<- AVAIL[[terra::nlyr(AVAIL)]]
    }
    eabyep_calc <- function(soilcp = scp, soilsat = ssat, avail = AVAIL, rain = prc[[1]], evap = ETMAX[[1]]){
      
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
                      Eratio       = eratio
                      # Logging      = logging
      )
      return(out)
    }
    ceabyep_calc <- compiler::cmpfun(eabyep_calc)
    
    watbal <- vector("list", terra::nlyr(ETMAX))
    for(j in 1:terra::nlyr(ETMAX)){
      water_balance <- ceabyep_calc(soilcp  = scp,
                                    soilsat = sst,
                                    avail   = AVAIL[[terra::nlyr(AVAIL)]],
                                    rain    = prc[[j]],
                                    evap    = ETMAX[[j]])
      # Update AVAIL with deep copy to avoid memory leaks
      AVAIL <<- terra::deepcopy(water_balance$Availability)
      # Store result and clean temporary objects
      watbal[[j]] <- water_balance
      rm(water_balance)
    }
    
    ERATIO <- watbal %>% purrr::map('Eratio') %>% terra::rast()
    
    # Calculate number of soil water stress days
    cvls <- matrix(data = c(-Inf, 0.5, 1), ncol = 3) # Classification values
    NDWS <- terra::classify(x = ERATIO, rcl = cvls, right = F) |> sum()
    terra::writeRaster(NDWS, outfile)
    terra::writeRaster(AVAIL, paste0(dirname(outfile),'/AVAIL-',yr,'-',mn,'.tif'))
  }
}

# ssps  <- args[1]
# gcms  <- args[2]
# prds  <- args[3]
#c('ACCESS-ESM1-5','EC-Earth3','INM-CM5-0','MPI-ESM1-2-HR','MRI-ESM2-0')
# Rscript fast_calc_NDWL50.R 'ssp245' 'ACCESS-ESM1-5' 'EC-Earth3'                   screen id = ndwl50
# Rscript fast_calc_NDWL50.R 'ssp245' 'INM-CM5-0' 'MPI-ESM1-2-HR'  'MRI-ESM2-0'     screen id = ndwl50_2

gcms <- c('ACCESS-ESM1-5','EC-Earth3','INM-CM5-0','MPI-ESM1-2-HR','MRI-ESM2-0')
ssps <- c('ssp126', 'ssp245', 'ssp370', 'ssp585')

#gcm <- gcms[1]
#ssp <- ssps[1]

#'2021_2040','2041_2060',
# Future setup
for (gcm in gcms) {
  for (ssp in ssps) {
    
    cmb <- paste0(ssp,'_',gcm)
    cat('To process -----> ', cmb, '\n')
    yrs <- 2021:2100
    mns <- c(paste0('0',1:9),10:12)
    stp <- base::expand.grid(yrs, mns) |> base::as.data.frame(); rm(yrs,mns)
    names(stp) <- c('yrs','mns')
    stp <- stp |>
      dplyr::arrange(yrs, mns) |>
      base::as.data.frame()
    
    pr_pth <- paste0(root,'/nex-gddp-cmip6/pr/', ssp, '/', gcm) # Precipitation
    tm_pth <- paste0(root,'/nex-gddp-cmip6/tasmin/', ssp, '/', gcm) # Minimum temperature
    tx_pth <- paste0(root,'/nex-gddp-cmip6/tasmax/', ssp, '/', gcm) # Maximum temperature
    sr_pth <- paste0(root,'/nex-gddp-cmip6/rsds/',ssp,'/',gcm) # Solar radiation
    out_dir <- paste0(root,'/nex-gddp-cmip6_indices/',cmb,'/NDWL50') # Output directory
    
    1:nrow(stp) |>
      purrr::map(.f = function(i){
        calc_ndws(yr = stp$yrs[i], mn = stp$mns[i])
        gc(verbose = T, full = T, reset = T)
        if (i%%5 == 0) {
          tmpfls <- list.files(tempdir(), full.names=TRUE)
          1:length(tmpfls) |> purrr::map(.f = function(k) {system(paste0("rm -f ", tmpfls[k]))})
        }
      })
    
  }
}
