

## Convert nc to tif 
## F Castro-Llanos

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, foreach, doParallel, tidyverse, glue)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# List the files ----------------------------------------------------------
root <- '/home/jovyan/common_data'
gcms <- c('ACCESS-ESM1-5','EC-Earth3','INM-CM5-0','MPI-ESM1-2-HR','MRI-ESM2-0')
ssps <- c('ssp126', 'ssp245', 'ssp370', 'ssp585')

# Function ----------------------------------------------------------------
nc2tif <- function(yr, mn){
  
  # yr  <- '2039'
  # mn  <- '01'
  
  ##
  cat('>>>>>> ', yr, mn, '\n')
  last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
  sr_pth <- paste0(root,'/ecmwf_agera5_cmip6_africa/solar_radiation_flux_',gcm,'_',ssp,'_',prd) 
  
  if(as.numeric(yr) > 2020 & mn == '02'){
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-28')), by = 'day')
  } else {
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
  }
  
  sr_fls <- paste0(sr_pth,'/Solar-Radiation-Flux_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.0.nc')
  sr_fls <- sr_fls[file.exists(sr_fls)]
  
  ## Read the raster ## 0.381, 2.688 
  system.time(sr <- rast(sr_fls))
  if(yr <= 2020) {srd <- srd/1000000}
  
  out_dir <- paste0(root, '/', 'ecmwf_agera5_cmip6_africa_tif', '/', 'solar_radiation_flux_', gcm, '_', ssp, '_', prd)
  dir_create(out_dir)
  out_fls <- paste0(out_dir, '/', gsub('.nc', '.tif', basename(sr_fls)))
  
  ## Clean raster
  system.time(
    clean_r <- rast(
      nrows = nrow(sr),
      ncols = ncol(sr),
      nlyrs = nlyr(sr),
      crs = crs(sr),
      extent = ext(sr),
      vals = values(sr)
    )
  )

  
  # system.time(
    # expr = {
      for(i in 1:length(out_fls)){terra::writeRaster(x = sr[[i]], filename = out_fls[i], overwrite = TRUE)}
    # } 
  # )
  cat('Done!\n')

}

# Loop --------------------------------------------------------------------

gcm <- 'ACCESS-ESM1-5'
ssp <- 'ssp126'
prd <- '2021_2040'

for (gcm in c('ACCESS-ESM1-5')) {
  for (ssp in c('ssp126', 'ssp245', 'ssp370', 'ssp585')) {
    for (prd in c('2021_2040','2041_2060','2061_2080','2081_2100')) {
      
      cmb <- paste0(ssp,'_',gcm,'_',prd)
      prd_num <- as.numeric(unlist(strsplit(x = prd, split = '_')))
      yrs <- prd_num[1]:prd_num[2]
      mns <- c(paste0('0',1:9),10:12)
      stp <- base::expand.grid(yrs, mns) |> base::as.data.frame(); rm(yrs,mns)
      names(stp) <- c('yrs','mns')
      stp <- stp |> dplyr::arrange(yrs, mns) |> base::as.data.frame()
      sr_pth <- paste0(root,'/ecmwf_agera5_cmip6_africa/solar_radiation_flux_',gcm,'_',ssp,'_',prd) # Solar radiation

      system.time(
        expr = {
          1:nrow(stp) |>
            purrr::map(.f = function(i){
              nc2tif(yr = stp$yrs[i], mn = stp$mns[i])
            })
        }
      )
    
    }
  }
}


# End ---------------------------------------------------------------------

# Función para procesar un solo raster
# process_single <- function(i, wrapped_raster, filename) {
#   library(terra)
#   # Deserializar el raster
#   rast <- unwrap(wrapped_raster[[i]])
#   # Escribir el raster
#   writeRaster(
#     x = rast,
#     filename = filename,
#     overwrite = TRUE,
#     wopt = list(
#       gdal = c("COMPRESS=LZW", "PREDICTOR=2"),
#       datatype = "FLT4S"
#     )
#   )
#   return(TRUE)
# }
# 
# # Serializar los rasters individuales
# wrapped_rasters <- map(1:nlyr(sr), function(i){r <- wrap(sr[[i]])})
# 
# # Configurar paralelización
# n_cores <- min(detectCores() - 1, length(out_fls))
# 
# # Ejecutar en paralelo usando mclapply (solo Linux/macOS)
# results <- mclapply(
#   seq_along(out_fls),
#   function(i){process_single(i, wrapped_rasters, out_fls[i])},
#   mc.cores = n_cores
# )
# 
# # Verificar resultados
# if (any(unlist(results) != TRUE)) {
#   warning("Algunos archivos no se procesaron correctamente")
# }
# 
# rm(wrapped_rasters, results)
