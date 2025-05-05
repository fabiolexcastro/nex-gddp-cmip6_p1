
## Check trend for ptot / tmin / tmax
## By: F Castro - Llanos
## March, 2025

# R options -----------------------------------------------------------------
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(fs,tidyverse,glue,terra,cowplot,gtable,gridExtra,gtools,patchwork,lubridate,rnaturalearthdata,rnaturalearth,raptr))

# Data ----------------------------------------------------------------------
pnts <- read_csv('./tbl/samplePoints.csv', show_col_types = FALSE)
pnts <- pnts[c(2, 3, 4, 9, 10),]
wrld <- ne_countries(returnclass = 'sf', scale = 50)
pnts <- mutate(pnts, iso = terra::extract(vect(wrld), pnts[,1:2])[,c('sov_a3')])

root <- '/home/jovyan/common_data/nex-gddp-cmip6_indices'
dirs <- as.character(dir_ls(root, type = 'directory'))
indx <- c('NTX30', 'NTX35', 'NDD')
dirs
map(dirs, dir_ls)

# Parameters 
ssps <- c('ssp126', 'ssp245', 'ssp370', 'ssp585')
gcms <- c('ACCESS-ESM1-5', 'EC-Earth3', 'INM-CM5-0', 'MPI-ESM1-2-HR', 'MRI-ESM2-0')

# Function for extract the values ------------------------------------------ 
to.check <- function(ind, ssp, gcm){
  
  # ind <- 'NTX30'
  # ssp <- 'ssp126'
  # gcm <- 'ACCESS-ESM1-5'
  
  ## To list the files
  cat('To process: ', ind, ' ', ssp, ' ', gcm, '\n')
  dir <- dir_ls(dirs) %>% grep(paste0(ssp, '_', gcm), ., value = T) %>% grep(ind, ., value = T) %>% as.character()
  fls <- dir_ls(dir) %>% as.character()
  dts <- basename(fls) %>% gsub('.tif$', '', .)
  
  ## To read as a raster
  rst <- terra::rast(fls)
  names(rst) <- dts
  
  ## To extract the values
  vls <- terra::extract(rst, pnts[,c(1:2)])
  
  ## Tidy the table
  vls <- cbind(pnts, vls) %>% gather(var, value, -c(x, y, id, iso, ID)) %>% as_tibble()
  vls <- mutate(vls, date = gsub(paste0(ind, '-'), '', var))
  vls <- mutate(vls, date = paste0(date, '-01'), date = as.Date(date, format = '%Y-%m-%d'))
  vls <- mutate(vls, gcme = gcm, sspe = ssp)
  
  ## To return 
  cat('Done!\n')
  return(vls)
  
}

# To apply the function ---------------------------------------------------- 

## NTx35
ntx35 <- map_dfr(.x = 1:length(ssps), .f = function(s){
  map_dfr(.x = 1:length(gcms), .f = function(g){
    to.check(ind = 'NTX35', ssp = ssps[s], gcm = gcms[g])
  })
})
write.csv(ntx35, './tbl/index-values/ntx35.csv', row.names = FALSE)

## NTx30
ntx30 <- map_dfr(.x = 1:length(ssps), .f = function(s){
  map_dfr(.x = 1:length(gcms), .f = function(g){
    to.check(ind = 'NTX30', ssp = ssps[s], gcm = gcms[g])
  })
})
write.csv(ntx30, './tbl/index-values/ntx30.csv', row.names = FALSE)

## NDD 
ndd <- map_dfr(.x = 1:length(ssps), .f = function(s){
  map_dfr(.x = 1:length(gcms), .f = function(g){
    to.check(ind = 'NDD', ssp = ssps[s], gcm = gcms[g])
  })
})
write.csv(ndd, './tbl/index-values/ndd.csv', row.names = FALSE)

# To make the graph  ------------------------------------------------------

make.graph <- function(tbl){
  
  tbl <- ntx30
  
  ##
  cat('To process!\n')
  idx <- pull(tbl, var) %>% str_split('-') %>% map_chr(1) %>% unique()
  
  ##
  rslt <- map_dfr(.x = ssps, .f = function(s){
    
    rslt <- map_dfr(.x = gcms, .f = function(g){
      
      dfm <- tbl %>% filter(sspe == s & gcme == g)
      
      rsl <- map_dfr(.x = unique(pull(dfm, iso)), .f = function(i){
        
        sbd <- dfm %>% 
          filter(iso == i) %>% 
          pull(value) %>% 
          ts(start = c(2021, 1), end = c(2100, 12), frequency = 12) %>% 
          decompose() %>% 
          .$trend
        sbd <- tibble(iso = i, gcm = g, ssp = s, date =  seq.Date(from = as.Date("2021-01-01"), by = "month", length.out = length(sbd)), trend = sbd)
        return(sbd)
        
      })
      
      return(rsl)
      
    })
    
  })
  
  ## 
  ggl <- ggplot(data = rslt, aes(x = date, y = trend, col = gcm)) + 
    geom_line() +
    facet_wrap(iso~ssp, ncol = 4, nrow = 5, scales = 'free_y') + 
    labs(x = 'Date', y = 'NDD', col = '') +
    theme_minimal() +
    theme(
      legend.position = 'bottom',
      axis.text.xy = element_text(angle = 90, hjust = 0.5), 
      strip.text = element_text(face = 'bold', hjust = 0.5)
    )
  
  ##
  ggsave(plot = ggl, filename = glue('./png/graphs_trend/{idx}.jpg'), units = 'in', width = 14, height = 12, dpi = 300, create.dir = TRUE)
  cat('Done!\n')
  
}






