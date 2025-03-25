

## Check trend for ptot / tmin / tmax
## By: F Castro - Llanos
## March, 2025

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(fs,tidyverse,terra,cowplot,gtable,gridExtra,gtools,patchwork,lubridate,rnaturalearthdata,rnaturalearth,raptr))

# Data
root <- '/home/jovyan/common_data/nex-gddp-cmip6_indices'
dirs <- as.character(dir_ls(root, type = 'directory'))
length(dirs) == 5 * 4 # 5 GCMs / 4 SSPs

# Datamask 
wrld <- ne_countries(returnclass = 'sf', scale = 50)
wrld <- vect(wrld)
mskr <- terra::rast('/home/jovyan/common_data/nex-gddp-cmip6_indices/ssp126_ACCESS-ESM1-5/PTOT/PTOT-2021-01.tif', lyr = 1)
mskr <- terra::crop(mskr, wrld)
mskr <- terra::mask(mskr, wrld)
mskr <- mskr * 0 

# Sample points
mskr.zone <- wrld %>% st_as_sf() %>% filter(sov_a3 %in% c('COL', 'ECU', 'PER') | region_un == 'Africa')
mskr.zone <- terra::crop(mskr, mskr.zone) %>% terra::mask(., mskr.zone)

pnts <- raptr::randomPoints(mskr.zone, n = 10, prob = FALSE)
pnts <- as.data.frame(pnts)
# pnts$id <- 1:nrow(pnts)

plot(mskr.zone)
points(pnts$x, pnts$y, pch = 16)
write.csv(pnts, './tbl/samplePoints.csv', row.names = F)

pnts <- read_csv('./tbl/samplePoints.csv')

## SSPs / GCMs / Years
ssps <- c('ssp126', 'ssp245', 'ssp370', 'ssp585')
gcms <- c('ACCESS-ESM1-5', 'EC-Earth3', 'INM-CM5-0', 'MPI-ESM1-2-HR', 'MRI-ESM2-0')
year <- 2021:2100

# Function to use 
to.check <- function(ssp, vrb){
  
  # ssp <- 'ssp126'
  # vrb <- 'TMIN'
  
  ## To list the directorie
  cat('To process: ', ssp, vrb, '\n')
  dire <- grep(ssp, dirs, value = T)
  dire <- paste0(dire, '/', vrb)
  
  ## Conditional
  if(any(dir.exists(dire))){
    
    cat('Continuous\n')
    
    fnal <- list()
    
    for(gcm in gcms){
      
      cat('-------------> ', gcm, '\n')
      fles <- grep(gcm, dire, value = T) |>
        dir_ls() |>
        as.character()
      
      dfrm <- list()
      
      for(yr in year){
        
        rstr <- grep(paste0('-', yr, '-'), fles, value = T)
        rstr <- rast(rstr)  
        names(rstr) <- paste0(vrb, '_', yr, '-', 1:12)
        
        vles <- terra::extract(rstr, pnts[,1:2])
        vles <- as.data.frame(cbind(pnts, vles))
        # vles <- mutate(vles, ID = 1:nrow(vles))
        vles <- dplyr::select(vles, -ID)
        # vles <- rename(vles, id = ID)
        
        dfrm[[as.character(yr)]] <- vles
        
      }      
      
      dfrm <- dfrm |> reduce(inner_join, by = c('x', 'y', 'id'))
      dfrm <- as_tibble(dfrm)
      dfrm <- gather(dfrm, var, value, -c(x, y, id))
      dfrm <- separate(data = dfrm, col = 'var', into = c('variable', 'date'), sep = '_')
      dfrm <- mutate(dfrm, year = str_sub(date, 1, 4), month = str_sub(date, 6, nchar(date)))
      dfrm <- mutate(dfrm, month = as.numeric(month))
      dfrm <- mutate(dfrm, month = ifelse(month < 10, paste0('0', month), as.character(month)))
      dfrm <- mutate(dfrm, date = paste0(year, '-', month, '-', '01'), date = as.Date(date, format = '%Y-%m-%d'))
      dfrm <- dplyr::select(dfrm, x, y, variable, date, value)
      dfrm <- mutate(dfrm, sspe = ssp, gcme = gcm, .before = 'x')
      
      fnal[[gcm]] <- dfrm
        
    }
    
    cat('------ Finish! --------\n')
    
    
  } else {
    cat('Folder doesnt exists')
  }
  
  return(fnal)

}

# To apply the function
rltn <- expand.grid(a.ssps = ssps, a.vars = c('PTOT', 'TMIN', 'TMAX')) 
rltn <- mutate(rltn, a.ssps = as.character(a.ssps), a.vars = as.character(a.vars))
rltn.prec <- filter(rltn, a.vars == 'PTOT')
rltn.tmin <- filter(rltn, a.vars == 'TMIN')
rltn.tmax <- filter(rltn, a.vars == 'TMAX')

# Prec --------------------------------------------------------------------

##
pnts.prec <- 1:nrow(rltn.prec) |> map(.f = function(k){to.check(ssp = rltn.prec$a.ssps[k], vrb = rltn.prec$a.vars[k])})
pnts.prec <- bind_rows(map(pnts.prec, bind_rows))
write.csv(pnts.prec, './tbl/prec_sample-vles.csv', row.names = F)
pnts.prec <- read_csv('./tbl/prec_sample-vles.csv', show_col_types = FALSE)

## To add the country 
nmes <- terra::extract(wrld, pnts[,c(1,2)])$admin
pnts <- mutate(pnts, country = nmes)
pnts.prec <- inner_join(pnts.prec, as.data.frame(pnts), by = c('x', 'y'))

## Function for making the graph
make.graph <- function(tble, ssp, gid){
  
  # tble <- pnts.prec; ssp <- 'ssp126'; gid <- 1
  
  ## To filter the dataset
  tbl <- filter(tble, id == gid & sspe == ssp)
  gcs <- unique(tbl$gcme)
  
  ## To make decompose funcion
  trnd <- map_dfr(.x = 1:length(gcs), .f = function(i){
    
    trn <- tbl %>%
      filter(gcme == gcs[i]) %>% 
      pull(value) %>% 
      ts(start = c(2021, 1), end = c(2100, 12), frequency = 12) %>% 
      decompose() %>%
      .$trend
    df_trend <- tibble(gcm = gcs[i], date = seq.Date(from = as.Date("2021-01-01"), by = "month",length.out = length(trn)), trend = as.numeric(trn))
    return(df_trend)
    
  })
  
  ## To build the graph
  ggl <- ggplot(data = trnd, aes(x = date, y = trend)) +
    geom_line()  +
    facet_wrap(.~gcm) + 
    ggtitle(label = paste0(unique(tbl$country))) +
    labs(y = unique(tbl$variable), x = '') +
    theme_minimal() +
    theme(
      strip.text = element_text(face = 'bold'), 
      plot.title = element_text(face = 'bold', hjust = 0.5),
      axis.text.y = element_text(angle = 90, hjust = 0.5)
    )
  
  ## To save the graph
  out <- './png/graphs_trend'
  ggsave(plot = ggl, filename = paste0(out, '/', 'trend_', gid, '_', unique(tbl$country), '_', unique(tbl$variable), '.jpg'), units = 'in', width = 9, height = 7, dpi = 300)
  
  
}

## To make the graphs
map(.x = 1:length(ssps), .f = function(s){
  map(.x = 1:10, .f = function(g){
    make.graph(tble = pnts.prec, ssp = ssps[s], gid = g)
  })
})

# Tmin ---------------------------------------------------------------------
pnts <- read_csv('./tbl/samplePoints.csv', show_col_types = FALSE)
pnts.tmin <- 1:nrow(rltn.tmin) |> map(.f = function(k){to.check(ssp = rltn.tmin$a.ssps[k], vrb = rltn.tmin$a.vars[k])})
pnts.tmin <- bind_rows(map(pnts.tmin, bind_rows))
write.csv(pnts.tmin, './tbl/tmin_sample-vles.csv', row.names = FALSE)

## To add the country 
nmes <- terra::extract(wrld, pnts[,c(1,2)])$admin
pnts <- mutate(pnts, country = nmes)
pnts.tmin <- inner_join(pnts.tmin, as.data.frame(pnts), by = c('x', 'y'))

## To make the graphs
map(.x = 1:length(ssps), .f = function(s){
  map(.x = 1:10, .f = function(g){
    make.graph(tble = pnts.tmin, ssp = ssps[s], gid = g)
  })
})

# Tmax --------------------------------------------------------------------
pnts <- read_csv('./tbl/samplePoints.csv', show_col_types = FALSE)
pnts.tmax <- 1:nrow(rltn.tmax) |> map(.f = function(k){to.check(ssp = rltn.tmax$a.ssps[k], vrb = rltn.tmax$a.vars[k])})
pnts.tmax <- bind_rows(map(pnts.tmax, bind_rows))
write.csv(pnts.tmax, './tbl/tmax_sample-vles.csv', row.names = FALSE)

## To add the country 
nmes <- terra::extract(wrld, pnts[,c(1,2)])$admin
pnts <- mutate(pnts, country = nmes)
pnts.tmax <- inner_join(pnts.tmax, as.data.frame(pnts), by = c('x', 'y'))

## To make the graphs
map(.x = 1:length(ssps), .f = function(s){
  map(.x = 1:10, .f = function(g){
    make.graph(tble = pnts.tmax, ssp = ssps[s], gid = g)
  })
})
