


## 
require(pacman)
pacman::p_load(terra, sf, fs, tidyverse, glue)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# List the files ----------------------------------------------------------
fles <- dir_ls('/home/jovyan/common_data/nex-gddp-cmip6')
fles <- as.character(fles)

## Precipitation 
dirs.pr <- grep('pr', fles, value = T) %>% dir_ls() %>% as.character()

pr1 <- dir_ls(dirs.pr[1]) %>% map(., dir_ls) %>% unlist() %>% as.character()
pr1 <- rast(pr1[1])

system.time(expr = {pr <- pr1 * 86400; pr <- terra::rotate(pr)})
system.time(expr = {pr <- terra::rotate(pr1); pr <- pr * 86400})


rs <- list()

system.time(
  expr = {
    
    for(i in 1:nlyr(pr1)){
      
      rs[[i]] <- pr1[[i]] * 86400
      rs[[i]] <- terra::rotate(rs[[i]])
      
    }    
  }  

)

