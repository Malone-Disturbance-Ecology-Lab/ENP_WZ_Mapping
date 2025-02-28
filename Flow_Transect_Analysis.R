
# This function is designed to turn vegetation zones into polygons.

library(terra)
library(tidyverse)
library(sf)

setwd( '/Volumes/MaloneLab/Research/Whitezone/Result/Emulation/Final/CDthreshold06/CD_closest_GE')
rast.files <- list.files( '/Volumes/MaloneLab/Research/Whitezone/Result/Emulation/Final/CDthreshold06/CD_closest_GE', pattern = "tif")

coast <- rast( rast.files) %>% terra::project("+proj=longlat +datum=WGS84")
terra::time(coast) <-  substr(rast.files, 8, 11) %>% as.numeric # get the time out of the list of file names

coast.factor <- coast %>% as.factor
# make all zeros NA
coast.factor [ coast.factor  == 0] <- NA

# Need individual rast layers
coast.1 <- coast.factor
coast.1[ coast.1  != 2 ] <- NA
coast.1[[1]] %>% plot

coast.2 <- coast.factor
coast.2[ coast.2  != 2 ] <- NA
coast.2[[1]] %>% plot

coast.3 <- coast.factor
coast.3[ coast.3  != 3 ] <- NA
coast.3[[1]] %>% plot

coast.1.sf <- polygonize_SpatRast(stack = coast.1)
coast.2.sf <- polygonize_SpatRast(stack = coast.2)
coast.3.sf <- polygonize_SpatRast(stack = coast.3)

ENP.transects <- makeTransect(line1 = coast.sf, line2 =  inland.sf, transect.n= 100)

Int.coastal.1 <- Trans_Intersection(ENP.transects,coast.1.sf  )
Int.coastal.2 <- Trans_Intersection(ENP.transects,coast.2.sf  )
Int.coastal.3 <- Trans_Intersection(ENP.transects,coast.3.sf  )

Int.coastal.1 %>% ggplot() + geom_sf() + facet_wrap(~time)

Int.coastal.2 %>% ggplot() + geom_sf() + facet_wrap(~time)

Int.coastal.3 %>% ggplot() + geom_sf() + facet_wrap(~time)



aoi <- read_sf('/Volumes/MaloneLab/Research/ENP/The Whitezone/WhiteZoneMap/SamplingGrid/AOI.kml')
aoi.line <-  st_cast(aoi,"LINESTRING")
aoi.point <- st_sample(aoi.line, size=30, type='regular') %>% st_cast( "POINT")