
##' @title This is the workflow to generate all the results 
##' @param line1 is the first lie
##' @param line2 is the second line
##' @param transect.n is the number of transects you want
##' @return sfc a simple feature of the transects and their lengths (length.org; Meters)
##' @author Sparkle L. Malone, Boya (Paul) Zhang
##' 
library(terra)
library(tidyverse)
library(sf)
library(sfheaders)
library(units)
library(dplyr)
library(ggplot2)
library(reshape2)
library(here)
library(rstudioapi)

rm(list = ls())
### This function gets current directory
get_directory <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file <- "--file="
  rstudio <- "RStudio"
  
  match <- grep(rstudio, args)
  if (length(match) > 0) {
    return(dirname(rstudioapi::getSourceEditorContext()$path))
  } else {
    match <- grep(file, args)
    if (length(match) > 0) {
      return(dirname(normalizePath(sub(file, "", args[match]))))
    } else {
      return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
    }
  }
}

path <- get_directory() ## get current directory path
setwd(path)
source( './function_MakeVegetationZoneLayers.R')
source( './function_Make_TransectLayers.R')
source( './function_for_plots.R')

rast.files <- list.files( '/Volumes/MaloneLab/Research/Whitezone/Result/Emulation/Final/CDthreshold06/CD_closest_GE', pattern = "tif")

setwd( '/Volumes/MaloneLab/Research/Whitezone/Result/Emulation/Final/CDthreshold06/CD_closest_GE')
coast <- rast( rast.files) %>% terra::project("+proj=longlat +datum=WGS84") ## need to remove all the xml files
terra::time(coast) <-  substr(rast.files, 8, 11) %>% as.numeric # get the time out of the list of file names

coast.factor <- coast %>% as.factor
# make all zeros NA
coast.factor [ coast.factor  == 0] <- NA

# Need individual rast layers
coast.1 <- coast.factor
coast.1[ coast.1  != 1 ] <- NA
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

#Int.coastal.1 %>% ggplot() + geom_sf() + facet_wrap(~time)
#Int.coastal.2 %>% ggplot() + geom_sf() + facet_wrap(~time)
#Int.coastal.3 %>% ggplot() + geom_sf() + facet_wrap(~time)


#################### Output the transect layer,  intersection transect layer, and the trend plot  ###########
### Transect layer
output_shapefile <- "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/transect.tif"
st_write(ENP.transects, "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/full_transect.shp")

### intersected transect layer
save_filtered_shapefile("/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/test_2024.shp", Int.coastal.1, 2024)

### trend plot
output_file <- "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/TZ_trend_.tif"
save_sf_plot_as_tiff(output_file, Int.coastal.3)

###################### Make Date Present Plot ##############################
set.seed(123)
# Example usage:
output_path <- "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/Transect_Presence_TZ_.tif"
save_transect_presence_plot(output_path, Int.coastal.3)

###################### Make Time Series Plot ##############################
output_path <- "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/total_length_by_time_TZ_.tiff"
save_total_length_plot(Int.coastal.3, output_path)
###################### Make Change Rate Plot ##############################
time_sequence <- 1984:2024
output_path <- "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/change_rate_std_plot_WZ_.tiff"
save_change_rate_std_plot(Int.coastal.2, output_path, time_sequence)

####################### Area analysis for the raster stack ##########
plot_land_cover_time_series("/Volumes/Malonelab/Research/Whitezone/Result/Emulation/Final/CDthreshold06/CD_closest_GE", 
                                                         "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/land_cover_time_series_.tif")
################### Frequency plot for transects Line by Line ##############
# transect_1 <- ENP.transects %>% filter(ID == 70)
# intersect_1 <- Int.coastal.1 %>% filter(ID == 70)
# 
# intersect_1_sf <- st_as_sf(intersect_1, wkt = "x")  # Assuming 'x' is in WKT format
# print(st_geometry(intersect_1_sf))  # This should give a list of geometries
# 
# intersect_1_sf$freq <- 1
# 
# intersect_freq <- intersect_1_sf %>%
#   group_by(ID) %>%
#   summarise(geometry = st_union(st_geometry(.)), freq = sum(freq), .groups = "drop") %>%
#   st_as_sf()  # Ensure the result is an sf object with proper geometry

#aoi <- read_sf('/Volumes/MaloneLab/Research/ENP/The Whitezone/WhiteZoneMap/SamplingGrid/AOI.kml')
#aoi.line <-  st_cast(aoi,"LINESTRING")
#aoi.point <- st_sample(aoi.line, size=30, type='regular') %>% st_cast( "POINT")
