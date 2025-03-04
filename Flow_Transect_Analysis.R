
# This function is designed to turn vegetation zones into polygons.

library(terra)
library(tidyverse)
library(sf)
library(sfheaders)
library(units)
library(dplyr)
library(ggplot2)

rm(list = ls())
source( '/Users/bz294/Documents/WhiteZone/Code/function_MakeVegetationZoneLayers.R')

source( '/Users/bz294/Documents/WhiteZone/Code/function_Make_TransectLayers.R')

setwd( '/Volumes/MaloneLab/Research/Whitezone/Result/Emulation/Final/CDthreshold06/CD_closest_GE')
rast.files <- list.files( '/Volumes/MaloneLab/Research/Whitezone/Result/Emulation/Final/CDthreshold06/CD_closest_GE', pattern = "tif")

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
output_file <- "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/VZ_trend.tif"
Int.coastal.1 %>% 
  ggplot() + 
  geom_sf() + 
  facet_wrap(~time) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()
  )
ggsave(output_file, plot = last_plot(), device = "tiff", dpi = 300)
#################### Output the first layer as shapefile ###########
# Assuming Int.coastal.1 is an sf object
output_shapefile <- "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/test_layer.tif"
# Extract rows where time equals 1984
filtered_data <- Int.coastal.1[Int.coastal.1$time == 1984, ]

# Check if filtered_data is still a valid sf object
if (inherits(filtered_data, "sf")) {
  # Output the filtered data as a shapefile
  st_write(filtered_data, "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/test_1984.shp")
} else {
  print("Filtered data is not a valid sf object.")
}


Int.coastal.2 %>% 
  ggplot() + 
  geom_sf() + 
  facet_wrap(~time) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()
  )

output_file <- "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/WZ_trend.tif"
Int.coastal.2 %>% 
  ggplot() + 
  geom_sf() + 
  facet_wrap(~time) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()
  )
ggsave(output_file, plot = last_plot(), device = "tiff", dpi = 300)

Int.coastal.3 %>% 
  ggplot() + 
  geom_sf() + 
  facet_wrap(~time) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()
  )

################### Frequency plot for transects Line by Line ##############
transect_1 <- ENP.transects %>% filter(ID == 100)
intersect_1 <- Int.coastal.1 %>% filter(ID == 100)

intersect_1_sf <- st_as_sf(intersect_1, wkt = "x")  # Assuming 'x' is in WKT format
print(st_geometry(intersect_1_sf))  # This should give a list of geometries

intersect_1_sf$freq <- 1

intersect_freq <- intersect_1_sf %>%
  group_by(ID) %>%
  summarise(geometry = st_union(st_geometry(.)), freq = sum(freq), .groups = "drop") %>%
  st_as_sf()  # Ensure the result is an sf object with proper geometry

# all_segments <- Int.coastal.1 %>%
#   group_by(ID, time) %>%  # Group by transect ID and year
#   summarise(geometry = st_union(st_geometry(.)), .groups = "drop")




aoi <- read_sf('/Volumes/MaloneLab/Research/ENP/The Whitezone/WhiteZoneMap/SamplingGrid/AOI.kml')
aoi.line <-  st_cast(aoi,"LINESTRING")
aoi.point <- st_sample(aoi.line, size=30, type='regular') %>% st_cast( "POINT")
