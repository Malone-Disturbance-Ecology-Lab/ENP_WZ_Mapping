##' @title This is the step 1 of the entire workflow 
##' @author Sparkle L. Malone, Boya (Paul) Zhang
##' @description Make buffer of transects and intersects with center points of pixels
##' @description Generating the changes of leading and trailing edges to calculate WZ expansion/retraction
##' 
##' 
rm(list = ls())
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
library(tools)
library(maps)
library(stringr)
library(leaflet)
library(ggspatial)

### Section 1 - Setting up and initial analysis
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

rast.files <- list.files(
  path = '/Volumes/MaloneLab/Research/Whitezone/Result/Emulation/Final/CDthreshold06/CD_closest_GE',
  pattern = "^raster.*\\.tif$",
  full.names = TRUE
)

#setwd( '/Volumes/MaloneLab/Research/Whitezone/Result/Emulation/Final/CDthreshold06/CD_closest_GE')
setwd( '/Volumes/MaloneLab/Research/Whitezone/Result/Emulation/Final/Roy')
coast <- rast(rast.files) %>%
  terra::project("+proj=longlat +datum=WGS84", method = "near") ## "near to make sure class value is not interpolated "
terra::time(coast) <-  substr(rast.files, 8, 11) %>% as.numeric # get the time out of the list of file names

ENP.transects <- makeTransect(line1 = coast.sf, line2 =  inland.sf, transect.n= 100)
# Calculate the centroid of each transect
ENP.transects <- ENP.transects %>%
  mutate(centroid = st_centroid(st_geometry(.)))

cent_coords <- st_coordinates(ENP.transects$centroid)

# Extract latitude and longitude from the centroid
ENP.transects <- ENP.transects %>%
  mutate(
    lon = cent_coords[, 1],
    lat = cent_coords[, 2]
  ) %>%
  select(-centroid)  # remove if centroid is not needed

########### Generating intersection between ENP.transects and study area ######
aoi <- st_read("/Volumes/MaloneLab/Research/Whitezone/Data/SEregion_reprojected_revise.shp")

aoi <- st_transform(aoi, st_crs(ENP.transects))
transects_in_polygon <- st_intersection(ENP.transects, aoi)
transects_in_polygon <- transects_in_polygon[, !(names(transects_in_polygon) == "id")]
st_write(transects_in_polygon, "/Volumes/MaloneLab/Research/Whitezone/Data/Intersected_Transects.shp", delete_layer = TRUE)

### Section 1 - calculating Intersections 
## Step 1 - get central points of the raster
coords <- as.data.frame(xyFromCell(coast, 1:ncell(coast)))
values <- as.data.frame(values(coast, na.rm = FALSE))
result_df <- cbind(coords, values)
colnames(result_df) <- c("lon", "lat", as.character(1984:2024))

## Step 2 - make buffer of transects
transects_utm <- st_transform(ENP.transects, crs = 32617)  # EPSG:32617 = UTM zone 17N
# Then generate your rectangle buffers using pixel size in meters
pixel_size_m <- 30  # Or whatever the pixel size is, in meters

buffer_polys_utm <- lapply(1:nrow(transects_utm), function(i) {
  create_rect_buffer_safe(transects_utm[i, ], width_meters = pixel_size_m)
})

buffer_sfc <- st_sfc(buffer_polys_utm, crs = 32617)

# Create sf object (add IDs if needed)
buffer_sf <- st_sf(transect_id = 1:length(buffer_polys_utm), geometry = buffer_sfc)
st_write(buffer_sf, "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/ENP_transect_rectangles.shp", delete_layer = TRUE)

buffer_polys_wgs84 <- lapply(buffer_polys_utm, function(p) st_transform(st_sfc(p, crs = 32617), crs = 4326)[[1]])

## Step 3 - Intersection the polygons with points
points_sf <- st_as_sf(result_df, coords = c("lon", "lat"), crs = 4326)
points_utm <- st_transform(points_sf, crs = 32617)
year_cols <- grep("^\\d{4}$", names(points_utm), value = TRUE)

# Step 4: Loop through transects and years
## This step takes 4.5 HOURS!!! Directly load the results.Rdata (below) instead
# results <- lapply(1:nrow(buffer_sf), function(i) {
#   poly <- buffer_sf$geometry[i]
#   transect_id <- buffer_sf$transect_id[i]
# 
#   do.call(rbind, lapply(year_cols, function(year_col) {
#     get_yearly_proportions(transect_id, poly, points_utm, year_col)
#   }))
# })
# save(results, file = "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/results.RData")
############# Output the NORTHERNEST WZ point for each transect ################
## it takes 4 hours!!
final_results <- extract_northernmost_points(buffer_sf, points_utm, year_cols)
south_results <- extract_northernmost_points(buffer_sf, points_utm, year_cols) ## 2 is for DZ

VZ_north_results <- extract_extreme_points(buffer_sf, points_sf, year_cols, direction = "north", 1) ## 1 is for VZ
VZ_south_results <- extract_extreme_points(buffer_sf, points_sf, year_cols, direction = "south", 1) ## 1 is for VZ

save(final_results, file = "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/north_point_Roy.RData")
save(south_results, file = "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/south_point_Roy.RData")

save(final_results, file = "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/north_point_Roy_VZ.RData")


save(final_results, file = "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/north_point_Roy.RData")
save(south_results, file = "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/south_point_Roy.RData")
load("/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/north_point.RData")
# Filter data for each year
points_1984 <- final_results %>% filter(Year == 1984)
points_1994 <- final_results %>% filter(Year == 1994)
points_2024 <- final_results %>% filter(Year == 2024)
points_1984_clean <- points_1984 %>% filter(!is.na(X) & !is.na(Y))
points_2024_clean <- points_2024 %>% filter(!is.na(X) & !is.na(Y))
points_1994_clean <- points_1994 %>% filter(!is.na(X) & !is.na(Y))

south_1984 <- south_results %>% filter(Year == 1984)
south_2024 <- south_results %>% filter(Year == 2024)
south_1984_clean <- south_1984 %>% filter(!is.na(X) & !is.na(Y))
south_2024_clean <- south_2024 %>% filter(!is.na(X) & !is.na(Y))

# Convert to sf object (assuming coordinates are in meters and CRS is known, e.g., EPSG:32617 for UTM zone 17N)
# Replace EPSG code with your actual CRS if different
sf_1984 <- st_as_sf(points_1984_clean, coords = c("X", "Y"), crs = 32617)
sf_2024 <- st_as_sf(points_2024_clean, coords = c("X", "Y"), crs = 32617)
sf_1994 <- st_as_sf(points_1994_clean, coords = c("X", "Y"), crs = 32617)

sf_1984_s <- st_as_sf(south_1984_clean, coords = c("X", "Y"), crs = 32617)
sf_2024_s <- st_as_sf(south_2024_clean, coords = c("X", "Y"), crs = 32617)

st_write(sf_1984, "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/points_1984.shp", delete_layer = TRUE)
st_write(sf_2024, "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/points_2024.shp", delete_layer = TRUE)

##### Output the lines connecting both point shapefiles ###
# (Assuming the column is named "Polygon_ID")
common_ids <- intersect(sf_1984$Polygon_ID, sf_2024$Polygon_ID)

############# signs of the connection line ###########
sf_1994$Polygon_ID <- as.integer(sf_1994$Polygon_ID)
sf_1984$Polygon_ID <- as.integer(sf_1984$Polygon_ID)
sf_2024$Polygon_ID <- as.integer(sf_2024$Polygon_ID)

sf_1984_s$Polygon_ID <- as.integer(sf_1984_s$Polygon_ID)
sf_2024_s$Polygon_ID <- as.integer(sf_2024_s$Polygon_ID)

line_sf <- calculate_signed_distance(
  sf_1984 = sf_1984,
  #sf_2024 = sf_2024,
  sf_2024 = sf_2024,
  output_dir = "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/",
  filename = "north_interior.shp"
)

line_sf_s <- calculate_signed_distance(
  sf_1984 = sf_1984_s,
  #sf_2024 = sf_2024,
  sf_2024 = sf_2024_s,
  output_dir = "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/",
  filename = "south_interior.shp"
)

line_sf$annual_rate_m <- line_sf$signed_distance_m / 40

################# Calculate the time series of interior boundary to coastline ######
### first intersect transect with polygon
southern <- get_southernmost_intersections(
  polygon_input = "/Users/bz294/Documents/WhiteZone/Points/Shapefiles/SEregion_reprojected_revise.shp",
  transects_input = ENP.transects,  # or a file path to transects shapefile
  id_field = "ID"  # change if your transect ID field has a different name
)

plot_avg_distance_to_coast(final_results, southern, output_path = "/Users/bz294/Documents/WhiteZone/Figure/avg_distance_to_coastline_2.tif")
