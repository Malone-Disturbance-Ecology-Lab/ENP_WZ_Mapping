
##' @title This is the workflow to generate all the results 
##' @param line1 is the first lie
##' @param line2 is the second line
##' @param transect.n is the number of transects you want
##' @return sfc a simple feature of the transects and their lengths (length.org; Meters)
##' @author Sparkle L. Malone, Boya (Paul) Zhang
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

path <- get_directory() ## get current directory path
setwd(path)
source( './function_MakeVegetationZoneLayers.R')
source( './function_Make_TransectLayers.R')
source( './function_for_plots.R')

rast.files <- list.files( '/Volumes/MaloneLab/Research/Whitezone/Result/Emulation/Final/CDthreshold06/CD_closest_GE', pattern = "tif")

setwd( '/Volumes/MaloneLab/Research/Whitezone/Result/Emulation/Final/CDthreshold06/CD_closest_GE')
coast <- rast(rast.files) %>%
  terra::project("+proj=longlat +datum=WGS84", method = "near") ## "near to make sure class value is not interpolated "
terra::time(coast) <-  substr(rast.files, 8, 11) %>% as.numeric # get the time out of the list of file names

ENP.transects <- makeTransect(line1 = coast.sf, line2 =  inland.sf, transect.n= 100)

### output the results as shapefiles
# coast.1.first <- coast.1[[1]]
# coast.1.sub.sf <- polygonize_SpatRast(stack = coast.1.first)
# coast.2.first <- coast.2[[1]]
# coast.2.sub.sf <- polygonize_SpatRast(stack = coast.2.first)
# coast.3.first <- coast.3[[1]]
# coast.3.sub.sf <- polygonize_SpatRast(stack = coast.3.first)
# st_write(coast.1.sub.sf, "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/coast_1_shapefile.shp")
# st_write(coast.2.sub.sf, "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/coast_2_shapefile.shp")
# st_write(coast.3.sub.sf, "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/coast_3_shapefile.shp")

######## calculating Ints ########################
## Step 1 - get central points of the raster
coords <- as.data.frame(xyFromCell(coast, 1:ncell(coast)))
values <- as.data.frame(values(coast, na.rm = FALSE))
result_df <- cbind(coords, values)
colnames(result_df) <- c("lon", "lat", as.character(1984:2024))
## Step 2 - make buffer of transects
#pixel_size <- res(coast)[1]  # x-resolution (assume square pixels, so x = y)
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

buffer_polys_wgs84 <- lapply(buffer_polys, function(p) st_transform(st_sfc(p, crs = 32617), crs = 4326)[[1]])

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
#save(results, file = "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/results.RData")

load("/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/results.RData")

# Step 5: Combine all into one data frame
summary_data_df <- do.call(rbind, results)
summary_data_df$time <- as.double(summary_data_df$time)
colnames(summary_data_df) <- c("ID", "time", "proportion_Int1", "proportion_Int2", "proportion_Int3")
write.csv(summary_data_df, "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/Results/summary_data_proportion.csv", row.names = FALSE)

#################### Output the transect layer,  intersection transect layer, and the trend plot  ###########
### Transect layer
output_shapefile <- "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/transect.tif"
#st_write(ENP.transects, "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/full_transect.shp")

### intersected transect layer
save_filtered_shapefile("/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/test_1984_wz.shp", Int.coastal.2, 1984)

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

###################### Summarize the annual salinity and water level ###########
# Set the path to your folder
folder_path <- "/Volumes/Malonelab/Research/Whitezone/Hydro_Salinity/Salinity/WZ_salinity_data"

# List all CSV files in the folder
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# calcualte the annual avg of salinity values for each site
all_annual_averages <- read_and_summarize_csv_files(csv_files)

top_5_sites <- all_annual_averages %>%
  group_by(site) %>%
  summarise(n_years = n()) %>%
  arrange(desc(n_years)) %>%
  slice_head(n = 5)

print(top_5_sites)
###################### Correlation analysis between transects and salinity ####
### start with long stations

# Step 1: Extract relevant data from all_annual_averages based on top_5_sites
top_5_site_names <- top_5_sites$site  # Get the site names from top_5_sites
all_site_names <- all_annual_averages$site

# Filter all_annual_averages for the top 5 sites
top_5_data <- all_annual_averages %>%
  filter(site %in% top_5_site_names) %>%
  select(site, year, annual_avg)

# Step 2: Filter summary_data_df for the years that match
matching_years <- top_5_data$year  # Get the years from top_5_sites

summary_data_filtered <- summary_data_df %>%
  filter(time %in% matching_years)

summary_data_filtered$time <- as.double(summary_data_filtered$time)
correlation_df <- compute_site_correlations(
  all_annual_averages = all_annual_averages,
  summary_data_filtered = summary_data_filtered,
  site_names = all_site_names
)

correlation_results_df <- bind_rows(correlation_df)
# View the correlation results
print(correlation_results_df)

###################### Map the salinity stations and correlations on a map ####################
csv_file <- "/Volumes/Malonelab/Research/Whitezone/Hydro_Salinity/Salinity/WZ_salinity_sites_latlon.csv"  # replace with your file path
data <- read.csv(csv_file)
data <- na.omit(data)  # Remove rows with N/A values

# Assuming the 2nd and 3rd columns are Lat and Lon in DDMMSS.SS format
data <- data %>% ## this data is the records for salinity and lat/lon
  mutate(
    lat_decimal = sapply(data[[3]], convert_to_decimal),  # Convert lat to decimal degrees
    lon_decimal = sapply(data[[4]], convert_to_decimal)  # Convert lon to decimal degrees
  )
data$lon_decimal <- -abs(data$lon_decimal)

####################### Plot it on interactive app using leaflet #############
# Define the color scale
station_of_interest <- "FLAB08"
## other sites are: Long stations: FLAB08, 14, 10, 11, 12, Short stations: C111W3, S197, C111JB, C111TC, C111MC, C111AC.
# Extract correlation values for site FLAB08_salinity
correlation_of_interest <- correlation_results_df %>%
  filter(site == paste0(station_of_interest, "_salinity")) %>%
  select(ID, correlation)

# Ensure 'correlation' is numeric and remove NAs
correlation_of_interest <- correlation_of_interest %>%
  mutate(correlation = as.numeric(correlation)) %>%
  filter(!is.na(correlation))

ENP.transects_with_corr <- ENP.transects %>%
  left_join(correlation_of_interest, by = "ID")

ENP.transects_with_corr_clean <- ENP.transects_with_corr[!is.na(ENP.transects_with_corr$correlation), ]

color_scale <- colorNumeric(palette = "YlGnBu", domain = ENP.transects_with_corr$correlation)
# Define the color scale based on correlation

point_of_interest <- data %>%
  filter(Station_ID == station_of_interest)

plot_correlation_map(
  data = data,
  point_of_interest = point_of_interest,
  ENP_transects = ENP.transects_with_corr_clean,
  station_of_interest = station_of_interest,
  color_scale = color_scale
)
####################### Water level DBHydro correlation analysis #################################
dbhydro <- read.csv("/Volumes/Malonelab/Research/Whitezone/Hydro_Salinity/Hydro/DBhydro/annual_summary.csv", stringsAsFactors = FALSE)

# Step 1: Extract relevant data from all_annual_averages based on top_5_sites
all_site_names <- dbhydro$source_file

# Filter all_annual_averages for the top 5 sites
top_5_data <- dbhydro %>%
  select(year, source_file, mean_value, sd_value)

# Step 2: Filter summary_data_df for the years that match
matching_years <- top_5_data$year  # Get the years from top_5_sites

total_sites <- length(all_site_names)
summary_data_df$time <- as.double(summary_data_df$time)
correlation_results_dbh <- compute_dbhydro_correlations(
  all_site_names = all_site_names,
  dbhydro = dbhydro,
  summary_data = summary_data_df
)

#correlation_results_dbh_df <- bind_rows(correlation_results_dbh) %>%
#  distinct()

#save(correlation_results_dbh_df, file = "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/Results/correlation_results_dbh_df.Rdata")
load("/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/Results/correlation_results_dbh_df.Rdata")

############# Plot correlation between WZ proportions and water level mean across transects ###
station_of_interest <- "S197_T"

# Extract correlation values for site 
correlation_of_interest <- correlation_results_dbh_df %>%
  filter(site == station_of_interest) %>%
  select(ID, correlation)

# Ensure 'correlation' is numeric and remove NAs
correlation_of_interest <- correlation_of_interest %>%
  mutate(correlation = as.numeric(correlation)) %>%
  filter(!is.na(correlation))

ENP.transects_with_corr <- ENP.transects %>%
  left_join(correlation_of_interest, by = "ID")

ENP.transects_with_corr_clean <- ENP.transects_with_corr[!is.na(ENP.transects_with_corr$correlation), ]

csv_file <- "/Volumes/Malonelab/Research/Whitezone/Hydro_Salinity/Hydro/DBhydro/WZ_hydro_sites_latlon.csv"  # replace with your file path
db_summary <- read.csv(csv_file)
db_summary <- na.omit(db_summary)  # Remove rows with N/A values

# Assuming the 2nd and 3rd columns are Lat and Lon in DDMMSS.SS format
db_summary <- db_summary %>% ## this data is the records for salinity and lat/lon
  mutate(
    lat_decimal = sapply(db_summary[[2]], convert_to_decimal),  # Convert lat to decimal degrees
    lon_decimal = sapply(db_summary[[3]], convert_to_decimal)  # Convert lon to decimal degrees
  )
db_summary$lon_decimal <- -abs(db_summary$lon_decimal)
db_summary <- db_summary %>%
  mutate(Station_ID = Name)

point_of_interest <- db_summary %>%
  filter(Name == station_of_interest)

plot_correlation_map(
  data = db_summary,
  point_of_interest = point_of_interest,
  ENP_transects = ENP.transects_with_corr_clean,
  station_of_interest = station_of_interest,
  color_scale = color_scale
)

####################### Distance analysis ###########################

