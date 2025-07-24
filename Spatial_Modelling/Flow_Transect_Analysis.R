
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

setwd( '/Volumes/MaloneLab/Research/Whitezone/Result/Emulation/Final/CDthreshold06/CD_closest_GE')
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
### output the results as shapefiles
## coast.1.first <- coast.1[[1]]
## coast.1.sub.sf <- polygonize_SpatRast(stack = coast.1.first)
## coast.2.first <- coast.2[[1]]
## coast.2.sub.sf <- polygonize_SpatRast(stack = coast.2.first)
## coast.3.first <- coast.3[[1]]
## coast.3.sub.sf <- polygonize_SpatRast(stack = coast.3.first)
## st_write(coast.1.sub.sf, "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/coast_1_shapefile.shp")
## st_write(coast.2.sub.sf, "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/coast_2_shapefile.shp")
## st_write(coast.3.sub.sf, "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/coast_3_shapefile.shp")

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
south_results <- extract_extreme_points(buffer_sf, points_sf, year_cols, direction = "south")
#save(final_results, file = "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/north_point.RData")
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
############## Experiment ###################################
## load the results, instead of running the four steps above. 
load("/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/results.RData")

# Step 5: Combine all into one data frame
summary_data_df <- do.call(rbind, results)
summary_data_df$year <- as.double(summary_data_df$year)
colnames(summary_data_df) <- c("ID", "time", "proportion_Int1", "proportion_Int2", "proportion_Int3")
write.csv(summary_data_df, "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/Results/summary_data_proportion.csv", row.names = FALSE)

#################### Output the transect layer,  intersection transect layer, and the trend plot  ###########
## Transect layer
output_shapefile <- "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/transect.tif"
#st_write(ENP.transects, "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/full_transect.shp")

## intersected transect layer
save_filtered_shapefile("/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/test_1984_wz.shp", Int.coastal.2, 1984)

## trend plot
output_file <- "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/TZ_trend_.tif"
#save_sf_plot_as_tiff(output_file, Int.coastal.3)
## plot the new trend based on raster, not shapefiles ###
selected_layers <- c(coast[[1]], coast[[21]], coast[[41]])
plot_trend_layers(selected_layers,
                  output_path = "/Users/bz294/Documents/WhiteZone/Figure/Concept_trend_raw.tif")

###################### Make Date Present Plot ##############################
set.seed(123)
# Example usage:
#output_path <- "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/Transect_Presence_TZ_.tif"
#save_transect_presence_plot(output_path, Int.coastal.3)

###################### Make Time Series Plot ##############################
#output_path <- "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/total_length_by_time_TZ_.tiff"
#save_total_length_plot(Int.coastal.3, output_path)
###################### Make Change Rate Plot ##############################
#time_sequence <- 1984:2024
#output_path <- "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/change_rate_std_plot_WZ_.tiff"
#save_change_rate_std_plot(Int.coastal.2, output_path, time_sequence)

####################### Area analysis for the raster stack ##########
plot_land_cover_time_series("/Volumes/Malonelab/Research/Whitezone/Result/Emulation/Final/CDthreshold06/CD_closest_GE", 
                                                         "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/land_cover_time_series_2.tif")

###################### Summarize the annual salinity and water level ###########
# Set the path to your folder
folder_path <- "/Volumes/Malonelab/Research/Whitezone/Hydro_Salinity/Salinity/WZ_salinity_data"

# List all CSV files in the folder
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# calcualte the annual avg of salinity values for each site
salinity_data <- read_and_summarize_csv_files(csv_files)

###################### Correlation analysis between transects and salinity ####
### start with long stations

# Step 1: Extract relevant data from salinity_data based on top_5_sites
all_site_names <- salinity_data$site
all_site_names <- all_site_names[!grepl("^C", all_site_names)] ## remove the sites starting from C
all_site_names <- all_site_names[!grepl("^S197", all_site_names)]

# Filter salinity_data for the top 5 sites
salinity_data_filtered <- salinity_data %>%
  filter(site %in% all_site_names) %>%
  select(site, year, annual_avg)

# Step 2: Filter summary_data_df for the years that match
matching_years <- salinity_data_filtered$year  # Get the years from top_5_sites

summary_data_filtered <- summary_data_df %>%
  filter(time %in% matching_years)

summary_data_filtered$time <- as.double(summary_data_filtered$time)
correlation_df <- compute_site_correlations(
  all_annual_averages = salinity_data,
  summary_data_filtered = summary_data_filtered,
  site_names = all_site_names
)

correlation_salinity_df <- bind_rows(correlation_df)
correlation_salinity_df <- correlation_salinity_df %>%
  distinct()

correlation_salinity_df <- correlation_salinity_df %>%
  mutate(Station_ID = str_remove(site, "_salinity"))

csv_file <- "/Volumes/Malonelab/Research/Whitezone/Hydro_Salinity/Salinity/WZ_salinity_sites_latlon.csv"  # replace with your file path
salinity_summary <- read.csv(csv_file)
salinity_summary <- na.omit(salinity_summary)  # Remove rows with N/A values

# Assuming the 2nd and 3rd columns are Lat and Lon in DDMMSS.SS format
salinity_summary <- salinity_summary %>% ## this data is the records for salinity and lat/lon
  mutate(
    lat_decimal = sapply(salinity_summary[[3]], convert_to_decimal),  # Convert lat to decimal degrees
    lon_decimal = sapply(salinity_summary[[4]], convert_to_decimal)  # Convert lon to decimal degrees
  )
salinity_summary$lon_decimal <- -abs(salinity_summary$lon_decimal)
out_dir <- "/Users/yourname/Documents/Shapefiles"
## output the salinity summary shapefile
export_flab_stations(salinity_summary, '/Volumes/MaloneLab/Research/Whitezone/Result/Spatial_Modelling/Results/')

correlation_salinity_df <- calculate_distance_join(correlation_salinity_df, salinity_summary, ENP.transects)
print(correlation_salinity_df)


###################### Map the salinity stations and correlations on a map ####################
####################### Plot it on interactive app using leaflet #############
# Define the color scale
station_of_interest <- "FLAB08"
## other sites are: Long stations: FLAB08, 14, 10, 11, 12, Short stations: C111W3, S197, C111JB, C111TC, C111MC, C111AC.
# Extract correlation values for site FLAB08_salinity
correlation_of_interest <- correlation_salinity_df %>%
  filter(site == paste0(station_of_interest, "_salinity")) %>%
  select(ID, correlation)

# Ensure 'correlation' is numeric and remove NAs
correlation_of_interest <- correlation_of_interest %>%
  mutate(correlation = as.numeric(correlation)) %>%
  filter(!is.na(correlation))

ENP.transects_with_corr <- transects_in_polygon %>%
  left_join(correlation_of_interest, by = "ID")

ENP.transects_with_corr_clean <- ENP.transects_with_corr[!is.na(ENP.transects_with_corr$correlation), ]

color_scale <- colorNumeric(palette = "YlGnBu", domain = ENP.transects_with_corr$correlation)
# Define the color scale based on correlation
salinity_summary <- salinity_summary %>% ungroup()
point_of_interest <- salinity_summary %>% filter(Station_ID == station_of_interest)

plot_correlation_map(
  data = salinity_summary,
  point_of_interest = point_of_interest,
  ENP_transects = ENP.transects_with_corr_clean,
  station_of_interest = station_of_interest,
  color_scale = color_scale
)

####################### Water level DBHydro correlation analysis #################################
dbhydro <- read.csv("/Volumes/Malonelab/Research/Whitezone/Hydro_Salinity/Hydro/DBhydro/annual_summary.csv", stringsAsFactors = FALSE)

# Step 1: Extract relevant data from salinity_data based on top_5_sites
all_site_names <- dbhydro$source_file
all_site_names <- all_site_names[!all_site_names %in% c("JBTS", "MBTS", "OYSTERCK")]
# Filter salinity_data for the top 5 sites
top_5_data <- dbhydro %>%
  select(year, source_file, mean_value, sd_value)

# Step 2: Filter summary_data_df for the years that match
matching_years <- top_5_data$year  # Get the years from top_5_sites

total_sites <- length(all_site_names)
summary_data_df$time <- as.double(summary_data_df$time)
correlation_dbh_df <- compute_dbhydro_correlations(
  all_site_names = all_site_names,
  dbhydro = dbhydro,
  summary_data = summary_data_df
)

correlation_dbh_df <- bind_rows(correlation_dbh_df) %>%
  distinct()

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

correlation_dbh_df <- correlation_dbh_df %>%
  rename(Station_ID = site)
correlation_dbh_df <- calculate_distance_join(correlation_dbh_df, db_summary, transects_in_polygon)
#save(correlation_results_dbh_df, file = "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/Results/correlation_results_dbh_df.Rdata")
#load("/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/Results/correlation_results_dbh_df.Rdata")

############# Plot correlation between WZ proportions and water level mean across transects ###
station_of_interest <- "NP-EP1"

# Extract correlation values for site 
correlation_of_interest <- correlation_dbh_df %>%
  filter(Station_ID == station_of_interest) %>%
  select(ID, correlation)

# Ensure 'correlation' is numeric and remove NAs
correlation_of_interest <- correlation_of_interest %>%
  mutate(correlation = as.numeric(correlation)) %>%
  filter(!is.na(correlation))

ENP.transects_with_corr <- transects_in_polygon %>%
  left_join(correlation_of_interest, by = "ID")

ENP.transects_with_corr_clean <- ENP.transects_with_corr[!is.na(ENP.transects_with_corr$correlation), ]

point_of_interest <- db_summary %>%
  filter(Name == station_of_interest)

plot_correlation_map(
  data = db_summary,
  point_of_interest = point_of_interest,
  ENP_transects = ENP.transects_with_corr_clean,
  station_of_interest = station_of_interest,
  color_scale = color_scale
)

########################## EDEN analysis ##############################
eden <- read.csv("/Volumes/Malonelab/Research/Whitezone/Hydro_Salinity/Hydro/EDEN/annual_summary.csv",
                 stringsAsFactors = FALSE,
                 fileEncoding = "latin1")

# Step 1: Extract relevant data from salinity_data based on top_5_sites
all_site_names <- eden$source_file

# Filter salinity_data for the top 5 sites
top_5_data <- eden %>%
  select(year, source_file, mean_value, sd_value)

# Step 2: Filter summary_data_df for the years that match
matching_years <- top_5_data$year  # Get the years from top_5_sites

total_sites <- length(all_site_names)

correlation_eden_df <- compute_dbhydro_correlations(
  all_site_names = all_site_names,
  dbhydro = eden,
  summary_data = summary_data_df
)

correlation_eden_df <- bind_rows(correlation_eden_df) %>%
  distinct()

station_of_interest <- "EPSW"

# Extract correlation values for site 
correlation_of_interest <- correlation_eden_df %>%
  filter(site == station_of_interest) %>%
  select(ID, correlation)

# Ensure 'correlation' is numeric and remove NAs
correlation_of_interest <- correlation_of_interest %>%
  mutate(correlation = as.numeric(correlation)) %>%
  filter(!is.na(correlation))

ENP.transects_with_corr <- transects_in_polygon %>%
  left_join(correlation_of_interest, by = "ID")

ENP.transects_with_corr_clean <- ENP.transects_with_corr[!is.na(ENP.transects_with_corr$correlation), ]
eden_summary <- eden_summary %>%
  rename(Name = WZ_sites)
point_of_interest <- eden_summary %>%
  filter(Name == station_of_interest)

csv_file <- "/Volumes/Malonelab/Research/Whitezone/Hydro_Salinity/Hydro/EDEN/WZ_hydro_latlon.csv"  # replace with your file path
eden_summary <- read.csv(csv_file)
eden_summary <- na.omit(eden_summary)  # Remove rows with N/A values
eden_summary <- eden_summary %>%
  filter(
    !startsWith(WZ_sites, "Stillwater") &
      !startsWith(WZ_sites, "CP") &
      !startsWith(WZ_sites, "NCL") &
      !startsWith(WZ_sites, "OL") &
      !startsWith(WZ_sites, "Trout") &
      !WZ_sites %in% c("G-3777", "EVER5A", "Joe_Bay_2E", "West_Highway_Creek", "EP1R", "McCormick_Creek", "Taylor_River_at_mouth")
  )
# Combine lat and lat_min into decimal degrees
eden_summary$lat_decimal <- eden_summary$lat + sapply(eden_summary$lat_min, convert_minsec)

# Combine lon and lon_min into decimal degrees (note: negative for West)
eden_summary$lon_decimal <- eden_summary$lon - sapply(eden_summary$lon_min, convert_minsec)

eden_summary <- eden_summary %>%
  mutate(Station_ID = WZ_sites)

plot_correlation_map(
  data = eden_summary,
  point_of_interest = point_of_interest,
  ENP_transects = ENP.transects_with_corr_clean,
  station_of_interest = station_of_interest,
  color_scale = color_scale
)

correlation_eden_df <- correlation_eden_df %>%
  rename(Station_ID = site)

correlation_eden_df <- calculate_distance_join(correlation_eden_df, eden_summary, transects_in_polygon)

################# Combine the db and eden hydro and output a shapefile #####
export_combined_hydro(
  db_summary = db_summary,
  eden_summary = eden_summary,
  out_path = "/Volumes/MaloneLab/Research/Whitezone/Result/Spatial_Modelling/Results/combined_hydro.shp"
)
####################### Distance analysis ###########################
# Make sure ID is treated as a factor (optional for coloring)
correlation_dbh_df <- correlation_dbh_df %>%
  filter(!Station_ID %in% c("JBTS", "MBTS", "OYSTERCK"))

correlation_salinity_df <- correlation_salinity_df %>%
  filter(!grepl("^C", site))

correlation_dbh_df <- correlation_dbh_df %>%
  mutate(ID = as.factor(ID))

correlation_salinity_df <- correlation_salinity_df %>%
  mutate(ID = as.factor(ID))
# Get all unique Station_IDs
station_ids <- unique(correlation_dbh_df$Station_ID)
station_salinity_ids <- unique(correlation_salinity_df$Station_ID)

tiff("/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/Results/correlation_vs_distance_by_station_salinity.tif", 
     width = 3000, height = 2000, res = 300)  # Adjust size/resolution as needed


# Generate the faceted scatter plot
ggplot(correlation_salinity_df, aes(x = distance_m, y = correlation)) +
  geom_point(alpha = 0.7, size = 1.5) +
  facet_wrap(~ Station_ID, scales = "free") +
  labs(
    title = "Correlation vs Distance (per Station)",
    x = "Distance (m)",
    y = "Correlation"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    panel.spacing = unit(1, "lines")
  )

# Close TIFF device to save file
dev.off()
##################### Transect based plot (not station based) #########
tiff("/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/Results/correlation_eden_by_transect.tif", width = 4000, height = 4000, res = 300)
p <- ggplot(correlation_eden_df, aes(x = distance_m, y = correlation)) +
  geom_point(alpha = 0.6, size = 1) +
  facet_wrap(~ ID, scales = "free", ncol = 10) +  # 10 columns layout
  labs(
    title = "Correlation vs Distance for Each Transect",
    x = "Distance (m)",
    y = "Correlation"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 7),
    plot.title = element_text(hjust = 0.5),
    panel.spacing = unit(0.3, "lines"),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 8)
  )

# Print plot (renders in RStudio plot pane or notebook)
print(p)
dev.off()

####################### Analyzing the trend of wl and salinity change time series ###
# Convert site to factor (optional but recommended for clear facets)
salinity_data$site <- as.factor(salinity_data$site)

dbhydro <- dbhydro %>%
  rename(site = source_file)

eden <- eden %>%
  rename(site = source_file)

dbhydro$site <- as.factor(dbhydro$site)
eden$site <- as.factor(eden$site)

eden_cm <- eden %>%
  mutate(
    mean_value = mean_value * 30.48,  # Convert to cm
    site = str_replace(site, "mouth", "Mouth")  # Capitalize 'mouth' to 'Mouth'
  )

eden_cm <- eden_cm %>%
  filter(site %in% eden_summary$WZ_sites)

dbhydro_cm <- dbhydro %>%
  mutate(mean_value = mean_value * 30.48)

salinity_data <- salinity_data %>%
  rename(mean_value = annual_avg)
# Estimate wrap width: adjust divisor to tune wrapping
plot_eden_faceted(eden_cm, wrap_width = 15)

dbhydro_cm_filtered <- dbhydro_cm %>%
  filter(!str_detect(site, "MBTS|JBTS|OYSTERCK"))
plot_dbhydro_annual_avg(dbhydro_cm_filtered)

salinity_data$site <- str_remove(salinity_data$site, "_salinity$")
salinity_data <- salinity_data %>% 
  rename(mean_value = annual_avg)
plot_salinity_faceted(salinity_data)
############### Trend analysis based on site measurements ###############
correlations <- eden %>%
  filter(!site %in% c("Joe_Bay_2E", "Upstream_Taylor")) %>%  # Exclude sites
  group_by(site) %>%
  summarise(
    correlation = cor(year, mean_value, use = "complete.obs"),
    .groups = "drop"
  )

ggplot(correlations, aes(x = reorder(site, correlation), y = correlation)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Correlation between Year and Mean Value by Site",
       x = "Site",
       y = "Correlation")

site_data_counts <- dbhydro %>%
  group_by(source_file) %>%
  summarise(
    count = sum(!is.na(mean_value)),
    .groups = "drop"
  )

print(site_data_counts)

###################### Modelling the WZ expansion mechanism! ############
model_result_salinity <- correlation_salinity_df %>%
  filter(!grepl("^C", Station_ID)) %>%  # Exclude Station_IDs starting with "C"
  group_by(ID) %>%
  summarise(
    mean_correlation = mean(correlation, na.rm = TRUE),
    std_correlation = sd(correlation, na.rm = TRUE),
    .groups = "drop"
  )

### combined correlation_dbh_df and correlation_eden_df
correlation_hydro_comb <- rbind(correlation_dbh_df, correlation_eden_df)

model_result_hydro <- correlation_hydro_comb %>%
  group_by(ID) %>%
  summarise(
    mean_correlation = mean(correlation, na.rm = TRUE),
    std_correlation = sd(correlation, na.rm = TRUE),
    .groups = "drop"
  )
## figure out whether WZ is expanding or not?
## for salinity using summary_data_filtered, for dbh using summary_data_df
trend_result <- summary_data_df %>% ## this figures out whether WZ is expanding
  group_by(ID) %>%
  summarise(
    slope = coef(lm(proportion_Int2 ~ time))[2],
    .groups = "drop"
  ) %>%
  mutate(WZ_trend_direction = case_when(
    slope > 0 ~ 1,
    slope < 0 ~ -1,
    TRUE ~ 0  # handles slope == 0
  )) %>%
  select(ID, WZ_trend_direction)

trend_result_VZ <- summary_data_df %>% ## this figures out whether WZ is expanding
  group_by(ID) %>%
  summarise(
    slope = coef(lm(proportion_Int1 ~ time))[2],
    .groups = "drop"
  ) %>%
  mutate(VZ_trend_direction = case_when(
    slope > 0 ~ 1,
    slope < 0 ~ -1,
    TRUE ~ 0  # handles slope == 0
  )) %>%
  select(ID, VZ_trend_direction)

trend_result_value <- summary_data_df %>%
  arrange(ID, time) %>%
  group_by(ID) %>%
  summarise(
    prop_Int1_1984 = first(proportion_Int1),
    prop_Int2_1984 = first(proportion_Int2),
    prop_Int3_1984 = first(proportion_Int3),
    trend_Int1 = last(proportion_Int1) - first(proportion_Int1),
    trend_Int2 = last(proportion_Int2) - first(proportion_Int2),
    trend_Int3 = last(proportion_Int3) - first(proportion_Int3),
    .groups = "drop"
  )

# Join the two data frames by ID
length_results <- transects_in_polygon %>%
  mutate(length.org = as.numeric(length.org)) %>%  # Drop units here
  left_join(trend_result_value, by = "ID") %>%
  mutate(
    length_change_Int1 = length.org * trend_Int1,
    length_change_Int2 = length.org * trend_Int2,
    length_change_Int3 = length.org * trend_Int3
  )

output_path <- "/Volumes/MaloneLab/Research/Whitezone/Result/Spatial_Modelling/transects_with_trend.shp"
st_write(length_results, output_path, delete_layer = TRUE)

############################ Model the liklihood of change in Int1 and Int2 ######
analyze_length_changes(length_results)

############################### The end ###########################################
# Ensure all ID columns are the same type (e.g., character)
trend_result <- trend_result %>%
  mutate(ID = as.character(ID))

trend_result_VZ <- trend_result_VZ %>%
  mutate(ID = as.character(ID))

model_result_salinity <- model_result_salinity %>%
  mutate(ID = as.character(ID))

model_result_hydro <- model_result_hydro %>%
  mutate(ID = as.character(ID))

# Join with custom suffixes
combined_result <- trend_result %>%
  left_join(trend_result_VZ, by = "ID") %>%
  left_join(model_result_salinity, by = "ID", suffix = c("", ".salinity")) %>%
  left_join(model_result_hydro, by = "ID", suffix = c("", ".dbh"))

combined_result <- combined_result %>%
  rename(
    mean_correlation.hydro = mean_correlation.dbh,
    std_correlation.hydro = std_correlation.dbh,
    mean_correlation.salinity = mean_correlation,
    std_correlation.salinity = std_correlation
  )

combined_result$ID <- as.numeric(as.character(combined_result$ID))
line_sf$Polygon_ID <- as.numeric(as.character(line_sf$Polygon_ID))
WZ_change_mechanism <- merge(combined_result, line_sf, by.x = "ID", by.y = "Polygon_ID")

WZ_change_mechanism <- WZ_change_mechanism %>%
  select(-starts_with("std"))

combined_result <- combined_result %>%
  left_join(line_sf %>% select(Polygon_ID, signed_distance_m),
            by = c("ID" = "Polygon_ID"))
####################### Checking the WZ mechanism variable #########################
plot_correlation_distance(WZ_change_mechanism, "/Volumes/MaloneLab/Research/Whitezone/Result/Spatial_Modelling/Results/")
# Gather for plotting correlation vs distance for salinity and hydro in one plot
df_long <- WZ_change_mechanism %>%
  select(mean_correlation.salinity, mean_correlation.hydro, distance_m) %>%
  pivot_longer(cols = c(mean_correlation.salinity, mean_correlation.hydro), 
               names_to = "correlation_type", values_to = "correlation_value")

p <- ggplot(df_long, aes(x = distance_m, y = correlation_value, color = correlation_type)) + 
  geom_point(alpha = 0.6) + 
  geom_smooth(method = "lm", se = TRUE) + 
  scale_color_manual(values = c("mean_correlation.salinity" = "blue", "mean_correlation.hydro" = "red"),
                     labels = c("Salinity", "Hydro")) + 
  theme_minimal(base_family = "Times") + 
  labs(x = "Distance (m)", y = "Mean Correlation", color = "Correlation Type",
       title = "Correlation of Salinity and Hydro with Distance") + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))

# Save as TIFF
ggsave("/Volumes/MaloneLab/Research/Whitezone/Result/Spatial_Modelling/Results/correlation_vs_distance_new.tif", plot = p, width = 8, height = 6, dpi = 300, units = "in", device = "tiff")

lm_salinity <- lm(mean_correlation.salinity ~ distance_m, data = WZ_change_mechanism)
lm_hydro <- lm(mean_correlation.hydro ~ distance_m, data = WZ_change_mechanism)

summary(lm_salinity)
summary(lm_hydro)

thrd <- 0
combined_result <- combined_result %>%
  mutate(mechanism = case_when(
    signed_distance_m > 0 & mean_correlation.salinity > 0 & mean_correlation.hydro > thrd ~ "SLR - WZ up",
    signed_distance_m > 0 & mean_correlation.salinity > 0 & mean_correlation.hydro < thrd ~ "Water Management - WZ up",
    signed_distance_m < 0 & mean_correlation.salinity > 0 & mean_correlation.hydro < thrd ~ "Water Management - WZ down",
    signed_distance_m < 0 & mean_correlation.salinity > 0 & mean_correlation.hydro > thrd ~ "SLR - WZ down",
    TRUE ~ "Unknown or NA"
  ))

################## Make plot for the mechanism map !!!! ####################
combined_result <- combined_result %>%
  mutate(ID = as.factor(ID))

ENP.transects <- ENP.transects %>%
  mutate(ID = as.factor(ID))

ENP.transects_with_corr <- ENP.transects %>%
  left_join(combined_result, by = "ID")

ENP.transects_with_corr_clean <- ENP.transects_with_corr[!is.na(ENP.transects_with_corr$WZ_trend_direction), ]
ENP.transects_with_corr_clean <- ENP.transects_with_corr_clean %>%
  mutate(ID = as.character(ID)) %>%
  left_join(
    line_sf %>%
      st_drop_geometry() %>%
      mutate(Polygon_ID = as.character(Polygon_ID)) %>%
      select(Polygon_ID, annual_rate_m),
    by = c("ID" = "Polygon_ID")
  )
ENP.transects_with_corr_clean %>%
  select(ID, annual_rate_m, geometry) %>%
  st_write("/Volumes/MaloneLab/Research/Whitezone/Result/Spatial_Modelling/Results/ENP_with_rate.shp", delete_dsn = TRUE)

plot_mechanism_map(
  ENP_transects = ENP.transects_with_corr_clean,
  output_dir = '/Volumes/MaloneLab/Research/Whitezone/Result/Spatial_Modelling/mechanism_0.shp'
)


ggplot(data = ENP.transects_with_corr_clean) +
  geom_sf(aes(color = mechanism), size = 1) +
  scale_color_manual(values = c('#d7191c','#fdae61','#abdda4','#2b83ba')) +
  theme_minimal(base_size = 15) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),       # remove grid lines
    axis.text = element_blank(),        # remove axis text
    axis.ticks = element_blank(),       # remove axis ticks
    axis.title = element_blank()        # remove axis titles
  ) +
  labs(color = "Mechanism")

ggsave("/Volumes/MaloneLab/Research/Whitezone/Result/Spatial_Modelling/mechanism_map.tif", width = 10, height = 8, dpi = 300)
############### Making plot for all the sites and transects - Data section ####
plot_station_sources(
  eden_summary = eden_summary,
  data = salinity_summary, 
  db_summary = db_summary, 
  ENP_transects = transects_in_polygon)

############### Plot WZ trend, signs of correlation with salinity and wl ###
plot_transect_signs(ENP.transects_with_corr_clean, "WZ_trend_direction", "WZ trend")

########### Output shapefiles for eden and dbhdyro #####
export_points_shapefile(
  df = eden_summary,
  lon_col = "lon_decimal",
  lat_col = "lat_decimal",
  output_path = "/Volumes/MaloneLab/Research/Whitezone/Result/Spatial_Modelling/Results/eden_summary_points.shp"
)
