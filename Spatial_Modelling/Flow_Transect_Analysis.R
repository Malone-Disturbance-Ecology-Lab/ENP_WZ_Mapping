
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

convert_to_decimal <- function(degrees_minutes_seconds) {
  # Extract degrees, minutes, and seconds
  degrees <- as.numeric(substr(degrees_minutes_seconds, 1, 2))  # First 2 digits (degrees)
  minutes <- as.numeric(substr(degrees_minutes_seconds, 3, 4))  # Next 2 digits (minutes)
  seconds <- as.numeric(substr(degrees_minutes_seconds, 5, 10))  # Last 5 digits (seconds with decimals)
  
  # Convert to decimal degrees
  decimal_degrees <- degrees + (minutes / 60) + (seconds / 3600)
  
  return(decimal_degrees)
}

polygonize_try <- function(stack) {
  all_polygons <- list()
  
  for (i in 1:nlyr(stack)) {
    cat("Processing layer:", i, "\n")
    rast <- stack[[i]]
    names(rast) <- "value"
    
    # Keep only values 1, 2, 3; set others including NA to NA
    valid_vals <- c(1, 2, 3)
    mask_rast <- rast
    vals <- values(rast)
    vals[!vals %in% valid_vals] <- NA
    values(mask_rast) <- vals
    
    # Polygonize the masked raster (dissolve=FALSE to keep separate polygons for each patch)
    polys <- as.polygons(mask_rast, dissolve=FALSE)
    
    # Filter out polygons with NA value (if any)
    polys <- polys[!is.na(polys$value), ]
    
    # Add time attribute if available
    polys$time <- time(rast)
    
    all_polygons[[i]] <- polys
  }
  
  # Combine all polygons into one sf object
  result <- do.call(rbind, lapply(all_polygons, st_as_sf))
  return(result)
}

# color_scale is for leaflet plotting
path <- get_directory() ## get current directory path
setwd(path)
source( './function_MakeVegetationZoneLayers.R')
source( './function_Make_TransectLayers.R')
source( './function_for_plots.R')

rast.files <- list.files( '/Volumes/MaloneLab/Research/Whitezone/Result/Emulation/Final/CDthreshold06/CD_closest_GE', pattern = "tif")

setwd( '/Volumes/MaloneLab/Research/Whitezone/Result/Emulation/Final/CDthreshold06/CD_closest_GE')
#coast <- rast( rast.files) %>% terra::project("+proj=longlat +datum=WGS84") ## need to remove all the xml files
coast <- rast(rast.files) %>%
  terra::project("+proj=longlat +datum=WGS84", method = "near")
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

coast.polygons <- polygonize_try(coast.factor)
merged_polygons <- coast.polygons %>%
  group_by(time, value) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

coast.1.sf.new <- merged_polygons %>% 
  filter(value == 1) %>%  # keep only value == 1
  st_as_sf()              # ensure it's an sf object

coast.2.sf.new <- merged_polygons %>% 
  filter(value == 2) %>%  # keep only value == 1
  st_as_sf()              # ensure it's an sf object

coast.3.sf.new <- merged_polygons %>% 
  filter(value == 3) %>%  # keep only value == 1
  st_as_sf()              # ensure it's an sf object

Int.coastal.1.new <- Trans_Intersection(ENP.transects,coast.1.sf.new  )
Int.coastal.2.new <- Trans_Intersection(ENP.transects,coast.2.sf.new  )
Int.coastal.3.new <- Trans_Intersection(ENP.transects,coast.3.sf.new  )
#plot(merged_polygons[[2]]["value"], main = paste("Polygons for Year", merged_by_time[[1]]$time[1]))
st_write(merged_polygons[3, ], "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/merged_year_2_3.shp", delete_layer = TRUE)

# Write it to a shapefile
# st_write(first_polygon, "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/third_polygon.shp", delete_layer = TRUE)
# coast.1.first <- coast.1[[1]]
# coast.1.sub.sf <- polygonize_SpatRast(stack = coast.1.first)
# coast.2.first <- coast.2[[1]]
# coast.2.sub.sf <- polygonize_SpatRast(stack = coast.2.first)
# coast.3.first <- coast.3[[1]]
# coast.3.sub.sf <- polygonize_SpatRast(stack = coast.3.first)
# st_write(coast.1.sub.sf, "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/coast_1_shapefile.shp")
# st_write(coast.2.sub.sf, "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/coast_2_shapefile.shp")
# st_write(coast.3.sub.sf, "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/coast_3_shapefile.shp")

#Int.coastal.1 %>% ggplot() + geom_sf() + facet_wrap(~time)
#Int.coastal.2 %>% ggplot() + geom_sf() + facet_wrap(~time)
#Int.coastal.3 %>% ggplot() + geom_sf() + facet_wrap(~time)


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
########################## Output info for driver analysis ###############
# Add a source column to identify the dataset
df1 <- Int.coastal.1 %>% mutate(source = "Int.coastal.1")
df2 <- Int.coastal.2 %>% mutate(source = "Int.coastal.2")
df3 <- Int.coastal.3 %>% mutate(source = "Int.coastal.3")
# Bind all rows together
all_data <- bind_rows(df1, df2, df3)

# Group by ID and time, then calculate total length per group
summary_data <- all_data %>%
  group_by(ID, time) %>%
  summarise(
    total_length = sum(length.new, na.rm = TRUE),
    length_Int1 = sum(length.new[source == "Int.coastal.1"], na.rm = TRUE),
    length_Int2 = sum(length.new[source == "Int.coastal.2"], na.rm = TRUE),
    length_Int3 = sum(length.new[source == "Int.coastal.3"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    proportion_Int1 = length_Int1 / total_length,
    proportion_Int2 = length_Int2 / total_length,
    proportion_Int3 = length_Int3 / total_length
  )

# View result
head(summary_data)
summary_data_df <- summary_data %>% st_drop_geometry()
write.csv(summary_data_df, "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/Results/summary_data.csv", row.names = FALSE)
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
station_of_interest <- "S197"
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
correlation_results_dbh <- compute_dbhydro_correlations(
  all_site_names = all_site_names,
  dbhydro = dbhydro,
  summary_data = summary_data
)

correlation_results_dbh_df <- bind_rows(correlation_results_dbh) %>%
  distinct()

############# Plot correlation between WZ proportions and water level mean across transects ###
station_of_interest <- "NP-EPSW"

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

################### Experiment with one site only ############################
# first_file <- csv_files[4]
# cat("Reading file:", first_file, "\n")
# raw_data <- read_csv(first_file, col_names = c("time", "value"))
# unique(raw_data)
# print(head(raw_data))
# 
# clean_data <- raw_data %>%
#   slice(-1) %>%
#   mutate(
#     time = ymd(time),              # Use ymd() for "2016-01-25" format
#     value = as.numeric(value),
#     year = year(time)
#   )
# 
# annual_avg <- clean_data %>%
#   group_by(year) %>%
#   summarise(annual_avg = mean(value, na.rm = TRUE)) %>%
#   ungroup()


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
