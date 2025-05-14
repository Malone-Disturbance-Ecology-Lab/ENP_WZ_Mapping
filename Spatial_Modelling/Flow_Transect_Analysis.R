
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
library(tools)
library(maps)
library(stringr)
library(leaflet)

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

# coast.1.first <- coast.1[[1]]
# coast.1.sub.sf <- polygonize_SpatRast(stack = coast.1.first)
# coast.2.first <- coast.2[[1]]
# coast.2.sub.sf <- polygonize_SpatRast(stack = coast.2.first)
# coast.3.first <- coast.3[[1]]
# coast.3.sub.sf <- polygonize_SpatRast(stack = coast.3.first)
# st_write(coast.1.sub.sf, "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/coast_1_shapefile.shp")
# st_write(coast.2.sub.sf, "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/coast_2_shapefile.shp")
# st_write(coast.3.sub.sf, "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/coast_3_shapefile.shp")

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
folder_path <- "/Users/bz294/Documents/WhiteZone/Data/Driver_Analysis/Salinity/WZ_salinity_data"

# List all CSV files in the folder
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# Loop through each CSV file and track failures
all_annual_averages <- lapply(csv_files, function(file) {
  cat("Reading file:", file, "\n")
  
  # Read raw data with 2 columns: time and value
  raw_data <- read_csv(file, col_names = c("time", "value"))
  
  # Remove header row if present
  data <- raw_data %>% slice(-1)
  
  # Try parsing dates
  data <- data %>%
    mutate(
      parsed_time = parse_date_time(time, orders = c("ymd", "mdy", "dmy")),
      value = as.numeric(value),
      site = file_path_sans_ext(basename(file))
    )
  
  # Check for parsing failures
  failed_rows <- data %>% filter(is.na(parsed_time))
  
  if (nrow(failed_rows) > 0) {
    cat("❌ Failed to parse dates in file:", file, "\n")
    print(failed_rows)
  }
  
  # If no valid dates, skip to next
  if (all(is.na(data$parsed_time))) {
    return(NULL)
  }
  
  # Calculate annual average
  annual_avg <- data %>%
    filter(!is.na(parsed_time)) %>%
    mutate(year = year(parsed_time)) %>%
    group_by(site, year) %>%
    summarise(annual_avg = mean(value, na.rm = TRUE), .groups = "drop")
  
  return(annual_avg)
}) %>%
  bind_rows()

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

# Filter all_annual_averages for the top 5 sites
top_5_data <- all_annual_averages %>%
  filter(site %in% top_5_site_names) %>%
  select(site, year, annual_avg)

# Step 2: Filter summary_data_df for the years that match
matching_years <- top_5_data$year  # Get the years from top_5_sites

summary_data_filtered <- summary_data_df %>%
  filter(time %in% matching_years)

correlation_results <- lapply(top_5_site_names, function(site_name) {
  # Extract annual_avg for the site
  site_data <- top_5_data %>%
    filter(site == site_name)
  
  # Extract the proportion_Int2 values for the matching years
  site_summary <- summary_data_filtered %>%
    filter(time %in% site_data$year) %>%
    select(ID, time, proportion_Int2)
  
  # Merge the site data (annual_avg) with the transect data (proportion_Int2) by year
  merged_data <- site_data %>%
    inner_join(site_summary, by = c("year" = "time"))
  
  # Initialize an empty list to store correlations for each transect
  site_correlation <- lapply(1:100, function(ID_value) {
    # Filter data for the given transect ID
    transect_data <- merged_data %>%
      filter(ID == ID_value)
    
    # Perform the correlation analysis for the current transect (ID)
    if (nrow(transect_data) > 1) {
      correlation_result <- cor(transect_data$annual_avg, transect_data$proportion_Int2, use = "complete.obs")
      return(data.frame(site = site_name, ID = ID_value, correlation = correlation_result))
    } else {
      return(NULL)  # Skip if there is insufficient data for correlation
    }
  })
  
  # Bind the list of correlations into a data frame for the site
  site_correlation_df <- bind_rows(site_correlation)
  return(site_correlation_df)
})

# Step 4: Combine the results into one data frame
correlation_results_df <- bind_rows(correlation_results)

# View the correlation results
print(correlation_results_df)

###################### Map the salinity stations and correlations on a map ####################
convert_to_decimal <- function(degrees_minutes_seconds) {
  # Extract degrees, minutes, and seconds
  degrees <- as.numeric(substr(degrees_minutes_seconds, 1, 2))  # First 2 digits (degrees)
  minutes <- as.numeric(substr(degrees_minutes_seconds, 3, 4))  # Next 2 digits (minutes)
  seconds <- as.numeric(substr(degrees_minutes_seconds, 5, 10))  # Last 5 digits (seconds with decimals)
  
  # Convert to decimal degrees
  decimal_degrees <- degrees + (minutes / 60) + (seconds / 3600)
  
  return(decimal_degrees)
}

csv_file <- "/Volumes/Malonelab/Research/Whitezone/Hydro_Salinity/Salinity/WZ_salinity_sites_latlon.csv"  # replace with your file path
data <- read.csv(csv_file)
data <- na.omit(data)  # Remove rows with N/A values

# Assuming the 2nd and 3rd columns are Lat and Lon in DDMMSS.SS format
data <- data %>%
  mutate(
    lat_decimal = sapply(data[[3]], convert_to_decimal),  # Convert lat to decimal degrees
    lon_decimal = sapply(data[[4]], convert_to_decimal)  # Convert lon to decimal degrees
  )
data$lon_decimal <- -abs(data$lon_decimal)

# Define the color scale for the correlation values
color_scale <- colorNumeric(palette = "YlGnBu", domain = ENP.transects_with_corr$correlation_value)


# Extract correlation values for site FLAB08_salinity
correlation_flab08 <- correlation_results_df %>%
  filter(site == "FLAB08_salinity") %>%
  select(ID, correlation)

# Ensure 'correlation' is numeric and remove NAs
correlation_flab08 <- correlation_flab08 %>%
  mutate(correlation = as.numeric(correlation)) %>%
  filter(!is.na(correlation))

ENP.transects_with_corr <- ENP.transects %>%
  left_join(correlation_flab08, by = "ID")

# Define the color scale based on correlation
color_scale <- colorNumeric(palette = "YlGnBu", domain = ENP.transects_with_corr$correlation)

flab08_point <- data %>%
  filter(Station_ID == "FLAB08_salinity")

# Plot map
leaflet(ENP.transects_with_corr) %>%
  addProviderTiles("Esri.WorldImagery") %>%
  
  # Add existing station markers
  addCircleMarkers(data = data,
                   lng = ~lon_decimal,
                   lat = ~lat_decimal,
                   color = "red",
                   radius = 5,
                   label = ~as.character(Station_ID),
                   labelOptions = labelOptions(noHide = TRUE, 
                                               direction = "top",
                                               textsize = "10px",
                                               fontweight = "bold")) %>%
  addCircleMarkers(data = flab08_point,
                   lng = ~lon_decimal,
                   lat = ~lat_decimal,
                   color = "gold",
                   radius = 10,
                   stroke = TRUE,
                   weight = 3,
                   fillOpacity = 1,
                   label = "FLAB08_salinity",
                   labelOptions = labelOptions(noHide = TRUE,
                                               direction = "top",
                                               textsize = "12px",
                                               fontweight = "bold"))
  # Add transect lines colored by correlation
  addPolylines(color = ~color_scale(correlation),
               weight = 3,
               opacity = 0.7) %>%
  
  # Add legend
  addLegend("bottomright",
            pal = color_scale,
            values = ~correlation,
            title = "Correlation",
            opacity = 1)


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
