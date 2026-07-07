##' @title This is the step2 of the workflow  
##' @author Boya (Paul) Zhang
##' @description first load the results of intersections with all the transects and proportion of land-cover class in time
##' @description Make plots, including total area changes in time (Figure 7)
##' @description Calculate transect-based correlation with salinity and wl, and make interactive plots
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