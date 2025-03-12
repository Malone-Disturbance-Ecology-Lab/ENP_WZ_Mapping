
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

output_shapefile <- "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/transect.tif"
st_write(ENP.transects, "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/full_transect.shp")
######################## End of shapefile output ##################################
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

###################### Make Date Present Plot ##############################
set.seed(123)
p <- ggplot(Int.coastal.3, aes(x = ID, y = time)) +
  geom_point(shape = "-", size = 8) +  # Use horizontal bars to represent presence
  theme_minimal() +
  labs(x = "Transect ID", y = "Year", title = "Transect TZ Presence by Year") +
  scale_y_continuous(breaks = seq(1984, 2024, by = 5)) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),  # White background
    plot.background = element_rect(fill = "white", color = NA),  # Ensure full white background
    plot.title = element_text(hjust = 0.5, face = "bold")  # Center title and make it bold
  )
# Save as TIFF
ggsave("/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/Transect_Presence_TZ.tif", plot = p, width = 8, height = 6, dpi = 300, device = "tiff")
###################### Make Time Series Plot ##############################
aggregated_data <- Int.coastal.3 %>%
  group_by(time) %>%
  summarise(total_length_new = sum(length.new, na.rm = TRUE))

# Create the plot
p <- ggplot(aggregated_data, aes(x = time, y = total_length_new)) +
  geom_line(color = "blue", size = 1) +  # Line plot for summation of 'length.new'
  geom_point(color = "red", size = 2) +  # Points for each year
  labs(x = "Year", y = "Summation of Transects Length", title = "Total Length by Year for All TZ Transects") +
  theme(
    panel.background = element_rect(fill = "white"),  # Set the panel background to white
    plot.background = element_rect(fill = "white"),   # Set the plot background to white
    plot.title = element_text(hjust = 0.5)            # Center the plot title
  )

# Show the plot
ggsave("/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/total_length_by_time_TZ.tiff", plot = p, device = "tiff", width = 8, height = 6, dpi = 300)

###################### Make Change Rate Plot ##############################

time_sequence <- 1984:2024

# Function to calculate change rate and standard deviation for each transect ID
calculate_change_rate_and_std <- function(df, time_sequence) {
  # Create an empty data frame to store the results
  result <- data.frame(ID = integer(0), change_rate = numeric(0), std_dev = numeric(0))
  
  # Loop over each unique transect ID
  for (transect_id in unique(df$ID)) {
    # Filter the data for the specific transect ID
    transect_data <- df %>% filter(ID == transect_id)
    
    # Create the time series (match the time sequence 1984-2024)
    time_series <- rep(NA, length(time_sequence))
    for (i in 1:length(time_sequence)) {
      year_data <- transect_data %>% filter(time == time_sequence[i])
      if (nrow(year_data) > 0) {
        time_series[i] <- year_data$length.new  # Replace `value` with your actual column name for the data
      }
    }
    
    # If there are enough data points (>= 2), fit a linear model
    if (sum(!is.na(time_series)) >= 2) {
      model <- lm(time_series ~ time_sequence)
      
      # Extract the change rate (slope) and standard deviation (residuals)
      change_rate <- coef(model)[2]
      std_dev <- sd(resid(model))
      
      # Append the results for this transect ID
      result <- rbind(result, data.frame(ID = transect_id, change_rate = change_rate, std_dev = std_dev))
    }
  }
  
  return(result)
}

change_rate_std_df <- calculate_change_rate_and_std(Int.coastal.3, time_sequence)
# Add the upper and lower bounds for the shaded area
change_rate_std_df <- change_rate_std_df %>%
  mutate(
    lower_bound = change_rate - std_dev / 2,
    upper_bound = change_rate + std_dev / 2
  )

# Create the plot
p = ggplot(change_rate_std_df, aes(x = ID)) +
  # Shaded area (ribbon) for change rate with standard deviation bounds
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), fill = "blue", alpha = 0.3) +
  # Plot the change rate as a line
  geom_line(aes(y = change_rate), color = "blue", size = 1) +
  # Add points for the change rate values
  geom_point(aes(y = change_rate), color = "blue", size = 1) +
  # Plot a horizontal line at y = 0
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  # Labels and theme
  labs(x = "Transect ID", y = "Change Rate & Standard Deviation", title = "Change Rate with Standard Deviation for Each Transect") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),  # Set background to white
    plot.background = element_rect(fill = "white"),   # Set plot background to white
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )


ggsave("/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/change_rate_std_plot_TZ.tiff", plot = p, device = "tiff", width = 8, height = 6, dpi = 300)

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





aoi <- read_sf('/Volumes/MaloneLab/Research/ENP/The Whitezone/WhiteZoneMap/SamplingGrid/AOI.kml')
aoi.line <-  st_cast(aoi,"LINESTRING")
aoi.point <- st_sample(aoi.line, size=30, type='regular') %>% st_cast( "POINT")
