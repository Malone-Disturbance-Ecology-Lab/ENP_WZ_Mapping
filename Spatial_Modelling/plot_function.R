# This function is designed to turn vegetation zones into polygons.
##' @author Boya Zhang, Sparkle L. Malone
##' 
library(terra)
library(tidyverse)
library(sf)
library(sfheaders)
library(units)
library(dplyr)
library(ggplot2)
library(reshape2)

### Plot trend map time series
save_sf_plot_as_tiff <- function(output_file, intersection_data) {
  # Create the ggplot with sf data and facet by time
  plot <- intersection_data %>% 
    ggplot() + 
    geom_sf() + 
    facet_wrap(~time) + 
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  # Save the plot as a TIFF file
  ggsave(output_file, plot = plot, device = "tiff", dpi = 300)
  
  message(paste("Plot saved as", output_file))
}

#### save shapefile
save_filtered_shapefile <- function(output_path, intersection_data, year) {
  # Extract rows where time equals the specified year
  filtered_data <- intersection_data[intersection_data$time == year, ]
  
  # Check if filtered_data is still a valid sf object
  if (inherits(filtered_data, "sf")) {
    # Output the filtered data as a shapefile
    st_write(filtered_data, output_path)
    message(paste("Shapefile saved to", output_path))
  } else {
    print("Filtered data is not a valid sf object.")
  }
}

#### Intersection presence
save_transect_presence_plot <- function(output_path, intersection_data) {
  # Create the plot
  p <- ggplot(intersection_data, aes(x = ID, y = time)) + 
    geom_point(shape = "-", size = 8) +  # Use horizontal bars to represent presence
    theme_minimal() + 
    labs(x = "Transect ID", y = "Year", title = "Transect TZ Presence by Year") + 
    scale_y_continuous(breaks = seq(1984, 2024, by = 5)) + 
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),  # White background
      plot.background = element_rect(fill = "white", color = NA),   # Ensure full white background
      plot.title = element_text(hjust = 0.5, face = "bold")  # Center title and make it bold
    )
  
  # Save the plot as a TIFF file
  ggsave(output_path, plot = p, width = 8, height = 6, dpi = 300, device = "tiff")
  
  message(paste("Plot saved as", output_path))
}

#### Total length of intersection time series
# Define the function
save_total_length_plot <- function(intersection_data, output_path) {
  # Aggregate the data by 'time' to get the total length
  aggregated_data <- intersection_data %>% 
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
  
  # Save the plot as a TIFF file
  ggsave(output_path, plot = p, device = "tiff", width = 8, height = 6, dpi = 300)
  
  message(paste("Plot saved as", output_path))
}

#### Change rate
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

# Define the function to calculate change rate and standard deviation, and plot it
save_change_rate_std_plot <- function(intersection_data, output_path, time_sequence) {
  # Calculate the change rate and standard deviation (assuming the function 'calculate_change_rate_and_std' exists)
  change_rate_std_df <- calculate_change_rate_and_std(intersection_data, time_sequence)
  
  # Add the upper and lower bounds for the shaded area
  change_rate_std_df <- change_rate_std_df %>%
    mutate(
      lower_bound = change_rate - std_dev / 2,
      upper_bound = change_rate + std_dev / 2
    )
  
  # Create the plot
  p <- ggplot(change_rate_std_df, aes(x = ID)) + 
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
      panel.grid.major = element_blank(),               # Remove major gridlines
      panel.grid.minor = element_blank(),               # Remove minor gridlines
      plot.title = element_text(hjust = 0.5)            # Center the plot title
    )
  
  # Save the plot as a TIFF file
  ggsave(output_path, plot = p, device = "tiff", width = 8, height = 6, dpi = 300)
  
  message(paste("Plot saved as", output_path))
}

## Area analysis
# Define the function
plot_land_cover_time_series <- function(raster_folder_path, output_path) {
  
  # Specify the folder containing your raster files
  raster_files <- list.files(path = raster_folder_path, pattern = "\\.tif$", full.names = TRUE)
  
  # Read all rasters into a SpatRaster stack
  raster_stack <- rast(raster_files)
  
  # Function to count pixels of each class
  count_classes <- function(r) {
    values <- as.vector(r)
    table(factor(values, levels = c(1, 2, 3)))  # Ensure all classes are counted
  }
  
  # Apply the function to each raster in the stack
  class_counts <- sapply(1:nlyr(raster_stack), function(i) count_classes(raster_stack[[i]]))
  
  # Convert to data frame
  time_series_df <- data.frame(
    Time = seq_along(raster_files),  # Time step index
    VZ = class_counts[1, ],
    WZ = class_counts[2, ],
    TZ = class_counts[3, ]
  )
  
  # Reshape for ggplot
  time_series_long <- melt(time_series_df, id.vars = "Time", variable.name = "Class", value.name = "PixelCount")
  
  time_series_long$Time <- time_series_long$Time - 1
  time_labels <- c("0" = "1984", "10" = "1994", "20" = "2004", "30" = "2014", "40" = "2024")
  
  # Create the plot
  p <- ggplot(time_series_long, aes(x = Time, y = PixelCount * 30 * 30 / 1000000, color = Class, group = Class)) +
    geom_line(size = 1) +
    geom_point() +
    labs(
      title = "Time Series of Land Cover Area",
      x = "Time Step",
      y = "Area (km²)",
      color = "Class"
    ) +
    theme_minimal(base_size = 14) +  # Set minimal theme
    scale_x_continuous(breaks = c(0, 10, 20, 30, 40), labels = time_labels) +  # Change X-axis labels
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),  # Center and bold title
      panel.background = element_rect(fill = "white", color = "white"),  # White background
      plot.background = element_rect(fill = "white", color = "white"), 
      legend.background = element_rect(fill = "white"), 
      legend.key = element_rect(fill = "white")
    )
  
  # Save the plot as a .tif file
  ggsave(output_path, plot = p, width = 8, height = 6, dpi = 300, device = "tiff")
  
  print(paste("Plot saved to:", output_path))
}

