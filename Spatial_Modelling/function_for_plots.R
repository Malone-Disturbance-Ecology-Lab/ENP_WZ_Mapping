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

################# get current directory #########
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


##### Generate buffer for transects ###
create_rect_buffer_safe <- function(transect, width_meters) {
  coords <- st_coordinates(transect)
  p1 <- coords[1, ]
  p2 <- coords[nrow(coords), ]
  
  dx <- p2["X"] - p1["X"]
  dy <- p2["Y"] - p1["Y"]
  length <- sqrt(dx^2 + dy^2)
  
  perp_x <- -dy / length
  perp_y <- dx / length
  
  hw <- width_meters / 2
  
  p1_left  <- c(p1["X"] + perp_x * hw, p1["Y"] + perp_y * hw)
  p1_right <- c(p1["X"] - perp_x * hw, p1["Y"] - perp_y * hw)
  p2_left  <- c(p2["X"] + perp_x * hw, p2["Y"] + perp_y * hw)
  p2_right <- c(p2["X"] - perp_x * hw, p2["Y"] - perp_y * hw)
  
  coords_mat <- rbind(p1_left, p2_left, p2_right, p1_right, p1_left)
  st_polygon(list(coords_mat))
}

#### Get yearly proportion of land cover ###
get_yearly_proportions <- function(transect_id, polygon, points, year_col) {
  inside_idx <- st_intersects(points, polygon, sparse = FALSE)[, 1]
  pts <- points[inside_idx, ]
  
  if (nrow(pts) == 0) return(NULL)
  
  vals <- st_drop_geometry(pts)[[year_col]]
  vals <- vals[!is.na(vals) & vals != 0]
  
  if (length(vals) == 0) {
    return(data.frame(
      transect_id = transect_id,
      year = year_col,
      prop_1 = NA,
      prop_2 = NA,
      prop_3 = NA
    ))
  }
  
  total <- length(vals)
  data.frame(
    transect_id = transect_id,
    year = year_col,
    prop_1 = sum(vals == 1) / total,
    prop_2 = sum(vals == 2) / total,
    prop_3 = sum(vals == 3) / total
  )
}

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

############### Salinity - read all the csv files and calculate annual avg ####
read_and_summarize_csv_files <- function(csv_files) {
  all_annual_averages <- lapply(csv_files, function(file) {
    cat("📂 Reading file:", file, "\n")
    
    # Read raw data with 2 columns: time and value
    raw_data <- read_csv(file, col_names = c("time", "value"), show_col_types = FALSE)
    
    # Remove header row if present
    data <- raw_data %>% slice(-1)
    
    # Try parsing dates and converting value
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
  
  return(all_annual_averages)
}

##### calculating correlations among sites ####
compute_site_correlations <- function(all_annual_averages, summary_data_filtered, site_names) {
  correlation_results <- lapply(site_names, function(site_name) {
    # Extract annual_avg for the site
    site_data <- all_annual_averages %>%
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
        correlation_result <- cor(
          transect_data$annual_avg,
          transect_data$proportion_Int2,
          use = "complete.obs"
        )
        return(data.frame(site = site_name, ID = ID_value, correlation = correlation_result))
      } else {
        return(NULL)  # Skip if there is insufficient data for correlation
      }
    })
    
    # Bind the list of correlations into a data frame for the site
    site_correlation_df <- bind_rows(site_correlation)
    return(site_correlation_df)
  })
  
  # Combine all sites into one data frame
  bind_rows(correlation_results)
}

### plot correlation maps
plot_correlation_map <- function(data, point_of_interest, ENP_transects, station_of_interest, color_scale) {
  leaflet(ENP_transects) %>%
    addProviderTiles("Esri.WorldImagery") %>%
    
    # Add station markers
    addCircleMarkers(
      data = data,
      lng = ~lon_decimal,
      lat = ~lat_decimal,
      color = "red",
      radius = 5,
      label = ~as.character(Station_ID),
      labelOptions = labelOptions(
        noHide = TRUE,
        direction = "top",
        textsize = "10px",
        fontweight = "bold"
      )
    ) %>%
    
    # Add station of interest as a gold star
    addCircleMarkers(
      data = point_of_interest,
      lng = ~lon_decimal,
      lat = ~lat_decimal,
      color = "gold",
      radius = 10,
      stroke = TRUE,
      weight = 3,
      fillOpacity = 1,
      label = paste0(station_of_interest),
      labelOptions = labelOptions(
        noHide = TRUE,
        direction = "top",
        textsize = "12px",
        fontweight = "bold"
      )
    ) %>%
    
    # Add transect lines colored by correlation
    addPolylines(
      data = ENP_transects,
      color = ~color_scale(correlation),
      weight = 3,
      opacity = 0.7,
      popup = ~paste0("Correlation: ", round(correlation, 3))
    ) %>%
    
    # Add legend
    addLegend(
      "bottomright",
      pal = color_scale,
      values = ENP_transects$correlation,
      title = "Correlation",
      opacity = 1
    )
}

#### calculate dbhydro wl correlations among sites and transects ##
compute_dbhydro_correlations <- function(all_site_names, dbhydro, summary_data) {
  total_sites <- length(all_site_names)
  
  correlation_results_dbh <- lapply(seq_along(all_site_names), function(i) {
    site_name <- all_site_names[[i]]
    
    # Progress tracking
    cat(sprintf("Progress: %d%% (%s)\n", round(100 * i / total_sites), site_name))
    
    # Extract mean_value for the site from dbhydro
    site_data <- dbhydro %>%
      filter(source_file == site_name) %>%
      distinct(year, .keep_all = TRUE)  # Keep only one row per year
    
    # Extract proportion_Int2 values for matching years
    site_summary <- summary_data %>%
      filter(time %in% site_data$year) %>%
      select(ID, time, proportion_Int2) %>%
      distinct(ID, time, .keep_all = TRUE)
    
    # Merge the site data with the transect data by year
    merged_data <- site_data %>%
      inner_join(site_summary, by = c("year" = "time"))
    
    # Calculate correlation for each ID
    site_correlation <- lapply(1:100, function(ID_value) {
      transect_data <- merged_data %>%
        filter(ID == ID_value)
      
      if (nrow(transect_data) > 1) {
        correlation_result <- cor(
          transect_data$mean_value,
          transect_data$proportion_Int2,
          use = "complete.obs"
        )
        return(data.frame(site = site_name, ID = ID_value, correlation = correlation_result))
      } else {
        return(NULL)
      }
    })
    
    # Combine correlations for the site
    bind_rows(site_correlation)
  })
  
  # Combine results across all sites
  bind_rows(correlation_results_dbh)
}