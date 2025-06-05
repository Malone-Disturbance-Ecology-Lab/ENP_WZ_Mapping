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
library(RColorBrewer)
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

############# transects buffer analysis and counting proprotion of pixels classes ##
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
######## Plot the conceptual figure for the trend of all three classes##
plot_trend_layers <- function(selected_layers, years = c(1984, 2004, 2024), output_path = "trend_plot.tif") {
  # Helper to suppress sparse class-2 pixels
  suppress_sparse_twos <- function(v) {
    center <- v[5]
    if (is.na(center) || center != 2) return(center)
    count_2 <- sum(v == 2, na.rm = TRUE)
    if (count_2 < 4) return(NA)
    return(2)
  }

  # Process each layer
  processed_layers <- list()
  for (i in 1:3) {
    layer <- selected_layers[[i]]
    layer[layer == 0] <- NA
    if (i == 1) {
      processed_layers[[i]] <- focal(layer, w = 3, fun = suppress_sparse_twos, na.policy = "omit", fillvalue = NA)
    } else {
      processed_layers[[i]] <- focal(layer, w = 3, fun = mean, na.policy = "omit", na.rm = TRUE)
    }
  }

  # Open TIFF device
  tiff(filename = output_path, width = 12, height = 4, units = "in", res = 300)

  # Layout for 3 plots + 1 legend
  layout(matrix(c(1, 2, 3, 4), nrow = 1), widths = c(1.1, 1.1, 1.1, 0.4))
  par(mar = c(1, 1, 2, 1), family = "Times")

  # Colors and breaks
  cols <- c("#2b83ba", "#d7191c", "#ffffbf")
  brks <- c(0.5, 1.5, 2.5, 3.5)

  # Plot each subplot
  for (i in 1:3) {
    plot(processed_layers[[i]],
         main = as.character(years[i]),
         col = cols,
         breaks = brks,
         axes = FALSE,
         legend = FALSE,
         box = FALSE)
  }

  # Shared legend in last panel
  plot.new()
  legend("center",
         legend = c("VZ", "WZ", "TZ"),
         fill = cols,
         title = "Class",
         cex = 1.4,
         bty = "n",
         x.intersp = 1.2,
         y.intersp = 1.5)

  dev.off()
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

####### calcualte distance between sites and transects ####
calculate_distance_join <- function(correlation_salinity_df, data, ENP_transects) {
  
  # Step 1: Join lat/lon from data to correlation_salinity_df by Station_ID
  data_subset <- data %>% select(Station_ID, lat_decimal, lon_decimal)
  
  correlation_salinity_df <- correlation_salinity_df %>%
    left_join(data_subset, by = "Station_ID")
  
  # Step 2: Join lat/lon from ENP_transects to correlation_salinity_df by ID
  ENP_subset <- ENP_transects %>% select(ID, lat, lon)
  
  correlation_salinity_df <- correlation_salinity_df %>%
    left_join(ENP_subset, by = "ID")
  
  # Step 3: Rename columns for clarity
  correlation_salinity_df <- correlation_salinity_df %>%
    rename(
      site_lat = lat_decimal,
      site_lon = lon_decimal,
      tran_center_lat = lat,
      tran_center_lon = lon
    )
  
  # Step 4: Convert site and transect lat/lon to sf POINT objects
  # Step 4: Remove rows with missing coordinates
  correlation_salinity_df <- correlation_salinity_df %>%
    filter(!is.na(site_lat), !is.na(site_lon), !is.na(tran_center_lat), !is.na(tran_center_lon))
  
  # Then convert to sf POINT objects
  sites_sf <- correlation_salinity_df %>%
    st_as_sf(coords = c("site_lon", "site_lat"), crs = 4326, remove = FALSE)
  
  transects_sf <- correlation_salinity_df %>%
    st_as_sf(coords = c("tran_center_lon", "tran_center_lat"), crs = 4326, remove = FALSE)
  
  # Step 5: Calculate distances (meters) row-wise
  dist_m <- st_distance(sites_sf, transects_sf, by_element = TRUE)
  
  # Step 6: Add distance as a numeric column
  correlation_salinity_df <- correlation_salinity_df %>%
    mutate(distance_m = as.numeric(dist_m))
  
  return(correlation_salinity_df)
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

# Define a cleaner DMS-to-decimal converter
convert_minsec <- function(minsec_str) {
  parts <- regmatches(minsec_str, gregexpr("[0-9.]+", minsec_str))[[1]]
  minutes <- as.numeric(parts[1])
  seconds <- as.numeric(parts[2])
  return(minutes / 60 + seconds / 3600)
}

##### read eden files ###########
read_and_summarize_eden_files <- function(csv_files) {
  summaries <- lapply(csv_files, function(file) {
    # Read CSV, skip first 4 rows, so header is row 5
    df <- read_csv(file, skip = 4, col_names = TRUE)
    
    # Summarize numeric columns by mean (you can customize this)
    summary <- df %>%
      summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
      mutate(source_file = file)
    
    return(summary)
  })
  
  # Combine all summaries into one data frame
  all_summary <- bind_rows(summaries)
  return(all_summary)
}
##### plot correlation maps ##################################
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

######## plot the mechanism for the figure ###################
plot_mechanism_map <- function(ENP_transects_with_corr_clean) {
  
  # Create a color palette for mechanisms (categorical)
  mechanism_levels <- unique(ENP_transects_with_corr_clean$mechanism)
  color_scale <- colorFactor(
    palette = brewer.pal(min(length(mechanism_levels), 12), "Paired"),
    domain = mechanism_levels
  )
  
  leaflet(ENP_transects_with_corr_clean) %>%
    addProviderTiles("Esri.WorldImagery") %>%
    
    # Add transect lines colored by mechanism
    addPolylines(
      data = ENP_transects_with_corr_clean,
      color = ~color_scale(mechanism),
      weight = 3,
      opacity = 0.7,
      popup = ~paste0("Mechanism: ", mechanism)
    ) %>%
    
    # Add legend for mechanisms
    addLegend(
      "bottomright",
      pal = color_scale,
      values = ENP_transects_with_corr_clean$mechanism,
      title = "Mechanism",
      opacity = 1
    )
}

############## Making a plot with all the sites and transects - Data section ####
plot_station_sources <- function(eden_summary, data, db_summary, ENP_transects) {
  
  legend_df <- data.frame(
    label = c("EDEN_WL", "DBhydro_Salinity", "DBhydro_WL", "Transects"),
    color = c("deeppink", "dodgerblue", "limegreen", "orange")
  )
  
  leaflet() %>%
    addProviderTiles("Esri.WorldImagery") %>%
    
    # EDEN markers
    addCircleMarkers(
      data = eden_summary,
      lng = ~lon_decimal,
      lat = ~lat_decimal,
      color = "deeppink",
      radius = 5,
      label = ~paste0("EDEN: ", Station_ID)
    ) %>%
    
    # DBhydro Salinity markers
    addCircleMarkers(
      data = data,
      lng = ~lon_decimal,
      lat = ~lat_decimal,
      color = "dodgerblue",
      radius = 5,
      label = ~paste0("Salinity: ", Station_ID)
    ) %>%
    
    # DBhydro WL markers
    addCircleMarkers(
      data = db_summary,
      lng = ~lon_decimal,
      lat = ~lat_decimal,
      color = "limegreen",
      radius = 5,
      label = ~paste0("WL: ", Station_ID)
    ) %>%
    
    # Transect lines with transparency
    addPolylines(
      data = ENP_transects,
      color = "orange",
      weight = 2,
      opacity = 0.4
    ) %>%
    
    # Updated legend
    addLegend(
      position = "bottomright",
      colors = legend_df$color,
      labels = legend_df$label,
      title = "Data Source",
      opacity = 1
    )
}

################## Plot station data ######################
plot_eden_faceted <- function(data, wrap_width = 15, output_file = "/Volumes/MaloneLab/Research/Whitezone/Result/Spatial_Modelling/Results/eden_faceted_plot.tif") {
  p <- ggplot(data, aes(x = year, y = mean_value)) +
    geom_point(size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, color = "red", size = 1) +
    facet_wrap(
      ~ site,
      scales = "free_y",
      labeller = labeller(site = label_wrap_gen(width = wrap_width))
    ) +
    labs(
      x = "Year",
      y = "Annual Average"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold", size = 9),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 16),
      plot.title = element_blank()
    )
  
  # Save the plot
  tiff(output_file, width = 3000, height = 3000, res = 300)
  print(p)
  dev.off()
}

plot_dbhydro_annual_avg <- function(data, output_file = "/Volumes/MaloneLab/Research/Whitezone/Result/Spatial_Modelling/Results/annual_dbhydro_avg_by_site.tif") {
  p <- ggplot(data, aes(x = year, y = mean_value)) +
    geom_point(size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, color = "red", size = 1) +
    facet_wrap(~ site, scales = "free_y") +
    scale_x_continuous(limits = c(1984, 2024), breaks = seq(1985, 2025, by = 10)) +
    labs(
      x = "Year",
      y = "Water Level (cm)"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold", size = 10),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      plot.title = element_text(hjust = 0.5)
    )
  
  tiff(output_file, width = 3000, height = 3000, res = 300)
  print(p)
  dev.off()
}

plot_salinity_faceted <- function(data, wrap_width = 15, output_file = "/Volumes/MaloneLab/Research/Whitezone/Result/Spatial_Modelling/Results/salinity_faceted_plot.tif") {
  # Filter out sites starting with "C"
  data_filtered <- data[!grepl("^C", data$site), ]
  
  # Create the plot
  p <- ggplot(data_filtered, aes(x = year, y = mean_value)) +
    geom_point(size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, color = "red", size = 1) +
    facet_wrap(
      ~ site,
      scales = "free_y",
      labeller = labeller(site = label_wrap_gen(width = wrap_width))
    ) +
    labs(
      x = "Year",
      y = "Annual Average"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold", size = 9),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 16),
      plot.title = element_blank()
    )
  
  # Save the plot to a TIFF file
  tiff(output_file, width = 3000, height = 3000, res = 300)
  print(p)
  dev.off()
}

################ Plot transects signs #########
plot_transect_signs <- function(ENP.transects_with_corr_clean, column_name, legend_title) {
  
  column_sym <- sym(column_name)
  
  data_filtered <- ENP.transects_with_corr_clean %>%
    filter(!is.na(!!column_sym), !!column_sym != 0) %>%
    mutate(
      sign_value = ifelse(!!column_sym > 0, "Positive", "Negative"),
      sign_color = ifelse(sign_value == "Positive", "blue", "red")
    )
  
  leaflet(data_filtered) %>%
    addProviderTiles("Esri.WorldImagery") %>%
    addPolylines(
      color = ~sign_color,
      weight = 3,
      opacity = 0.8,
      popup = ~paste0(column_name, ": ", get(column_name))
    ) %>%
    addLegend(
      position = "bottomright",
      colors = c("blue", "red"),
      labels = c("Positive", "Negative"),
      title = legend_title,
      opacity = 1
    )
}

export_points_shapefile <- function(df, lon_col, lat_col, output_path) {
  # Convert to sf POINT object
  sf_obj <- st_as_sf(
    df,
    coords = c(lon_col, lat_col),
    crs = 4326,
    remove = FALSE
  )
  
  # Write to shapefile, overwrite if exists
  st_write(sf_obj, output_path, delete_layer = TRUE)
  
  message(paste("Shapefile saved to", output_path))
}
