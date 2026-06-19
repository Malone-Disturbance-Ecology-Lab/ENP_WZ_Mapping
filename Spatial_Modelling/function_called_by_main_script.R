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

######################## Plot rates ####################################
plot_rates <- function(rate_results,
                       length_col = "length_change_Int2",  # column to use
                       utm_crs = 32617,
                       target_max_length = 180) {
  # 1. Project to UTM
  rates_proj <- st_transform(rate_results, utm_crs)
  
  # 2. Scale length based on chosen length column
  col_sym <- sym(length_col)
  max_rate <- max(pull(rates_proj, !!col_sym), na.rm = TRUE)
  
  rates_proj <- rates_proj %>%
    mutate(scaled_length = (!!col_sym) * target_max_length)
  
  # 3. Function to rescale a linestring
  rescale_line <- function(line_geom, new_length) {
    coords <- st_coordinates(line_geom)
    start_pt <- coords[1, 1:2]
    end_pt <- coords[nrow(coords), 1:2]
    direction <- end_pt - start_pt
    current_length <- sqrt(sum(direction^2))
    if (current_length == 0) return(line_geom)
    scaled_vector <- direction * (new_length / current_length)
    new_end <- start_pt + scaled_vector
    st_sfc(st_linestring(rbind(start_pt, new_end)), crs = st_crs(line_geom))
  }
  
  # 4. Rescale all lines
  rates_proj <- rates_proj %>%
    rowwise() %>%
    mutate(new_geom = rescale_line(x, scaled_length)) %>%
    ungroup() %>%
    st_set_geometry("new_geom")
  
  # 5. Compute reference points for connecting line
  rates_proj <- rates_proj %>%
    rowwise() %>%
    mutate(ref_point = {
      coords <- st_coordinates(new_geom)
      if (scaled_length >= 0) {
        # Positive: leading point (start of line)
        st_sfc(st_point(coords[1, 1:2]), crs = st_crs(new_geom))
      } else {
        # Negative: leading point (start of line)
        st_sfc(st_point(coords[1, 1:2]), crs = st_crs(new_geom))
      }
    }) %>%
    ungroup()
  
  # 6. Create connection line
  ref_coords <- st_coordinates(rates_proj$ref_point)
  connection_line <- st_sfc(st_linestring(ref_coords), crs = st_crs(rates_proj))
  
  # Return a list: rescaled transects + connecting line
  return(list(
    transects = rates_proj,
    connection_line = st_sf(geometry = connection_line)
  ))
}


###################### Plot transects ######################################
plot_and_save_transects <- function(results, connection_line, subset_transects_utm, filename = "transects_plot.tif", width = 6, height = 6, dpi = 300) {
  # Expecting results as a list: list(transects = sf, connection_line = sf)
  
  p <- ggplot() +
    geom_sf(data = results$x, color = "gray", size = 0.8) +
    geom_sf(data = connection_line, color = "black", size = 1) +
    geom_sf(data = subset_transects_utm, color = "blue", linewidth = 2, linetype = "dashed") +
    ggtitle("Rescaled Transects with Connection Line") +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  # Save as .tif
  ggsave(filename, plot = p, device = "tiff", width = width, height = height, dpi = dpi, bg = "white")
  
  message("Saved plot to: ", filename)
}

###################### obtaining the interior points for the time series #####
# Define main function
extract_northernmost_points <- function(buffer_sf, points_sf, year_cols) {
  
  # Nested helper function
  find_northernmost_point_value2 <- function(polygon, year_col, points_sf) {
    # Points inside polygon
    pts_in_poly <- points_sf[st_within(points_sf, polygon, sparse = FALSE), ]
    
    if (nrow(pts_in_poly) == 0) {
      return(data.frame(
        Year = year_col,
        Polygon_ID = polygon$ID,
        X = NA,
        Y = NA
      ))
    }
    
    # Filter points with land cover == 2 for that year
    pts_in_poly <- pts_in_poly %>% filter(.data[[year_col]] == 2)
    
    if (nrow(pts_in_poly) == 0) {
      return(data.frame(
        Year = year_col,
        Polygon_ID = polygon$ID,
        X = NA,
        Y = NA
      ))
    }
    
    coords <- st_coordinates(pts_in_poly)
    max_y_idx <- which.max(coords[, "Y"])
    
    data.frame(
      Year = year_col,
      Polygon_ID = polygon$ID,
      X = coords[max_y_idx, "X"],
      Y = coords[max_y_idx, "Y"]
    )
  }
  
  # Add ID to buffer if missing
  if (!"ID" %in% colnames(buffer_sf)) {
    buffer_sf$ID <- seq_len(nrow(buffer_sf))
  }
  
  # Initialize results
  results_list <- list()
  total_iters <- nrow(buffer_sf) * length(year_cols)
  pb <- txtProgressBar(min = 0, max = total_iters, style = 3)
  iter <- 0
  
  # Loop over polygons and years
  for (i in seq_len(nrow(buffer_sf))) {
    polygon <- buffer_sf[i, ]
    for (year in year_cols) {
      res <- find_northernmost_point_value2(polygon, year, points_sf)
      results_list[[length(results_list) + 1]] <- res
      
      iter <- iter + 1
      setTxtProgressBar(pb, iter)
    }
  }
  
  close(pb)
  
  # Combine and return
  final_results <- do.call(rbind, results_list)
  return(final_results)
}

########## the function provides options on extracting north or south points ####
extract_extreme_points <- function(buffer_sf, points_sf, year_cols, direction = "north") {
  
  # Nested helper function
  find_extreme_point <- function(polygon, year_col, points_sf, direction) {
    # Points inside polygon
    pts_in_poly <- points_sf[st_within(points_sf, polygon, sparse = FALSE), ]
    
    if (nrow(pts_in_poly) == 0) {
      return(data.frame(
        Year = year_col,
        Polygon_ID = polygon$ID,
        X = NA,
        Y = NA
      ))
    }
    
    # Filter points with land cover == 2 for that year
    pts_in_poly <- pts_in_poly %>% filter(.data[[year_col]] == 2)
    
    if (nrow(pts_in_poly) == 0) {
      return(data.frame(
        Year = year_col,
        Polygon_ID = polygon$ID,
        X = NA,
        Y = NA
      ))
    }
    
    coords <- st_coordinates(pts_in_poly)
    
    idx <- if (direction == "north") {
      which.max(coords[, "Y"])
    } else if (direction == "south") {
      which.min(coords[, "Y"])
    } else {
      stop("Invalid direction. Choose either 'north' or 'south'.")
    }
    
    data.frame(
      Year = year_col,
      Polygon_ID = polygon$ID,
      X = coords[idx, "X"],
      Y = coords[idx, "Y"]
    )
  }
  if (st_crs(points_sf) != st_crs(buffer_sf)) {
    points_sf <- st_transform(points_sf, st_crs(buffer_sf))
  }
  
  # Add ID to buffer if missing
  if (!"ID" %in% colnames(buffer_sf)) {
    buffer_sf$ID <- seq_len(nrow(buffer_sf))
  }
  
  # Initialize results
  results_list <- list()
  total_iters <- nrow(buffer_sf) * length(year_cols)
  pb <- txtProgressBar(min = 0, max = total_iters, style = 3)
  iter <- 0
  
  # Loop over polygons and years
  for (i in seq_len(nrow(buffer_sf))) {
    polygon <- buffer_sf[i, ]
    for (year in year_cols) {
      res <- find_extreme_point(polygon, year, points_sf, direction)
      results_list[[length(results_list) + 1]] <- res
      
      iter <- iter + 1
      setTxtProgressBar(pb, iter)
    }
  }
  
  close(pb)
  
  # Combine and return
  final_results <- do.call(rbind, results_list)
  return(final_results)
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
  cols <- c("#00441b", "#b08b5b", "#a1d996")
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


##### output the transect segment of interior shift ####
calculate_signed_distance <- function(sf_1984, sf_2024, output_dir, filename) {
  # Step 1: Get common IDs
  common_ids <- intersect(sf_1984$Polygon_ID, sf_2024$Polygon_ID)
  
  # Step 2: Filter both datasets
  sf_1984_filtered <- sf_1984 %>%
    filter(Polygon_ID %in% common_ids) %>%
    filter(!is.na(Polygon_ID))
  
  sf_2024_filtered <- sf_2024 %>%
    filter(Polygon_ID %in% common_ids) %>%
    mutate(
      X_2024 = st_coordinates(.)[, 1],
      Y_2024 = st_coordinates(.)[, 2]
    )
  print(common_ids)
  # Step 3: Merge with coordinates
  merged <- sf_1984_filtered %>%
    rename(geometry_1984 = geometry) %>%
    mutate(
      X_1984 = st_coordinates(geometry_1984)[, 1],
      Y_1984 = st_coordinates(geometry_1984)[, 2]
    ) %>%
    inner_join(
      sf_2024_filtered %>%
        st_drop_geometry() %>%
        select(Polygon_ID, X_2024, Y_2024),
      by = "Polygon_ID"
    ) %>%
    mutate(
      geometry_2024 = st_sfc(
        lapply(seq_len(nrow(.)), function(i) {
          st_point(c(.$X_2024[i], .$Y_2024[i]))
        }),
        crs = st_crs(sf_1984)
      ),
      north_flag = ifelse(Y_2024 > Y_1984, 1, -1)
    )
  
  # Step 4: Create lines
  lines <- st_sfc(
    lapply(1:nrow(merged), function(i) {
      pt1 <- c(merged$X_1984[i], merged$Y_1984[i])
      pt2 <- c(merged$X_2024[i], merged$Y_2024[i])
      st_linestring(rbind(pt1, pt2))
    }),
    crs = st_crs(merged$geometry_1984)
  )
  
  # Step 5: Construct sf with attributes
  line_sf <- st_sf(
    Polygon_ID = merged$Polygon_ID,
    north_flag = merged$north_flag,
    geometry = lines
  )
  
  # Step 6: Calculate distances
  line_sf$distance_m <- as.numeric(st_length(line_sf))
  line_sf$signed_distance_m <- line_sf$distance_m * line_sf$north_flag
  
  # Step 7: Summary stats
  signed_avg_distance <- mean(line_sf$signed_distance_m)
  max_expansion_distance <- max(line_sf$distance_m[line_sf$north_flag == 1])
  max_retraction_distance <- max(line_sf$distance_m[line_sf$north_flag == -1])
  
  cat("Signed average distance (m):", round(signed_avg_distance, 2), "\n")
  cat("Maximum expansion distance (m):", round(max_expansion_distance, 2), "\n")
  cat("Maximum retraction distance (m):", round(max_retraction_distance, 2), "\n")
  
  # Step 8: Write to shapefile
  out_path <- file.path(output_dir, filename)
  st_write(line_sf, out_path, delete_layer = TRUE)
  cat("Shapefile written to:", out_path, "\n")
  
  # Optionally return the sf object
  return(line_sf)
}

#### obtain southern point of the intersects between study area and transects 

get_southernmost_intersections <- function(polygon_input, transects_input, id_field = "ID") {
  # Load polygon if a file path
  if (is.character(polygon_input)) {
    polygon <- st_read(polygon_input)
  } else {
    polygon <- polygon_input
  }
  
  # Load transects if a file path
  if (is.character(transects_input)) {
    transects <- st_read(transects_input)
  } else {
    transects <- transects_input
  }
  
  # Ensure same CRS
  if (st_crs(transects) != st_crs(polygon)) {
    transects <- st_transform(transects, st_crs(polygon))
  }
  
  # Extract polygon boundary
  polygon_boundary <- st_boundary(polygon)
  
  # Intersect transects with polygon boundary
  intersections <- st_intersection(transects, polygon_boundary)
  
  # Cast intersection geometries to POINT (to get individual points)
  intersections_points <- st_cast(intersections, "POINT")
  
  # Extract coordinates and add as columns
  intersections_points <- intersections_points %>%
    mutate(
      X = st_coordinates(.)[,1],
      Y = st_coordinates(.)[,2]
    )
  
  # Select southernmost intersection per transect (minimum Y)
  southern_points <- intersections_points %>%
    group_by(across(all_of(id_field))) %>%
    slice_min(order_by = Y, n = 1) %>%
    ungroup()
  
  # Plot
  p <- ggplot() +
    geom_sf(data = polygon, fill = NA, color = "black", size = 0.8) +
    geom_sf(data = transects, color = "blue", size = 0.7) +
    geom_sf(data = southern_points, color = "red", size = 3) +
    theme_minimal() +
    labs(title = "Southernmost Intersection Points of Transects with Polygon Boundary",
         subtitle = "Polygon (black), Transects (blue), Southern Points (red)") +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12)
    )
  
  list(
    southern_points = southern_points,
    plot = p
  )
}

### Plot the distance to the coast plot ###
plot_avg_distance_to_coast <- function(final_results, southern, output_path = "avg_distance_to_coastline.tif") {
  
  # Prepare southern points data
  southern_df <- southern$southern_points %>%
    st_drop_geometry() %>%
    mutate(ID = as.character(ID)) %>%
    rename(X_south = X, Y_south = Y)
  
  # Prepare final_results
  final_results <- final_results %>%
    mutate(Polygon_ID = as.character(Polygon_ID))
  
  # Join southern points to final_results
  final_results_with_south <- final_results %>%
    left_join(southern_df, by = c("Polygon_ID" = "ID")) %>%
    mutate(
      distance = ifelse(
        is.na(X_south) | is.na(Y_south), NA,
        sqrt((X - X_south)^2 + (Y - Y_south)^2)
      )
    )
  
  # Compute average distance per year
  avg_dist_per_year <- final_results_with_south %>%
    group_by(Year) %>%
    summarise(avg_distance = mean(distance, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(Year = as.numeric(as.character(Year)))
  
  # Create the plot
  p <- ggplot(avg_dist_per_year, aes(x = Year, y = avg_distance)) +
    geom_line(color = "steelblue", size = 1) +
    geom_point(color = "darkred", size = 2) +
    theme_minimal(base_family = "Times New Roman", base_size = 14) +
    labs(
      x = "Year",
      y = "Avg Distance to Coast (m)"
    ) +
    scale_x_continuous(breaks = seq(1984, 2024, by = 5)) +
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 18),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  # Save the plot
  ggsave(
    filename = output_path,
    plot = p,
    width = 8,
    height = 6,
    dpi = 300,
    units = "in",
    device = "tiff",
    bg = "white"
  )
  
  # Return the summary table and plot object
  return(list(
    avg_dist_per_year = avg_dist_per_year,
    plot = p
  ))
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
plot_land_cover_time_series <- function(raster_folder_path, output_path, extend_years = NULL) {
  raster_files <- list.files(path = raster_folder_path, pattern = "\\.tif$", full.names = TRUE)
  raster_stack <- rast(raster_files)
  
  count_classes <- function(r) {
    values <- as.vector(r)
    table(factor(values, levels = c(1, 2, 3)))
  }
  
  class_counts <- sapply(1:nlyr(raster_stack), function(i) count_classes(raster_stack[[i]]))
  
  time_series_df <- data.frame(
    Time = seq_along(raster_files),
    VZ = class_counts[1, ],
    WZ = class_counts[2, ],
    TZ = class_counts[3, ]
  )
  
  time_series_df$WZ_area_km2 <- time_series_df$WZ * 30 * 30 / 1e6
  time_series_df$VZ_area_km2 <- time_series_df$VZ * 30 * 30 / 1e6
  time_series_df$TZ_area_km2 <- time_series_df$TZ * 30 * 30 / 1e6
  time_series_df$Year <- 1984:(1983 + length(raster_files))
  
  # Project future years if needed
  if (!is.null(extend_years)) {
    extend_to <- max(extend_years)
    last_observed_year <- max(time_series_df$Year[time_series_df$Year <= 2024])
    last_values <- time_series_df[time_series_df$Year == last_observed_year, ]
    
    # Annual change rates (km²)
    
    vz_rate <- 0.4922
    wz_rate <- 0.4589
    tz_rate <- -0.9511 ## annual rate of the predicted 25 and 50 yrs scenarios
    
    for (y in (last_observed_year + 1):extend_to) {
      prev_row <- tail(time_series_df, 1)
      new_row <- data.frame(
        Time = max(time_series_df$Time) + 1,
        VZ = NA, WZ = NA, TZ = NA,
        WZ_area_km2 = prev_row$WZ_area_km2 + wz_rate,
        VZ_area_km2 = prev_row$VZ_area_km2 + vz_rate,
        TZ_area_km2 = prev_row$TZ_area_km2 + tz_rate,
        Year = y
      )
      time_series_df <- rbind(time_series_df, new_row)
    }
  }
  
  # Reshape for plotting
  # ---- Reshape for plotting ----
  df_long <- melt(
    time_series_df[, c("Year", "WZ_area_km2", "VZ_area_km2", "TZ_area_km2")],
    id.vars = "Year",
    variable.name = "Class",
    value.name = "Area_km2"
  )
  df_long$Class <- gsub("_area_km2", "", df_long$Class)
  df_long$LineType <- ifelse(df_long$Year > 2024, "Projected", "Observed")
  
  # Make Class a factor to avoid level issues
  df_long$Class <- factor(df_long$Class, levels = c("VZ", "WZ", "TZ"))
  
  # ---- Base colors ----
  cols <- c("VZ" = "#00441b", "WZ" = "#b08b5b", "TZ" = "#a1d996")
  
  # ---- Remove NAs for plotting ----
  df_long_plot <- df_long[!is.na(df_long$Area_km2), ]
  
  p <- ggplot(df_long, aes(x = Year, y = Area_km2, color = Class, linetype = LineType)) +
    geom_line(size = 1.2) +
    geom_point(data = subset(df_long, LineType == "Observed"), size = 3) +
    geom_smooth(
      data = subset(df_long_plot,
                    Class %in% c("VZ") &
                      Year >= 2014 & Year <= 2024 &
                      LineType == "Observed"),
      aes(group = Class),
      method = "lm",
      se = FALSE,
      linewidth = 1.6,
      linetype = "solid",
      color = "grey"   # same color for both lines
    ) +
    scale_color_manual(values = cols) +
    scale_linetype_manual(values = c("Observed" = "solid", "Projected" = "dashed"), guide = "none") +
    scale_x_continuous(
      breaks = seq(1984, 2024, by = 10), 
      limits = c(1984, 2024)
    ) +
    #coord_cartesian(ylim = c(20, 100)) +
    #coord_cartesian(ylim = c(20, 110)) +
    coord_cartesian(ylim = c(0, 45)) +
    labs(
      title = "Predictions of Land Cover Area",
      x = "Year",
      y = "Area (km²)",
      color = "Class"
    ) +
    theme_minimal(base_size = 18, base_family = "Times New Roman") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", color = "black", size = 28),
      axis.title = element_text(color = "black", size = 28),
      axis.text = element_text(color = "black", size = 28),
      legend.title = element_text(color = "black", size = 28),
      legend.text = element_text(color = "black", size = 28),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white")
    )
  
  ggsave(output_path, plot = p, width = 8, height = 6, dpi = 300, device = "tiff")
  
  print(paste("Plot saved to:", output_path))
  
  # Print change rates
  get_rate <- function(start_year, end_year, df, col) {
    s <- df[df$Year == start_year, col]
    e <- df[df$Year == end_year, col]
    rate <- (e - s) / (end_year - start_year)
    return(rate)
  }
  
  cat("VZ Change Rate (km²/year):\n")
  cat("1984–2021:", round(get_rate(1984, 2021, time_series_df, "VZ_area_km2"), 2), "\n")
  cat("2000–2021:", round(get_rate(2000, 2021, time_series_df, "VZ_area_km2"), 2), "\n")
  cat("2010–2021:", round(get_rate(2010, 2021, time_series_df, "VZ_area_km2"), 2), "\n")
  cat("2014–2024:", round(get_rate(2014, 2024, time_series_df, "VZ_area_km2"), 2), "\n")
  
  # ---- Calculate total area by year ----
  time_series_df$Total_area_km2 <- time_series_df$VZ_area_km2 + 
    time_series_df$WZ_area_km2 + 
    time_series_df$TZ_area_km2
  
  cat("\nTotal Land Cover Area (km²) by Year:\n")
  print(time_series_df[, c("Year", "Total_area_km2")])
  
  # Example: print total area for a specific year (e.g., 2024)
  target_year <- 2022
  total_area_target <- time_series_df$Total_area_km2[time_series_df$Year == target_year]
  cat(paste("\nTotal area in", target_year, ":", round(total_area_target, 2), "km²\n"))
  
  # ---- Calculate total change in area for each class ----
  calculate_change <- function(df, start_year, end_year) {
    start_row <- df[df$Year == start_year, ]
    end_row <- df[df$Year == end_year, ]
    
    change_VZ <- end_row$VZ_area_km2 - start_row$VZ_area_km2
    change_WZ <- end_row$WZ_area_km2 - start_row$WZ_area_km2
    change_TZ <- end_row$TZ_area_km2 - start_row$TZ_area_km2
    
    cat("\nChange in area (km²) from", start_year, "to", end_year, ":\n")
    cat("VZ:", round(change_VZ, 2), "km²\n")
    cat("WZ:", round(change_WZ, 2), "km²\n")
    cat("TZ:", round(change_TZ, 2), "km²\n")
  }
  
  # Example usage: change from first to last observed year
  observed_start <- 2014
  observed_end <- 2024
  
  calculate_change(time_series_df, observed_start, observed_end)
  
  # ---- Print 2024 class areas ----
  target_year <- 2024
  year_row <- time_series_df[time_series_df$Year == target_year, ]
  
  cat("\nLand Cover Areas in", target_year, ":\n")
  cat("VZ:", round(year_row$VZ_area_km2, 2), "km²\n")
  cat("WZ:", round(year_row$WZ_area_km2, 2), "km²\n")
  cat("TZ:", round(year_row$TZ_area_km2, 2), "km²\n")
}
###################### Extended plot ################################################
# BAU
# Class 1 (VZ): Rate: 0.37 pixels/year)
# Class 2 (WZ):  Rate: 0.33 pixels/year)
# Class 3 (TZ): Rate: -0.70 pixels/year)

### pessimistic: Class 1 (VZ): Rate: 0.16 pixels/year)
### Class 2 (WZ): Rate: 0.50 pixels/year)
### Class 3 (TZ): Rate: -0.65 pixels/year)

## Optimistic
# Class 1 (VZ): Rate: 0.76 pixels/year)
# Class 2 (WZ): Rate: 0.06 pixels/year)
# Class 3 (TZ): Rate: -0.82 pixels/year)




plot_land_cover_time_series_extend <- function(raster_folder_path, output_prefix, extend_years = NULL) {

  raster_files <- list.files(path = raster_folder_path, pattern = "\\.tif$", full.names = TRUE)
  raster_stack <- rast(raster_files)
  
  count_classes <- function(r) {
    values <- as.vector(r)
    table(factor(values, levels = c(1, 2, 3)))
  }
  
  class_counts <- sapply(1:nlyr(raster_stack), function(i) count_classes(raster_stack[[i]]))
  
  time_series_df <- data.frame(
    Time = seq_along(raster_files),
    VZ = class_counts[1, ],
    WZ = class_counts[2, ],
    TZ = class_counts[3, ]
  )
  
  time_series_df$WZ_area_km2 <- time_series_df$WZ * 30 * 30 / 1e6
  time_series_df$VZ_area_km2 <- time_series_df$VZ * 30 * 30 / 1e6
  time_series_df$TZ_area_km2 <- time_series_df$TZ * 30 * 30 / 1e6
  time_series_df$Year <- 1984:(1983 + length(raster_files))
  
  # Define scenarios with their rates
  scenarios <- list(
    BAU = list(vz_rate = 0.37, wz_rate = 0.33, tz_rate = -0.70),
    Pessimistic = list(vz_rate = 0.16, wz_rate = 0.50, tz_rate = -0.65),
    Optimistic = list(vz_rate = 0.76, wz_rate = 0.06, tz_rate = -0.82)
  )
  
  # Loop through each scenario
  for (sc in names(scenarios)) {
    ts_df <- time_series_df
    
    if (!is.null(extend_years)) {
      extend_to <- max(extend_years)
      last_observed_year <- max(ts_df$Year[ts_df$Year <= 2024])
      last_values <- ts_df[ts_df$Year == last_observed_year, ]
      
      vz_rate <- scenarios[[sc]]$vz_rate
      wz_rate <- scenarios[[sc]]$wz_rate
      tz_rate <- scenarios[[sc]]$tz_rate
      
      for (y in (last_observed_year + 1):extend_to) {
        prev_row <- tail(ts_df, 1)
        new_row <- data.frame(
          Time = max(ts_df$Time) + 1,
          VZ = NA, WZ = NA, TZ = NA,
          WZ_area_km2 = prev_row$WZ_area_km2 + wz_rate,
          VZ_area_km2 = prev_row$VZ_area_km2 + vz_rate,
          TZ_area_km2 = prev_row$TZ_area_km2 + tz_rate,
          Year = y
        )
        ts_df <- rbind(ts_df, new_row)
      }
    }
    
    # Reshape for plotting
    df_long <- melt(ts_df[, c("Year", "WZ_area_km2", "VZ_area_km2", "TZ_area_km2")],
                    id.vars = "Year",
                    variable.name = "Class", value.name = "Area_km2")
    df_long$Class <- gsub("_area_km2", "", df_long$Class)
    df_long$LineType <- ifelse(df_long$Year > 2024, "Projected", "Observed")
    
    cols <- c("VZ" = "#00441b", "WZ" = "#b08b5b", "TZ" = "#a1d996")
    
    p <- ggplot(df_long, aes(x = Year, y = Area_km2, color = Class, linetype = LineType)) +
      geom_line(size = 1.2) +
      geom_point(data = subset(df_long, LineType == "Observed"), size = 3) +
      scale_color_manual(values = cols) +
      scale_linetype_manual(values = c("Observed" = "solid", "Projected" = "dashed"), guide = "none") +
      labs(
        x = "Year",
        y = "Area (km²)",
        color = "Class"
      ) +
      theme_minimal(base_size = 24, base_family = "Times New Roman") +
      theme(
        axis.title = element_text(color = "black", size = 24),
        axis.text = element_text(color = "black", size = 24),
        legend.title = element_text(color = "black", size = 24),
        legend.text = element_text(color = "black", size = 24),
        panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white", color = "white"),
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white")
      )
    
    # Save output for this scenario
    output_path <- paste0(output_prefix, "_", sc, ".tiff")
    ggsave(output_path, plot = p, width = 8, height = 6, dpi = 300, device = "tiff")
    
    print(paste("Plot saved for", sc, "scenario at:", output_path))
    
    # ---- Print final VZ area ----
    final_year <- max(ts_df$Year)
    final_vz_area <- ts_df$TZ_area_km2[ts_df$Year == final_year]
    cat(paste0("Scenario: ", sc, " | Final TZ area (km²) in ", final_year, ": ", round(final_vz_area, 3), "\n"))
  }
  
  
}
################ Plot interior boundary shift ###################################
plot_boundary_and_basemap <- function(shapefile_path, output_tif_path, transect_shapefile_path = NULL) {
  # 1. Read shapefile
  shp <- st_read(shapefile_path)
  
  # 2. Transform to Web Mercator
  shp_3857 <- st_transform(shp, 3857)
  
  # 3. Create a class column for legend
  shp_3857$Class <- ifelse(shp_3857$sgnd_d_ > 0, "Landward", "Seaward")
  
  # 4. Assign colors to classes
  line_colors <- c("Landward" = "red", "Seaward" = "cyan")
  
  # 5. Get basemap
  basemap <- get_tiles(shp_3857, provider = "Esri.WorldImagery", crop = TRUE)
  
  # 6. Base plot
  p <- ggplot() +
    layer_spatial(basemap) +
    geom_sf(data = shp_3857, aes(color = Class), size = 2) +
    scale_color_manual(values = line_colors) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 12, family = "Times New Roman")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3)))
  
  # 7. Overlay transects if provided
  if (!is.null(transect_shapefile_path)) {
    transects <- st_read(transect_shapefile_path)
    transects_3857 <- st_transform(transects, 3857)
    
    p <- p +
      geom_sf(
        data = transects_3857,
        color = "blue",       # or choose another color
        size = 1,
        linetype = "dashed"
      )
  }
  
  # 8. Display
  print(p)
  
  # 9. Save to file with white background
  ggsave(
    output_tif_path,
    plot = p,
    width = 8,
    height = 6,
    dpi = 300,
    bg = "white"
  )
  
  message("Plot saved to: ", output_tif_path)
  
  # 10. Summary table
  summary_df <- shp_3857 %>%
    mutate(sgnd_div40 = sgnd_d_ / 40) %>%
    mutate(range_group = case_when(
      Plyg_ID < 37 ~ "<37",
      Plyg_ID >= 37 & Plyg_ID <= 61 ~ "37-61",
      Plyg_ID > 61 ~ ">61"
    )) %>%
    group_by(range_group) %>%
    summarise(
      mean_sgnd_div40 = mean(sgnd_div40, na.rm = TRUE),
      median_sgnd_div40 = median(sgnd_div40, na.rm = TRUE),
      count = n()
    )
  
  print(summary_df)
}

#################### Report average rates for each subregion #####################
calculate_length_change_averages <- function(df, value_col, custom_range = c(27, 37)) {
  
  # ensure ID is numeric
  df <- df %>% mutate(ID = as.numeric(ID))
  
  # define main ID groups
  df <- df %>%
    mutate(ID_group = case_when(
      ID >= 1 & ID <= 37   ~ "1-37",
      ID >= 38 & ID <= 60  ~ "38-60",
      ID >= 61 & ID <= 100 ~ "61-100",
      TRUE ~ NA_character_
    ))
  
  # calculate averages for main groups
  avg_main <- df %>%
    group_by(ID_group) %>%
    summarise(avg_length_change = mean(.data[[value_col]], na.rm = TRUE)) %>%
    ungroup()
  
  # calculate average for custom range
  avg_custom <- df %>%
    filter(ID >= custom_range[1] & ID <= custom_range[2]) %>%
    summarise(ID_group = paste0(custom_range[1], "-", custom_range[2]),
              avg_length_change = mean(.data[[value_col]], na.rm = TRUE))
  
  # combine results
  result <- bind_rows(avg_main, avg_custom)
  
  return(result)
}

################# Calculate the average rate for WZ (41 years) and VZ (11 years) ##
calculate_avg_weighted_diff <- function(option = 1,
                                        year1 = 2013,
                                        year2 = 2024,
                                        layer = "proportion_Int1",
                                        summary_data_df,
                                        transects_in_polygon) {
  # Validate layer input
  if (!layer %in% c("proportion_Int1", "proportion_Int2")) {
    stop("Invalid layer. Choose either 'proportion_Int1' or 'proportion_Int2'.")
  }
  
  # Filter summary_data_df based on region option
  summary_filtered <- switch(
    as.character(option),
    "1" = summary_data_df[summary_data_df$ID >= 61, ],
    "2" = summary_data_df[summary_data_df$ID >= 39 & summary_data_df$ID < 61, ],
    "3" = summary_data_df[summary_data_df$ID < 39, ],
    stop("Invalid option. Choose 1, 2, or 3.")
  )
  
  # Subset for both years
  prop_year1 <- subset(summary_filtered, time == year1, select = c("ID", layer))
  prop_year2 <- subset(summary_filtered, time == year2, select = c("ID", layer))
  
  # Rename columns for merging
  names(prop_year1)[2] <- "prop_year1"
  names(prop_year2)[2] <- "prop_year2"
  
  # Merge by ID
  prop_diff_df <- merge(prop_year1, prop_year2, by = "ID")
  
  # Calculate difference
  prop_diff_df$prop_diff <- prop_diff_df$prop_year2 - prop_diff_df$prop_year1
  
  # Merge with transect lengths
  prop_diff_df <- merge(prop_diff_df, transects_in_polygon, by = "ID")
  
  # Multiply by length
  prop_diff_df$weighted_diff <- prop_diff_df$prop_diff * prop_diff_df$length.org
  
  # Return average weighted difference divided by number of years
  return(mean(prop_diff_df$weighted_diff, na.rm = TRUE) / (year2 - year1))
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

######### export salinity shapefile ##########
export_flab_stations <- function(data, out_dir) {
  # Step 1: Filter Station_IDs that start with "FLAB" and exclude "FLAB16"
  salinity_summary_filtered <- data %>%
    filter(grepl("^FLAB", Station_ID) & Station_ID != "FLAB16")
  
  # Step 2: Convert to sf object using lat/lon
  salinity_sf <- st_as_sf(
    salinity_summary_filtered,
    coords = c("lon_decimal", "lat_decimal"),
    crs = 4326  # WGS 84
  )
  
  # Step 3: Write to shapefile in the specified output directory
  shapefile_path <- file.path(out_dir, "FLAB_stations.shp")
  st_write(salinity_sf, shapefile_path, delete_dsn = TRUE)
  
  message("Shapefile saved to: ", shapefile_path)
}

############# export hydro shapefile ############################
export_combined_hydro <- function(db_summary, eden_summary, out_path) {
  # Step 1: Select relevant columns
  db_clean <- db_summary %>%
    select(Station_ID, lat_decimal, lon_decimal)
  
  eden_clean <- eden_summary %>%
    select(Station_ID, lat_decimal, lon_decimal)
  
  # Step 2: Combine the two dataframes
  combined_summary <- bind_rows(db_clean, eden_clean)
  
  # Step 3: Remove rows with missing coordinates
  combined_summary <- combined_summary %>%
    filter(!is.na(lat_decimal) & !is.na(lon_decimal))
  
  # Step 4: Convert to sf object
  combined_sf <- st_as_sf(
    combined_summary,
    coords = c("lon_decimal", "lat_decimal"),
    crs = 4326  # WGS 84
  )
  
  # Step 5: Write to shapefile
  st_write(combined_sf, out_path, delete_dsn = TRUE)
  
  message("Shapefile written to: ", out_path)
}
##### calculating correlations among sites by region ####
compute_site_correlations_by_region <- function(all_annual_averages,
                                                summary_data,
                                                west_sites,
                                                central_sites,
                                                east_sites,
                                                proportion_col) {
  
  # Ensure summary_data is a tibble
  summary_data <- as_tibble(summary_data)
  
  correlation_results <- lapply(1:100, function(ID_value) {
    
    # Determine region
    if (ID_value > 61) {
      region <- "west"
      site_names <- west_sites
    } else if (ID_value >= 37 && ID_value <= 61) {
      region <- "central"
      site_names <- central_sites
    } else {
      region <- "east"
      site_names <- east_sites
    }
    
    # Year range
    year_filter <- 1991:2024
    
    # Regional mean
    regional_mean <- all_annual_averages %>%
      dplyr::filter(site %in% site_names, year %in% year_filter) %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(mean_salinity = mean(annual_avg, na.rm = TRUE), .groups = "drop")
    
    # Transect data
    transect_data <- summary_data %>%
      dplyr::filter(ID == ID_value, time %in% regional_mean$year) %>%
      dplyr::select(ID, time, value = dplyr::all_of(proportion_col))
    
    # Merge and compute correlation
    merged <- regional_mean %>%
      dplyr::inner_join(transect_data, by = c("year" = "time"))
    
    if (nrow(merged) > 1) {
      correlation_result <- cor(merged$mean_salinity,
                                merged$value,
                                use = "complete.obs")
      return(data.frame(ID = ID_value,
                        region = region,
                        correlation = correlation_result))
    } else {
      return(NULL)
    }
  })
  
  dplyr::bind_rows(correlation_results)
}

####### calcualte distance between sites and transects ####
calculate_distance_join <- function(correlation_salinity_df, data, ENP_transects) {
  
  correlation_salinity_df <- correlation_salinity_df %>%
    mutate(ID = as.integer(ID))
  
  ENP_transects <- ENP_transects %>%
    mutate(ID = as.integer(ID))
  # Step 1: Join lat/lon from data to correlation_salinity_df by Station_ID
  data_subset <- data %>% select(Station_ID, lat_decimal, lon_decimal)
  
  correlation_salinity_df <- correlation_salinity_df %>%
    left_join(data_subset, by = "Station_ID")
  
  # Step 2: Join lat/lon from ENP_transects to correlation_salinity_df by ID
  ENP_transects <- ENP_transects %>%
    st_centroid() %>%
    mutate(
      lon = st_coordinates(.)[,1],
      lat = st_coordinates(.)[,2]
    )
  
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
compute_dbhydro_correlations <- function(all_site_names, wl, summary_data) {
  total_sites <- length(all_site_names)
  
  correlation_results_dbh <- lapply(seq_along(all_site_names), function(i) {
    site_name <- all_site_names[[i]]
    
    # Progress tracking
    cat(sprintf("Progress: %d%% (%s)\n", round(100 * i / total_sites), site_name))
    
    # Extract mean_value for the site from dbhydro
    site_data <- wl %>%
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
############ Model the length change versus initial condition #########
analyze_length_changes <- function(length_results) {
  # Fit linear models
  model_Int1 <- lm(length_change_Int1 ~ prop_Int1_1984 + prop_Int2_1984 + prop_Int3_1984, data = length_results)
  model_Int2 <- lm(length_change_Int2 ~ prop_Int1_1984 + prop_Int2_1984 + prop_Int3_1984, data = length_results)
  
  # Summarize the models
  summary_Int1 <- summary(model_Int1)
  summary_Int2 <- summary(model_Int2)
  
  # Print summaries
  cat("=== Linear Model: Length Change in Int1 ===\n")
  print(summary_Int1)
  
  cat("\n=== Linear Model: Length Change in Int2 ===\n")
  print(summary_Int2)
  
  # Basic interpretation logic
  cat("\n=== Interpretation ===\n")
  
  coef1 <- coef(model_Int1)
  coef2 <- coef(model_Int2)
  
  if (!is.na(coef1["prop_Int1_1984"]) && coef1["prop_Int1_1984"] > 0) {
    cat("Higher initial proportion of Int1 in 1984 tends to increase length change in Int1.\n")
  } else if (!is.na(coef1["prop_Int1_1984"])) {
    cat("Higher initial proportion of Int1 in 1984 tends to decrease length change in Int1.\n")
  }
  
  if (!is.na(coef1["prop_Int2_1984"]) && coef1["prop_Int2_1984"] > 0) {
    cat("Higher initial proportion of Int2 in 1984 tends to increase length change in Int1.\n")
  } else if (!is.na(coef1["prop_Int2_1984"])) {
    cat("Higher initial proportion of Int2 in 1984 tends to decrease length change in Int1.\n")
  }
  
  if (!is.na(coef2["prop_Int1_1984"]) && coef2["prop_Int1_1984"] > 0) {
    cat("Higher initial proportion of Int1 in 1984 tends to increase length change in Int2.\n")
  } else if (!is.na(coef2["prop_Int1_1984"])) {
    cat("Higher initial proportion of Int1 in 1984 tends to decrease length change in Int2.\n")
  }
  
  if (!is.na(coef2["prop_Int2_1984"]) && coef2["prop_Int2_1984"] > 0) {
    cat("Higher initial proportion of Int2 in 1984 tends to increase length change in Int2.\n")
  } else if (!is.na(coef2["prop_Int2_1984"])) {
    cat("Higher initial proportion of Int2 in 1984 tends to decrease length change in Int2.\n")
  }
  
  if (!is.na(coef2["prop_Int3_1984"]) && coef2["prop_Int3_1984"] > 0) {
    cat("Higher initial proportion of Int3 in 1984 tends to increase length change in Int2.\n")
  } else if (!is.na(coef2["prop_Int3_1984"])) {
    cat("Higher initial proportion of Int3 in 1984 tends to decrease length change in Int2.\n")
  }
}


######## plot the mechanism for the figure ###################

plot_mechanism_map <- function(ENP_transects_with_corr_clean, output_dir, options = leafletOptions(zoomControl = FALSE)) {
  # Use 'x' as geometry if needed; here assuming original geometry is used
  
  # Select only ID, mechanism, and geometry columns
  sf_for_export <- ENP_transects_with_corr_clean %>%
    select(ID, mechanism, geometry)
  
  # Export shapefile
  out_path <- file.path(output_dir)
  st_write(sf_for_export, out_path, delete_layer = TRUE)
  cat("Shapefile written to:", out_path, "\n")
  
  # Prepare color palette
  mechanism_levels <- unique(ENP_transects_with_corr_clean$mechanism)
  custom_colors <- c('#d7191c','#fdae61','#abdda4','#2b83ba')[1:length(mechanism_levels)]
  
  color_scale <- colorFactor(
    palette = custom_colors,
    domain = mechanism_levels
  )
  
  # Build leaflet map
  leaflet(ENP_transects_with_corr_clean) %>%
    addProviderTiles("Esri.WorldImagery") %>%
    addPolylines(
      color = ~color_scale(mechanism),
      weight = 3,
      opacity = 0.7,
      popup = ~paste0("ID: ", ID, "<br>",
                      "Mechanism: ", mechanism)
    ) %>%
    addLegend(
      "bottomright",
      pal = color_scale,
      values = ENP_transects_with_corr_clean$mechanism,
      title = "Mechanism",
      opacity = 1,  # ensure color blocks are opaque
      labFormat = labelFormat(),        # keep label formatting default
      # This makes legend background solid white
      # (note: background style may depend on browser or leaflet version)
      # Add solid background for the legend container
      # Leaflet's addLegend does not allow direct background setting, but can be controlled via CSS or:
      # Use leaflet.extras::addControl if full custom style is needed
      group = NULL,
      position = "bottomright"
    ) %>%
    htmlwidgets::onRender("
      function(el, x) {
        var legend = el.querySelector('.leaflet-control .legend');
        if (legend) {
          legend.style.background = 'white';
          legend.style.opacity = '1';
        }
      }
    ")
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

########### Plot correlation versus distance ######
plot_correlation_distance <- function(data, output_dir) {
  # Create the plot
  p <- ggplot(data, aes(x = mean_correlation.salinity,
                        y = mean_correlation.hydro,
                        color = distance_m)) +
    geom_point(size = 3) +
    scale_color_gradientn(colors = c('#feedde', '#fdbe85', '#fd8d3c', '#e6550d', '#a63603')) +
    theme_minimal(base_family = "Times") +
    theme(
      panel.background = element_rect(fill = "black", color = NA),
      plot.background = element_rect(fill = "black", color = NA),
      panel.grid.major = element_line(color = "gray30"),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "white", size = 14),
      axis.title = element_text(color = "white", size = 16, face = "bold"),
      plot.title = element_text(color = "white", size = 20, face = "bold", hjust = 0.5),
      legend.background = element_rect(fill = "black"),
      legend.text = element_text(color = "white", size = 14),
      legend.title = element_text(color = "white", size = 16, face = "bold")
    ) +
    labs(
      x = "Mean Correlation with Salinity",
      y = "Mean Correlation with Hydro",
      color = "Distance (m)",
      title = "Relationship Between Correlations and Distance"
    )
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Define file path
  file_path <- file.path(output_dir, "correlation_distance_plot.tif")
  
  # Save the plot
  ggsave(filename = file_path, plot = p, width = 8, height = 6, dpi = 300, device = "tiff")
}
