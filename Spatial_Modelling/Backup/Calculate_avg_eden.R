
##' @title This is the workflow to calculate annual average water level
##' @author Boya (Paul) Zhang
##' 

library(dplyr)
library(readr)

folder_path <- "/Volumes/Malonelab/Research/Whitezone/Hydro_Salinity/Hydro/EDEN/data"

# List all CSV files
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# Read each CSV, using the 34th row as header
data_list <- lapply(csv_files, function(file) {
  read.csv(file, skip = 4, stringsAsFactors = FALSE)
})

# Combine into one dataframe (assumes column names match across files)
combined_df <- dplyr::bind_rows(data_list)

summary_list <- lapply(csv_files, function(file) {
  df <- read.csv(file, skip = 4, stringsAsFactors = FALSE)
  
  # Check if at least two columns exist
  if (ncol(df) < 2) return(NULL)
  
  # Rename the first two columns for convenience
  names(df)[1:2] <- c("timestamp_col", "value_col")
  
  # Parse the first column as date/time
  df$timestamp_col <- parse_date_time(df$timestamp_col, orders = c("ymd", "mdy", "dmy", "Ymd"))
  
  # Filter out rows with invalid dates or NA in value column
  df <- df[!is.na(df$timestamp_col) & !is.na(df$value_col), ]
  
  # Extract year from parsed date
  df$year <- year(df$timestamp_col)
  
  # Summarize mean and sd by year
  summary <- df %>%
    group_by(year) %>%
    summarise(
      mean_value = mean(as.numeric(value_col), na.rm = TRUE),
      sd_value = sd(as.numeric(value_col), na.rm = TRUE),
      .groups = "drop"
    )
  
  summary$source_file <- basename(file)
  return(summary)
})

annual_summary_df <- bind_rows(summary_list)
annual_summary_df$`source_file` <- sub("\\.csv$", "", annual_summary_df$`source_file`)

eden <- read.csv("/Volumes/Malonelab/Research/Whitezone/Hydro_Salinity/Hydro/EDEN/WZ_hydro_latlon.csv",
                 stringsAsFactors = FALSE,
                 fileEncoding = "latin1")

# Combine lat and lat_min into decimal degrees
eden$lat_decimal <- eden$lat + sapply(eden$lat_min, convert_minsec)

# Combine lon and lon_min into decimal degrees (note: negative for West)
eden$lon_decimal <- eden$lon - sapply(eden$lon_min, convert_minsec)

annual_summary_with_coords <- annual_summary_df %>%
  left_join(eden, by = c("source_file" = "WZ_sites"))

annual_summary_with_coords <- annual_summary_with_coords %>%
  select(-lat, -lon, -lat_min, -lon_min)

write.csv(annual_summary_with_coords, "/Volumes/Malonelab/Research/Whitezone/Hydro_Salinity/Hydro/EDEN/annual_summary.csv", row.names = FALSE)
