
##' @title This is the workflow to calculate annual average water level
##' @author Boya (Paul) Zhang
##' 

library(dplyr)
library(readr)

folder_path <- "/Volumes/Malonelab/Research/Whitezone/Hydro_Salinity/Hydro/DBhydro/WZ_sites_data/"

# List all CSV files
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# Read each CSV, using the 34th row as header
data_list <- lapply(csv_files, function(file) {
  read.csv(file, skip = 33, stringsAsFactors = FALSE)
})

# Combine into one dataframe (assumes column names match across files)
combined_df <- dplyr::bind_rows(data_list)

summary_list <- lapply(csv_files, function(file) {
  df <- read.csv(file, skip = 33, stringsAsFactors = FALSE)
  
  if (!all(c("TIMESTAMP", "VALUE") %in% names(df))) return(NULL)
  
  # Parse the TIMESTAMP column
  df$TIMESTAMP <- parse_date_time(df$TIMESTAMP, orders = c("ymd", "mdy", "dmy", "Ymd"))
  
  # Filter out rows with invalid dates or NA in VALUE
  df <- df[!is.na(df$TIMESTAMP) & !is.na(df$VALUE), ]
  
  df$year <- year(df$TIMESTAMP)
  
  summary <- df %>%
    group_by(year) %>%
    summarise(
      mean_value = mean(VALUE, na.rm = TRUE),
      sd_value = sd(VALUE, na.rm = TRUE),
      .groups = "drop"
    )
  
  summary$source_file <- basename(file)
  return(summary)
})

annual_summary_df <- bind_rows(summary_list)
annual_summary_df$`source_file` <- sub("\\.csv$", "", annual_summary_df$`source_file`)

dbh <- read.csv("/Volumes/Malonelab/Research/Whitezone/Hydro_Salinity/Hydro/DBhydro/WZ_hydro_sites_latlon.csv", stringsAsFactors = FALSE)
dbh <- dbh %>% ## this data is the records for salinity and lat/lon
  mutate(
    lat_decimal = sapply(dbh[[2]], convert_to_decimal),  # Convert lat to decimal degrees
    lon_decimal = sapply(dbh[[3]], convert_to_decimal)  # Convert lon to decimal degrees
  )
dbh$lon_decimal <- -abs(dbh$lon_decimal)

annual_summary_with_coords <- annual_summary_df %>%
  left_join(dbh, by = c("source_file" = "Name"))

annual_summary_with_coords <- annual_summary_with_coords %>%
  select(-Lat, -Lon)

write.csv(annual_summary_with_coords, "/Volumes/Malonelab/Research/Whitezone/Hydro_Salinity/Hydro/DBhydro/annual_summary.csv", row.names = FALSE)
