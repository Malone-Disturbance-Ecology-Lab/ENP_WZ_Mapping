##### Automatically interact with the DBhydro website to get lat and lon 
##### Given Station names
##### Author: Boya (Paul) Zhang

library(httr)
library(rvest)
library(stringr)

get_lat_lon <- function(station_name) {
  url <- "https://my.sfwmd.gov/dbhydroplsql/show_dbkey_info_detail.show_station_characters"
  
  # Simulate form submission using POST request
  response <- POST(
    url,
    body = list(
      v_station = station_name,  # name of form field
      v_target_code = "all"      # other required field
    ),
    encode = "form"
  )
  
  # Parse HTML
  page <- content(response, as = "parsed")
  
  # Extract text of the table containing Latitude and Longitude
  table_text <- page %>%
    html_nodes("table") %>%
    .[[1]] %>%                   # This assumes lat/lon is in the 3rd table
    html_text()
  
  # Use regex to find Latitude and Longitude (in ddmmss.sss format)
  lat <- str_extract(table_text, "Latitude\\s+:\\s+\\d{2,3}[0-5]{1}\\d{1}\\d{2}\\.\\d+")
  lon <- str_extract(table_text, "Longitude\\s+:\\s+-?\\d{2,3}[0-5]{1}\\d{1}\\d{2}\\.\\d+")
  
  # Clean up
  lat_clean <- str_extract(lat, "\\d{2,3}[0-5]{1}\\d{1}\\d{2}\\.\\d+")
  lon_clean <- str_extract(lon, "-?\\d{2,3}[0-5]{1}\\d{1}\\d{2}\\.\\d+")
  
  return(data.frame(Station = station_name, Latitude = lat_clean, Longitude = lon_clean))
}

stations <- c("FLAB11")
latlon_df <- do.call(rbind, lapply(stations, get_lat_lon))
print(latlon_df)