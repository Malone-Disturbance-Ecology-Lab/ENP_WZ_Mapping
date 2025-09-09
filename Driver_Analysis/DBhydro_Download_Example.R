# initial data exploration for DBHYDRO sites &
# dataframe creation for long term water level and salinity analysis (Q1)

# last updated: PZ, 5/24/2024

rm(list = ls())
getwd()
# load packages

#devtools::install_github("ropensci/dbhydroR")
#install.packages("dbhydroR")
library(dbhydroR)
library(tidyverse)        # everything else
library(terra)            # Spatial analysis to find elevation for TSPh1a and TSPh7


# coastal site ENPHC water level
# https://my.sfwmd.gov/dbhydroplsql/show_dbkey_info.date_selection?v_js_flag=Y&v_db_request_id=8977677&v_parameter_string=&v_dbkey=39307&v_frequency=&v_sdate=19930729&v_edate=20240524&v_datum=

# example:
# flab11.sal <- get_wq(station_id = c("FLAB11"), date_min = "1990-01-01",
#                      date_max = "2025-05-01", test_name = "SALINITY")

# Get salinity data
data <- read.csv("/Users/bz294/Documents/WhiteZone/Data/Driver_Analysis/Salinity/List_of_Salinity.csv")
filtered_data <- data[data[[2]] == "Yes", ]
station_names <- filtered_data[[1]]
salinity_data <- list()

# Loop through each station
for (i in seq_along(station_names)) {
  salinity_data[[i]] <- get_wq(
    station_id = station_names[i],
    date_min = "1990-01-01",
    date_max = "2025-05-01",
    test_name = "SALINITY"
  )
}

# Optionally, name the list elements by station name
names(salinity_data) <- station_names

for (i in seq_along(salinity_data)) {
  station <- names(salinity_data)[i]
  df <- salinity_data[[i]]
  write.csv(df, file = paste0(station, "_salinity_east.csv"), row.names = FALSE)
}

station_info <- get_station_metadata("FLAB11")



## Get water level data
#get_hydro





