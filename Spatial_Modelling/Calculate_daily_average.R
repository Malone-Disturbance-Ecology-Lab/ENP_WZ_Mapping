
##' @title This is the workflow to calculate  
##' @param line1 is the first lie
##' @param line2 is the second line
##' @param transect.n is the number of transects you want
##' @return sfc a simple feature of the transects and their lengths (length.org; Meters)
##' @author Sparkle L. Malone, Boya (Paul) Zhang
##' 

library(dplyr)
library(readr)

# Set the folder path
folder_path <- "/Volumes/Malonelab/Research/Whitezone/Hydro_Salinity/Hydro/DBhydro/WZ_sites_data/"  # Replace with your actual path

# List all CSV files in the folder
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# Read and combine all CSVs, skipping the first 34 rows
combined_df <- lapply(csv_files, function(file) {
  read_csv(file, skip = 34, col_types = cols(.default = "c"))
}) %>%
  bind_rows()