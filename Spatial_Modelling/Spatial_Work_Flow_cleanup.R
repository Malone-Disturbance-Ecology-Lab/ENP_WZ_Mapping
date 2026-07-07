
##' @title This is the workflow to generate all the results 
##' @param line1 is the first lie
##' @param line2 is the second line
##' @param transect.n is the number of transects you want
##' @return sfc a simple feature of the transects and their lengths (length.org; Meters)
##' @author Sparkle L. Malone, Boya (Paul) Zhang
##' 
##' 
rm(list = ls())
library(terra)
library(tidyverse)
library(sf)
library(sfheaders)
library(units)
library(dplyr)
library(ggplot2)
library(reshape2)
library(here)
library(rstudioapi)
library(tools)
library(maps)
library(stringr)
library(leaflet)
library(ggspatial)
library(purrr)
library(scales)
library(maptiles)
library(lme4)
library(Matrix)
library(prettymapr)  # Needed for basemap tiles
library(stats)

### Section 1 - Setting up and initial analysis ####
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

path <- get_directory() ## get current directory path
setwd(path)
source( './function_MakeVegetationZoneLayers.R')
source( './function_Make_TransectLayers.R')
source( './function_called_by_main_script.R')

rast.files <- list.files(
  path = '/Volumes/MaloneLab/Research/Whitezone/Result/Emulation/Final/CDthreshold06/CD_closest_GE',
  pattern = "^raster.*\\.tif$",
  full.names = TRUE
)

### Section 2 - Read in raster stack and transects and make plots ####
### Section 2.1 - read in raster stack ###
setwd( '/Volumes/MaloneLab/Research/Whitezone/Result/Emulation/Final/CDthreshold06/CD_closest_GE')

coast <- rast(rast.files) %>%
  terra::project("+proj=longlat +datum=WGS84", method = "near") ## "near to make sure class value is not interpolated "
terra::time(coast) <-  substr(rast.files, 8, 11) %>% as.numeric # get the time out of the list of file names

ENP.transects <- makeTransect(line1 = coast.sf, line2 =  inland.sf, transect.n= 100)
# Calculate the centroid of each transect
ENP.transects <- ENP.transects %>%
  mutate(centroid = st_centroid(st_geometry(.)))

cent_coords <- st_coordinates(ENP.transects$centroid)

# Extract latitude and longitude from the centroid
ENP.transects <- ENP.transects %>%
  mutate(
    lon = cent_coords[, 1],
    lat = cent_coords[, 2]
  ) %>%
  dplyr::select(-centroid)  # remove if centroid is not needed

### Section 2.2 - plot time series of land cover change ###
# plot_land_cover_time_series("/Volumes/Malonelab/Research/Whitezone/Result/Emulation/Final/CDthreshold06/CD_closest_GE/Subregion_new/Region_2", 
#                             "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/land_cover_time_series_eastern.tif")

plot_land_cover_time_series("/Volumes/Malonelab/Research/Whitezone/Result/Emulation/Final/Roy/",
                            "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/roy_whole.tif", NULL, 20, 105, 1984)

plot_land_cover_time_series("/Volumes/Malonelab/Research/Whitezone/Result/Emulation/Final/Roy/Subregion_new/Region_0/", 
                            "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/roy_western.tif", NULL, 0, 45, 1984)

## Need to move this part - this belongs to projection #######
# plot_land_cover_time_series_extend(
#   raster_folder_path = "/Volumes/Malonelab/Research/Whitezone/Result/Emulation/Final/CDthreshold06/CD_closest_GE",
#   #raster_folder_path = "/Volumes/MaloneLab/Research/Whitezone/Result/Emulation/Final/CDthreshold06/CD_closest_GE/Subregion_new/Region_2",
#   output_prefix = "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/land_cover_time_series_extended.tif",
#   extend_years = c(2049, 2074)
# )

plot_land_cover_time_series_extend(
  raster_folder_path = "/Volumes/Malonelab/Research/Whitezone/Result/Emulation/Final/Roy",
  output_prefix = "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/Roy_Results/land_cover_time_series_extended.tif",
  extend_years = c(2049, 2074)
)

plot_boundary_and_basemap(
  shapefile_path = "/Volumes/MaloneLab/Research/Whitezone/Result/Spatial_Modelling/north_interior.shp",
  output_tif_path = "/Users/bz294/Documents/WhiteZone/Figure/Backup/interior_boundary.tif"
)

plot_boundary_and_basemap(
  shapefile_path = "/Volumes/MaloneLab/Research/Whitezone/Result/Spatial_Modelling/north_interior.shp",
  output_tif_path = "/Users/bz294/Documents/WhiteZone/Figure/Backup/interior_boundary.tif",
  transect_shapefile_path = "/Volumes/MaloneLab/Research/Whitezone/Spatial_Data/transects_37_61.shp"
)

######## Plot the boundary lines ###########################
## Section 2.3 Output the transects for splitting the SESE # 
## Transects 37 and 63
ENP.transects <- st_as_sf(ENP.transects)

# filter for ID 37 and 63
subset_transects <- ENP.transects %>%
  filter(ID %in% c(37, 61))

subset_transects_utm <- st_transform(subset_transects, 32617)
# write to shapefile
st_write(subset_transects, 
         "/Volumes/MaloneLab/Research/Whitezone/Spatial_Data/transects_37_61.shp",
         delete_dsn = TRUE)  # overwrite if exists

############## Section 3 Prepare the land-cover proportion data (summary_data_df) and transects ####
## load the results, instead of running the four steps above. 
#load("/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/results.RData")

load("/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/results_Roy.RData")
# Step 3.1: Combine all land cover change data into one data frame
summary_data_df <- do.call(rbind, results)
summary_data_df$year <- as.double(summary_data_df$year)
colnames(summary_data_df) <- c("ID", "time", "proportion_Int1", "proportion_Int2", "proportion_Int3")

## Step 3.2 extract transects ##
aoi <- st_read("/Volumes/MaloneLab/Research/Whitezone/Data/SEregion_reprojected_revise.shp")
aoi <- st_transform(aoi, st_crs(ENP.transects))
transects_in_polygon <- st_intersection(ENP.transects, aoi)
transects_in_polygon <- transects_in_polygon[, !(names(transects_in_polygon) == "id")]

#################### Section 4. Rates calculation #####################################
#### 4.1 calculate rates for DZ and VZ
temp <- summary_data_df %>%
  arrange(ID, time) %>%
  group_by(ID) %>%
  summarise(
    trend_Int1 = last(proportion_Int1) - first(proportion_Int1),
    trend_Int2 = last(proportion_Int2) - first(proportion_Int2),
    trend_Int3 = last(proportion_Int3) - first(proportion_Int3),
    .groups = "drop"
  )
  # summarise(
  #   prop_Int1_2014 = proportion_Int1[time == 2014],
  #   prop_Int2_2014 = proportion_Int2[time == 2014],
  #   prop_Int3_2014 = proportion_Int3[time == 2014],
  #   trend_Int1 = last(proportion_Int1) - proportion_Int1[time == 2014],
  #   trend_Int2 = last(proportion_Int2) - proportion_Int2[time == 2014],
  #   trend_Int3 = last(proportion_Int3) - proportion_Int3[time == 2014],
  #   .groups = "drop"
  # )


# Join the two data frames by ID
rate_results <- transects_in_polygon %>%
  mutate(length.org = as.numeric(length.org)) %>%  # Drop units here
  left_join(temp, by = "ID") %>%
  mutate(
    length_change_Int1 = length.org * trend_Int1 / 40,
    length_change_Int2 = length.org * trend_Int2 / 40,
    length_change_Int3 = length.org * trend_Int3 / 40
  )

output_path <- "/Volumes/MaloneLab/Research/Whitezone/Result/Spatial_Modelling/transects_with_trend_.shp"
st_write(rate_results, output_path, delete_layer = TRUE)

## 4.2 Plot the scaled rates for visualization ####
result <- plot_rates(rate_results, length_col = "length_change_Int2")
rescaled_transects <- result$transects
connection_line <- result$connection_line

### Simulation ####
# set.seed(123)  # for reproducibility
# rate_results <- rate_results %>%
#   mutate(
#     simulation = ifelse(
#       ID < 61,
#       rnorm(n(), mean = 40, sd = 12),   # positive
#       rnorm(n(), mean = -30, sd = 10)    # negative
#     )
#   )
# result <- plot_rates(rate_results, length_col = "simulation")
# rescaled_simulation <- result$transects
# # Add a random value to each transect (0–1) to map color
# rescaled_simulation <- rescaled_simulation %>%
#   mutate(rand_color = runif(n()))
###################

result_vz <- plot_rates(rate_results, length_col = "length_change_Int1")
rescaled_vz <- result_vz$transects

plot_and_save_transects(rate_results, connection_line, subset_transects_utm, filename = "/Users/bz294/Documents/WhiteZone/Figure/Backup/transects.tiff")

# Define output path
output_path <- "/Users/bz294/Documents/WhiteZone/Figure/Backup/simulated_rates.tif"
output_path <- "/Users/bz294/Documents/WhiteZone/Figure/Backup/simulated_rates_nocolor.tif"
# Create the plot

# Plot
p <- ggplot() +
  #geom_sf(data = rescaled_simulation, aes(color = ID), size = 0.8) +
  geom_sf(data = rescaled_simulation, color="gray50", size = 0.8) +
  geom_sf(data = connection_line, color = "black", size = 1) +
  ggtitle("Rescaled Transects with Connection Line") +
  # scale_color_gradient(
  #   low = "#7f3f00",    # deep/dark
  #   high = "#f7f1e0",   # light
  #   limits = c(1, 100), # map IDs 1–100 to the gradient
  #   breaks = c(20, 57.5, 87.5), # midpoint for visual reference
  #   labels = c("<40", "40-75", ">75"),
  #   na.value = "gray"
  # ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

################# Generate a new simulation for the VZ #############
# set.seed(123)  # for reproducibility
# rate_results <- rate_results %>%
#   mutate(
#     simulation = ifelse(
#       ID < 61,
#       rnorm(n(), mean = -30, sd = 15),   # positive
#       rnorm(n(), mean = 40, sd = 18)    # negative
#     )
#   )
# result <- plot_rates(rate_results, length_col = "simulation")
# rescaled_simulation <- result$transects
# # Add a random value to each transect (0–1) to map color
# rescaled_simulation <- rescaled_simulation %>%
#   mutate(rand_color = runif(n()))

# p <- ggplot() +
#   geom_sf(data = rescaled_simulation, aes(color = ID), size = 0.8) +
#   geom_sf(data = connection_line, color = "black", size = 1) +
#   ggtitle("Rescaled Transects with Connection Line") +
#   scale_color_gradient(
#     low = "#00441b",    # deep/dark
#     high = "#a1d996",   # light
#     limits = c(1, 100), # map IDs 1–100 to the gradient
#     breaks = c(20, 57.5, 87.5), # midpoint for visual reference
#     labels = c("<40", "40-75", ">75"),
#     na.value = "gray"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text = element_blank(),
#     axis.title = element_blank(),
#     axis.ticks = element_blank(),
#     panel.grid = element_blank(),
#     panel.background = element_rect(fill = "white", color = NA),
#     plot.background = element_rect(fill = "white", color = NA)
#   )
# # Save as TIFF
# ggsave(
#   filename = output_path,
#   plot = p,
#   device = "tiff",
#   width = 6,      # in inches
#   height = 6,     # in inches
#   dpi = 300       # resolution
# )


#### 4.3 regional rates for DZ and VZ ######
calculate_length_change_averages(rate_results, "length_change_Int2", custom_range = c(27, 37))
calculate_length_change_averages(rate_results, "length_change_Int1", custom_range = c(27, 37))
 
####Section 5 - Generating interior & exterior boundary of WZ shift 1984-2024######

# Section 5.1 Generating interior boundary shift line shapefiles
#load("/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/north_point.RData")
#load("/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/north_point_Roy.RData")

load("/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/north_point_Roy_VZ.RData")
#load("/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/south_point_Roy.RData")
#final_results <- south_results

#dz_north <- load("/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/north_point_Roy.RData")
points_1984 <- final_results %>% filter(Year == 1984)
points_1994 <- final_results %>% filter(Year == 1994)
points_2024 <- final_results %>% filter(Year == 2024)

points_1984_clean <- points_1984 %>% filter(!is.na(X) & !is.na(Y))
points_2024_clean <- points_2024 %>% filter(!is.na(X) & !is.na(Y))
points_1994_clean <- points_1994 %>% filter(!is.na(X) & !is.na(Y))

sf_1984 <- st_as_sf(points_1984_clean, coords = c("X", "Y"), crs = 32617)
sf_2024 <- st_as_sf(points_2024_clean, coords = c("X", "Y"), crs = 32617)
sf_1994 <- st_as_sf(points_1994_clean, coords = c("X", "Y"), crs = 32617)

st_write(sf_1984, "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/points_1984_Roy_VZ.shp", delete_layer = TRUE)
st_write(sf_2024, "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/points_2024_Roy_VZ.shp", delete_layer = TRUE)

sf_1994$Polygon_ID <- as.integer(sf_1994$Polygon_ID)
sf_1984$Polygon_ID <- as.integer(sf_1984$Polygon_ID)
sf_2024$Polygon_ID <- as.integer(sf_2024$Polygon_ID)

line_sf <- calculate_signed_distance(
  sf_1984 = sf_1984,
  #sf_2024 = sf_2024,
  sf_2024 = sf_2024,
  output_dir = "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/",
  filename = "south_interior_Roy_VZ.shp"
)

line_sf$annual_rate_m <- line_sf$signed_distance_m / 40

# Section 5.2 Generating exterior boundary shift line shapefiles
south_1984 <- south_results %>% filter(Year == 1984)
south_2024 <- south_results %>% filter(Year == 2024)
south_1984_clean <- south_1984 %>% filter(!is.na(X) & !is.na(Y))
south_2024_clean <- south_2024 %>% filter(!is.na(X) & !is.na(Y))
sf_1984_s <- st_as_sf(south_1984_clean, coords = c("X", "Y"), crs = 32617)
sf_2024_s <- st_as_sf(south_2024_clean, coords = c("X", "Y"), crs = 32617)
sf_1984_s$Polygon_ID <- as.integer(sf_1984_s$Polygon_ID)
sf_2024_s$Polygon_ID <- as.integer(sf_2024_s$Polygon_ID)
line_sf_s <- calculate_signed_distance(
  sf_1984 = sf_1984_s,
  #sf_2024 = sf_2024,
  sf_2024 = sf_2024_s,
  output_dir = "/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/",
  filename = "south_interior.shp"
)

line_sf_s$annual_rate_m <- line_sf_s$signed_distance_m / 40

##### calculate distance to coast ###

#==================================================
# 1. Convert ENP transects from WGS84 to UTM 17N
#==================================================

ENP.transects_utm <- st_read(
  "/Volumes/Malonelab/Research/Whitezone/Data/Intersected_Transects.shp",
  quiet = TRUE
) |>
  st_transform(32617)


#==================================================
# 2. Remove rows with missing coordinates
#==================================================

points_1984_clean <- points_1984 %>%
  filter(!is.na(X), !is.na(Y))


#==================================================
# 3. Convert points to sf
#    (also UTM 17N)
#==================================================

points_sf <- st_as_sf(
  points_2024_clean,
  coords = c("X", "Y"),
  crs = 32617
)


#==================================================
# 4. Function to extract southern-most point
#==================================================

get_south_point <- function(line) {
  
  coords <- st_coordinates(line)
  
  # minimum northing = southern-most
  south_coord <- coords[which.min(coords[,2]), 1:2]
  
  return(st_point(south_coord))
}


#==================================================
# 5. Extract southernmost point of every transect
#==================================================

south_transects <- st_sf(
  ID = ENP.transects_utm$ID,
  geometry = st_sfc(
    lapply(
      st_geometry(ENP.transects_utm),
      get_south_point
    ),
    crs = st_crs(ENP.transects_utm)
  )
)


#==================================================
# 6. Match Polygon_ID with transect ID
#==================================================

idx <- match(
  points_sf$Polygon_ID,
  south_transects$ID
)


# remove points without corresponding transects
valid <- !is.na(idx)

points_sf <- points_sf[valid, ]

south_geom <- st_geometry(south_transects)[idx[valid]]


#==================================================
# 7. Calculate distance (meters)
#==================================================

points_sf$distance_south_m <- as.numeric(
  st_distance(
    points_sf,
    south_geom,
    by_element = TRUE
  )
)


#==================================================
# 8. Average distance
#==================================================

mean_distance <- mean(
  points_sf$distance_south_m,
  na.rm = TRUE
)

print(mean_distance)


#==================================================
# 9. Back to dataframe if needed
#==================================================

points_1984_clean <- points_sf %>%
  st_drop_geometry()


#==================================================
# 8. Convert back to dataframe if needed
#==================================================

points_1984_clean <- points_sf %>%
  st_drop_geometry()


###################### Section 6 - processing salinity data ###########
### Section 6.1 read in the salinity data
salinity_data <- read.csv("/Volumes/MaloneLab/Research/Whitezone/Hydro_Salinity/Salinity/salinity_annual_avg_full.csv")
selected_sites <- c(
  "FLAB12_salinity",
  "FLAB11_salinity",
  "FLAB10_salinity",
  "ENPTR_salinity",
  "ENPLM_salinity",
  "ENPTC_salinity",
  "FLAB08_salinity",
  "FLAB07_salinity",
  "ENPDK_salinity",
  "ENPBS_salinity",
  "ENPLB_salinity",
  "ENPLS_salinity",
  "ENPHC_salinity",
  "FLAB04_salinity"
)

all_site_names <- salinity_data$site[salinity_data$site %in% selected_sites]

# Section 6.2: Filter summary_data_df for the years that match ###
matching_years <- salinity_data$year  ## 1991 onwards 

summary_data <- summary_data_df %>%
  filter(time %in% matching_years)

summary_data$time <- as.double(summary_data$time)

#### 6.3 Assign the sites to a region (west, central, and east) ##
west <- c(
  "FLAB12_salinity",
  "FLAB11_salinity",
  "FLAB10_salinity",
  "ENPTR_salinity",
  "ENPLM_salinity",
  "ENPTC_salinity"
)

central <- c(
  "FLAB08_salinity",
  "FLAB07_salinity",
  "ENPDK_salinity",
  "ENPBS_salinity",
  "ENPLB_salinity",
  "ENPLS_salinity",
  "ENPHC_salinity"
)

east <- c(
  "FLAB04_salinity"
)
######### Statistical model on the absolute values of salinity among subregions##
salinity_data <- salinity_data %>%
  mutate(subregion = case_when(
    site %in% west ~ "West",
    site %in% central ~ "Central",
    site %in% east ~ "East",
    TRUE ~ NA_character_  # if any site is not in the lists
  ))

# Compute mean salinity per site
site_means <- salinity_data %>%
  group_by(site, subregion) %>%
  summarise(mean_salinity = mean(annual_avg, na.rm = TRUE), .groups = "drop")

# Compare across subregions using ANOVA
aov_model <- aov(mean_salinity ~ subregion, data = site_means)
summary(aov_model)

# Post-hoc pairwise comparisons (Tukey test)
TukeyHSD(aov_model)

####### Statistical analysis on the rate differences among the groups ##
summary_data$ID_group <- cut(summary_data$ID,
                             breaks = c(-Inf, 36, 61, Inf),
                             labels = c("<37","37-61",">61"))

# Fit model
model <- lmer(proportion_Int1 ~ time * ID_group + (1|ID), data = summary_data)
summary(model)
##### Section 6.4 Calculate correlations by region #################
correlation_salinity_dz <- compute_site_correlations_by_region(
  all_annual_averages = salinity_data,
  summary_data = summary_data,
  west_sites = west,
  central_sites = central,
  east_sites = east,
  proportion_col = "proportion_Int2"  # STRING
)

##### Section 6.4 Join with the lat and lon info #################
csv_file <- "/Volumes/Malonelab/Research/Whitezone/Hydro_Salinity/Salinity/WZ_salinity_sites_latlon.csv"  # replace with your file path
salinity_summary <- read.csv(csv_file)
salinity_summary <- na.omit(salinity_summary)  # Remove rows with N/A values

# Assuming the 2nd and 3rd columns are Lat and Lon in DDMMSS.SS format
salinity_summary <- salinity_summary %>% ## this data is the records for salinity and lat/lon
  mutate(
    lat_decimal = sapply(salinity_summary[[3]], convert_to_decimal),  # Convert lat to decimal degrees
    lon_decimal = sapply(salinity_summary[[4]], convert_to_decimal)  # Convert lon to decimal degrees
  )
salinity_summary$lon_decimal <- -abs(salinity_summary$lon_decimal)
export_flab_stations(salinity_summary, '/Volumes/MaloneLab/Research/Whitezone/Result/Spatial_Modelling/Results/')

#### 6.5 Plot DZ correlation with salinity based on rates ######
# Join correlation to the transects
rescaled_transects_colored <- rescaled_transects %>%
  left_join(correlation_salinity_dz, by = "ID")

# Extract coordinates
rescaled_coords_dz <- st_coordinates(rescaled_transects_colored) %>%
  as.data.frame()

# Now bind attributes back from the sf object
rescaled_coords_dz <- rescaled_coords_dz %>%
  dplyr::mutate(row_id = L1) %>%   # match rows to features
  dplyr::left_join(
    rescaled_transects_colored %>% 
      tibble::rowid_to_column("row_id") %>%   # row_id aligns with L1
      st_drop_geometry(),                     # drop geometry but keep attributes
    by = "row_id"
  )

transect_plot_dz <- ggplot() +
    geom_path(
      data = rescaled_coords_dz,
      aes(X, Y, group = ID, color = correlation),
      linewidth = 3
    ) +
    geom_sf(data = connection_line, color = "black", linewidth = 2) +
    geom_sf(
      data = subset_transects_utm,
      color = "black",
      linewidth = 3,
      linetype = "11"
    ) +
    scale_color_gradientn(
      colors = c("#c6dbef", "#6baed6", "#08519c"),
      values = scales::rescale(c(-0.6, 0, 0.6)),
      limits = c(-0.6, 0.6),
      breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
      oob = scales::squish,
      na.value = "grey80",
      name = "Correlation"
    ) +
    coord_sf(datum = NA, expand = FALSE) +
    theme_classic() +
    theme(
      # white backgrounds
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      
      # remove grids
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # remove axes
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      
      # remove border
      panel.border = element_blank(),
      
      # legend
      legend.position = "bottom",
      legend.title = element_text(size = 18, family = "Times"),
      legend.text = element_text(size = 18, family = "Times"),
      legend.key.height = unit(1.2, "cm"),
      legend.key.width = unit(0.5, "cm")
    ) +
    guides(
      color = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        barwidth = unit(15, "cm"),
        barheight = unit(0.8, "cm")
      )
    )

ggsave(
  filename = "/Users/bz294/Documents/WhiteZone/Figure/Backup/rescaled_dz_plot_roy.tif",
  plot = transect_plot_dz,   # use the named plot object
  device = "tiff",
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"            # ensures white background
)

############# Plot the transect land-cover data ######
# Make sure year is numeric
summary_data <- summary_data %>%
  mutate(time = as.numeric(as.character(time)))

# Filter for ID 86 and years 1991-2024
df_86 <- summary_data %>%
  filter(ID == 86, time >= 1991, time <= 2024)

# Plot the time series
ggplot(df_86, aes(x = time, y = proportion_Int2)) +
  geom_line(color = "blue") +  # muted yellow
  geom_point(color = "red") +
  labs(title = "Time Series of Proportion_Int2 for ID 86 (1991-2024)",
       x = "Year",
       y = "Proportion_Int2") +
  scale_x_continuous(breaks = c(1991, 2002, 2013, 2024)) +
  theme_minimal()

##### Plot for the time series of annual mean across the eastern subregion ##
# Filter for East subregion and years 1991-2024
west_df <- salinity_data %>%
  filter(subregion == "West", year >= 1991, year <= 2024)

# Calculate yearly average across all sites
annual_mean <- west_df %>%
  group_by(year) %>%
  summarise(mean_annual = mean(annual_avg, na.rm = TRUE))

# Plot the time series
# Create the plot and assign it to a variable
ggplot(annual_mean, aes(x = year, y = mean_annual)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  labs(title = "Annual Mean Salinity in Western Subregion",
       x = "Year",
       y = "Mean Annual Salinity (PSU)") +
  scale_x_continuous(breaks = c(1991, 2002, 2013, 2024)) +
  theme_minimal()



#### 6.6 Plot VZ correlation with salinity based on rates ######
# Join correlation to the transects

correlation_salinity_vz <- compute_site_correlations_by_region(
  all_annual_averages = salinity_data,
  summary_data = summary_data,
  west_sites = west,
  central_sites = central,
  east_sites = east,
  proportion_col = "proportion_Int1"  # STRING
)

rescaled_vz_colored <- rescaled_vz %>%
  left_join(correlation_salinity_vz, by = "ID")

# Extract coordinates
rescaled_coords_vz <- st_coordinates(rescaled_vz_colored) %>%
  as.data.frame()

# Now bind attributes back from the sf object
rescaled_coords_vz <- rescaled_coords_vz %>%
  dplyr::mutate(row_id = L1) %>%   # match rows to features
  dplyr::left_join(
    rescaled_vz_colored %>% 
      tibble::rowid_to_column("row_id") %>%   # row_id aligns with L1
      st_drop_geometry(),                     # drop geometry but keep attributes
    by = "row_id"
  )

# transect_plot_vz <- ggplot() +
#   geom_path(data = rescaled_coords_vz, aes(X, Y, group = ID, color = correlation), size = 3) +
#   geom_sf(data = connection_line, color = "black", size = 2) +
#   geom_sf(data = subset_transects_utm, color = "blue", linewidth = 2, linetype = "dashed") +
#   scale_color_gradient(low = "#33a354",high = "#00441b", na.value = "gray") +
#   coord_sf(datum = NA) +
#   theme_minimal() +
#   theme(
#     panel.grid = element_blank(),
#     axis.text = element_blank(),
#     axis.title = element_blank(),
#     axis.ticks = element_blank(),
#     legend.title = element_text(size = 18, family = "Times"),
#     legend.text = element_text(size = 18, family = "Times"),
#     legend.key.height = unit(1.2, "cm"),  # make legend taller (optional)
#     legend.key.width = unit(0.5, "cm")    # adjust width (optional)
#   )

transect_plot_vz <- ggplot() +
  geom_path(data = rescaled_coords_vz, aes(X, Y, group = ID, color = correlation), size = 3) +
  geom_sf(data = connection_line, color = "black", size = 2) +
  geom_sf(data = subset_transects_utm, color = "black", linewidth = 3, linetype = "11") +
  scale_color_gradientn(
    #colors = c("#7f3f00", "#d9d9d9", "#33a354"),  # negative → neutral → positive (green)
    colors = c("#c6dbef", "#6baed6", "#08519c"),
    values = rescale(c(
      min(rescaled_coords_vz$correlation, na.rm = TRUE),
      0,
      max(rescaled_coords_vz$correlation, na.rm = TRUE)
    )),
    na.value = "gray",
    name = "Correlation"
  ) +
  coord_sf(datum = NA) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 18, family = "Times"),
    legend.text = element_text(size = 18, family = "Times"),
    legend.key.height = unit(1.2, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.position = "bottom"  # move legend to bottom
  ) +
  guides(
    color = guide_colorbar(
      direction = "horizontal",  # horizontal legend
      title.position = "top",    # title above the bar
      barwidth = 15,             # width of the color bar
      barheight = 1              # height of the color bar
    )
  )

transect_plot_vz <- ggplot() +
  geom_path(
    data = rescaled_coords_vz,
    aes(X, Y, group = ID, color = correlation),
    size = 3
  ) +
  geom_sf(
    data = connection_line,
    color = "black",
    size = 2
  ) +
  geom_sf(
    data = subset_transects_utm,
    color = "black",
    linewidth = 3,
    linetype = "11"
  ) +
  scale_color_gradientn(
    colors = c("#c6dbef", "#6baed6", "#08519c"),
    limits = c(-0.6, 0.6),
    values = scales::rescale(c(-0.6, 0, 0.6)),
    breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
    labels = c("-0.6", "-0.3", "0", "0.3", "0.6"),
    na.value = "gray",
    name = "Correlation",
    oob = scales::squish
  ) +
  coord_sf(datum = NA) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 18, family = "Times"),
    legend.text = element_text(size = 14, family = "Times"),
    legend.key.height = grid::unit(0.8, "cm"),
    legend.key.width = grid::unit(2.5, "cm"),
    legend.position = "bottom"
  ) +
  guides(
    color = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 12,
      barheight = 0.8
    )
  )
ggsave(
  filename = "/Users/bz294/Documents/WhiteZone/Figure/Backup/rescaled_vz_plot_roy.tif",
  plot = transect_plot_vz,   # use the named plot object
  device = "tiff",
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"            # ensures white background
)


########## Statistical test on subregional differences of correlation ######
# 1. Remove rows with duplicated IDs (remove the first occurrence)
rescaled_coords_clean <- rescaled_coords_vz %>%
  group_by(ID) %>%
  mutate(row_number = row_number()) %>%
  filter(row_number != 1) %>%  # remove the first row for each ID
  ungroup() %>%
  dplyr::select(-row_number)

# 2. Create grouping variable based on ID ranges
rescaled_coords_clean <- rescaled_coords_clean %>%
  mutate(group = case_when(
    ID < 37 ~ "<37",
    ID >= 37 & ID <= 61 ~ "37-61",
    ID > 61 ~ ">61"
  ))

# Check the grouping
table(rescaled_coords_clean$group)

# 3. Perform ANOVA to test for differences in correlation among groups
anova_result <- aov(correlation ~ group, data = rescaled_coords_clean)
summary(anova_result)

# 4. Post-hoc pairwise comparisons (Tukey HSD)
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)


####################### Section 8. Prepare for data analysis #######
salinity_data$site <- as.factor(salinity_data$site)
salinity_data <- salinity_data %>%
  rename( annual_avg = mean_value)
salinity_data <- salinity_data %>%
  rename( site = Station_ID)
salinity_data$site <- str_remove(salinity_data$site, "_salinity$")

wl <- wl %>%
  rename(site = source_file)
wl$site <- as.factor(wl$site)
wl <- wl %>%
  mutate(
    mean_value = mean_value * 30.48,  # Convert to cm
  )


########### Section 7 Water level correlation analysis #################################
###### Section 7.1 Read the wl csv file ########
wl <- read.csv("/Volumes/MaloneLab/Research/Whitezone/Hydro_Salinity/Hydro/Merged/annual_summary_merged.csv", stringsAsFactors = FALSE)
all_site_names <- wl$source_file
# Filter salinity_data for the top 5 sites
wl_data <- wl %>%
  select(year, source_file, mean_value, sd_value)

# Step 2: Filter summary_data_df for the years that match
matching_years <- wl_data$year  # Get the years from top_5_sites

total_sites <- length(all_site_names)
summary_data_df$time <- as.double(summary_data_df$time)
correlation_dbh_df <- compute_dbhydro_correlations(
  all_site_names = all_site_names,
  wl = wl,
  summary_data = summary_data_df
)


##### Section 7.2 read in the lat and lon info and join the data frame ######
csv_file <- "/Volumes/MaloneLab/Research/Whitezone/Hydro_Salinity/Hydro/Merged/WZ_hydro_sites_latlon_merged.csv"  # replace with your file path
db_summary <- read.csv(csv_file)
db_summary <- db_summary[!duplicated(db_summary[[1]]), ]

# Assuming the 2nd and 3rd columns are Lat and Lon in DDMMSS.SS format
db_summary <- db_summary %>% ## this data is the records for salinity and lat/lon
  mutate(
    lat_decimal = sapply(db_summary[[2]], convert_to_decimal),  # Convert lat to decimal degrees
    lon_decimal = sapply(db_summary[[3]], convert_to_decimal)  # Convert lon to decimal degrees
  )
db_summary$lon_decimal <- -abs(db_summary$lon_decimal)
db_summary <- db_summary %>%
  mutate(Station_ID = Name)



# Define the lists
west <- c("NP-P37", "NP-146", "NP-OL", "Mud_Creek_at_mouth", "SWEVER5A", "NP-CHP")
central <- c("EVER7", "CT50R", "EVER6", "NP-EPSW", "EPSW", "CT27R", "NP-EP1", "EVER8",
             "HIGHWAY_CR", "ENPLS")
east <- c("SWEVER1", "S197_H", "S197_T", "S20_T", "S20_H")

# Combine all valid names
all_names <- c(west, central, east)

# 1. Remove rows with Name not in the lists
db_summary <- db_summary %>%
  filter(Name %in% all_names)

# 2. Check if all names in the lists are present in db_summary
missing_names <- setdiff(all_names, db_summary$Name)

if(length(missing_names) == 0){
  cat("All names are included in db_summary.\n")
} else {
  cat("These names are missing from db_summary:\n")
  print(missing_names)
}
write.csv(db_summary, "/Volumes/MaloneLab/Research/Whitezone/Hydro_Salinity/Hydro/Merged/db_summary_filtered.csv", row.names = FALSE)
############# Make plot of salinity & wl gauges ##############################

# Save
ggsave(
  filename = "/Users/bz294/Documents/WhiteZone/Figure/Sup/points_map.tif",
  plot = p,
  width = 8,
  height = 6,
  dpi = 600,
  bg = "white"
)


model_result_salinity <- correlation_salinity_df %>%
  filter(!grepl("^C", Station_ID)) %>%  # Exclude Station_IDs starting with "C"
  group_by(ID) %>%
  summarise(
    mean_correlation = mean(correlation, na.rm = TRUE),
    std_correlation = sd(correlation, na.rm = TRUE),
    .groups = "drop"
  )

### combined correlation_dbh_df and correlation_eden_df
correlation_hydro_comb <- rbind(correlation_dbh_df, correlation_eden_df)

model_result_hydro <- correlation_hydro_comb %>%
  group_by(ID) %>%
  summarise(
    mean_correlation = mean(correlation, na.rm = TRUE),
    std_correlation = sd(correlation, na.rm = TRUE),
    .groups = "drop"
  )

####### calculate area for a polygon shapefile ####
# Load required package
library(sf)

# 1. Read the shapefile
# Replace "path/to/your/shapefile.shp" with the actual file path
shp <- st_read("/Users/bz294/Documents/WhiteZone/Points/Shapefiles/SEregion_reprojected_revise.shx")

# 2. Make sure the shapefile has a projected coordinate system (in meters)
# If it's in geographic coordinates (degrees), reproject to a suitable projection (e.g., UTM)
if (st_is_longlat(shp)) {
  shp <- st_transform(shp, crs = 32617) # Example: UTM zone 17N (change to your region)
}

# 3. Calculate the area of each polygon (in square meters)
shp$area_m2 <- st_area(shp)

# 4. Optionally, calculate area in hectares
shp$area_ha <- as.numeric(shp$area_m2) / 10000

# 5. View results
print(shp)




########### Output shapefiles for eden and dbhdyro #####
export_points_shapefile(
  df = eden_summary,
  lon_col = "lon_decimal",
  lat_col = "lat_decimal",
  output_path = "/Volumes/MaloneLab/Research/Whitezone/Result/Spatial_Modelling/Results/eden_summary_points.shp"
)
