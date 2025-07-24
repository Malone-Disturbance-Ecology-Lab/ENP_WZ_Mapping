##' @title This is the workflow to generate all the results 
##' @author  Boya (Paul) Zhang
##' @description Modelling driver mechanism of WZ expansion
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


####################### Distance analysis ###########################
# Make sure ID is treated as a factor (optional for coloring)
correlation_dbh_df <- correlation_dbh_df %>%
  filter(!Station_ID %in% c("JBTS", "MBTS", "OYSTERCK"))

correlation_salinity_df <- correlation_salinity_df %>%
  filter(!grepl("^C", site))

correlation_dbh_df <- correlation_dbh_df %>%
  mutate(ID = as.factor(ID))

correlation_salinity_df <- correlation_salinity_df %>%
  mutate(ID = as.factor(ID))
# Get all unique Station_IDs
station_ids <- unique(correlation_dbh_df$Station_ID)
station_salinity_ids <- unique(correlation_salinity_df$Station_ID)

tiff("/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/Results/correlation_vs_distance_by_station_salinity.tif", 
     width = 3000, height = 2000, res = 300)  # Adjust size/resolution as needed


# Generate the faceted scatter plot
ggplot(correlation_salinity_df, aes(x = distance_m, y = correlation)) +
  geom_point(alpha = 0.7, size = 1.5) +
  facet_wrap(~ Station_ID, scales = "free") +
  labs(
    title = "Correlation vs Distance (per Station)",
    x = "Distance (m)",
    y = "Correlation"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    panel.spacing = unit(1, "lines")
  )

# Close TIFF device to save file
dev.off()
##################### Transect based plot (not station based) #########
tiff("/Volumes/Malonelab/Research/Whitezone/Result/Spatial_Modelling/Results/correlation_eden_by_transect.tif", width = 4000, height = 4000, res = 300)
p <- ggplot(correlation_eden_df, aes(x = distance_m, y = correlation)) +
  geom_point(alpha = 0.6, size = 1) +
  facet_wrap(~ ID, scales = "free", ncol = 10) +  # 10 columns layout
  labs(
    title = "Correlation vs Distance for Each Transect",
    x = "Distance (m)",
    y = "Correlation"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 7),
    plot.title = element_text(hjust = 0.5),
    panel.spacing = unit(0.3, "lines"),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 8)
  )

# Print plot (renders in RStudio plot pane or notebook)
print(p)
dev.off()

####################### Analyzing the trend of wl and salinity change time series ###
# Convert site to factor (optional but recommended for clear facets)
salinity_data$site <- as.factor(salinity_data$site)

dbhydro <- dbhydro %>%
  rename(site = source_file)

eden <- eden %>%
  rename(site = source_file)

dbhydro$site <- as.factor(dbhydro$site)
eden$site <- as.factor(eden$site)

eden_cm <- eden %>%
  mutate(
    mean_value = mean_value * 30.48,  # Convert to cm
    site = str_replace(site, "mouth", "Mouth")  # Capitalize 'mouth' to 'Mouth'
  )

eden_cm <- eden_cm %>%
  filter(site %in% eden_summary$WZ_sites)

dbhydro_cm <- dbhydro %>%
  mutate(mean_value = mean_value * 30.48)

salinity_data <- salinity_data %>%
  rename(mean_value = annual_avg)
# Estimate wrap width: adjust divisor to tune wrapping
plot_eden_faceted(eden_cm, wrap_width = 15)

dbhydro_cm_filtered <- dbhydro_cm %>%
  filter(!str_detect(site, "MBTS|JBTS|OYSTERCK"))
plot_dbhydro_annual_avg(dbhydro_cm_filtered)

salinity_data$site <- str_remove(salinity_data$site, "_salinity$")
salinity_data <- salinity_data %>% 
  rename(mean_value = annual_avg)
plot_salinity_faceted(salinity_data)
############### Trend analysis based on site measurements ###############
correlations <- eden %>%
  filter(!site %in% c("Joe_Bay_2E", "Upstream_Taylor")) %>%  # Exclude sites
  group_by(site) %>%
  summarise(
    correlation = cor(year, mean_value, use = "complete.obs"),
    .groups = "drop"
  )

ggplot(correlations, aes(x = reorder(site, correlation), y = correlation)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Correlation between Year and Mean Value by Site",
       x = "Site",
       y = "Correlation")

site_data_counts <- dbhydro %>%
  group_by(source_file) %>%
  summarise(
    count = sum(!is.na(mean_value)),
    .groups = "drop"
  )

print(site_data_counts)

###################### Modelling the WZ expansion mechanism! ############
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
## figure out whether WZ is expanding or not?
## for salinity using summary_data_filtered, for dbh using summary_data_df
trend_result <- summary_data_df %>% ## this figures out whether WZ is expanding
  group_by(ID) %>%
  summarise(
    slope = coef(lm(proportion_Int2 ~ time))[2],
    .groups = "drop"
  ) %>%
  mutate(WZ_trend_direction = case_when(
    slope > 0 ~ 1,
    slope < 0 ~ -1,
    TRUE ~ 0  # handles slope == 0
  )) %>%
  select(ID, WZ_trend_direction)

trend_result_VZ <- summary_data_df %>% ## this figures out whether WZ is expanding
  group_by(ID) %>%
  summarise(
    slope = coef(lm(proportion_Int1 ~ time))[2],
    .groups = "drop"
  ) %>%
  mutate(VZ_trend_direction = case_when(
    slope > 0 ~ 1,
    slope < 0 ~ -1,
    TRUE ~ 0  # handles slope == 0
  )) %>%
  select(ID, VZ_trend_direction)

trend_result_value <- summary_data_df %>%
  arrange(ID, time) %>%
  group_by(ID) %>%
  summarise(
    prop_Int1_1984 = first(proportion_Int1),
    prop_Int2_1984 = first(proportion_Int2),
    prop_Int3_1984 = first(proportion_Int3),
    trend_Int1 = last(proportion_Int1) - first(proportion_Int1),
    trend_Int2 = last(proportion_Int2) - first(proportion_Int2),
    trend_Int3 = last(proportion_Int3) - first(proportion_Int3),
    .groups = "drop"
  )

# Join the two data frames by ID
length_results <- transects_in_polygon %>%
  mutate(length.org = as.numeric(length.org)) %>%  # Drop units here
  left_join(trend_result_value, by = "ID") %>%
  mutate(
    length_change_Int1 = length.org * trend_Int1,
    length_change_Int2 = length.org * trend_Int2,
    length_change_Int3 = length.org * trend_Int3
  )

output_path <- "/Volumes/MaloneLab/Research/Whitezone/Result/Spatial_Modelling/transects_with_trend.shp"
st_write(length_results, output_path, delete_layer = TRUE)

############################ Model the liklihood of change in Int1 and Int2 ######
analyze_length_changes(length_results)

############################### The end ###########################################
# Ensure all ID columns are the same type (e.g., character)
trend_result <- trend_result %>%
  mutate(ID = as.character(ID))

trend_result_VZ <- trend_result_VZ %>%
  mutate(ID = as.character(ID))

model_result_salinity <- model_result_salinity %>%
  mutate(ID = as.character(ID))

model_result_hydro <- model_result_hydro %>%
  mutate(ID = as.character(ID))

# Join with custom suffixes
combined_result <- trend_result %>%
  left_join(trend_result_VZ, by = "ID") %>%
  left_join(model_result_salinity, by = "ID", suffix = c("", ".salinity")) %>%
  left_join(model_result_hydro, by = "ID", suffix = c("", ".dbh"))

combined_result <- combined_result %>%
  rename(
    mean_correlation.hydro = mean_correlation.dbh,
    std_correlation.hydro = std_correlation.dbh,
    mean_correlation.salinity = mean_correlation,
    std_correlation.salinity = std_correlation
  )

combined_result$ID <- as.numeric(as.character(combined_result$ID))
line_sf$Polygon_ID <- as.numeric(as.character(line_sf$Polygon_ID))
WZ_change_mechanism <- merge(combined_result, line_sf, by.x = "ID", by.y = "Polygon_ID")

WZ_change_mechanism <- WZ_change_mechanism %>%
  select(-starts_with("std"))

combined_result <- combined_result %>%
  left_join(line_sf %>% select(Polygon_ID, signed_distance_m),
            by = c("ID" = "Polygon_ID"))
####################### Checking the WZ mechanism variable #########################
plot_correlation_distance(WZ_change_mechanism, "/Volumes/MaloneLab/Research/Whitezone/Result/Spatial_Modelling/Results/")
# Gather for plotting correlation vs distance for salinity and hydro in one plot
df_long <- WZ_change_mechanism %>%
  select(mean_correlation.salinity, mean_correlation.hydro, distance_m) %>%
  pivot_longer(cols = c(mean_correlation.salinity, mean_correlation.hydro), 
               names_to = "correlation_type", values_to = "correlation_value")

p <- ggplot(df_long, aes(x = distance_m, y = correlation_value, color = correlation_type)) + 
  geom_point(alpha = 0.6) + 
  geom_smooth(method = "lm", se = TRUE) + 
  scale_color_manual(values = c("mean_correlation.salinity" = "blue", "mean_correlation.hydro" = "red"),
                     labels = c("Salinity", "Hydro")) + 
  theme_minimal(base_family = "Times") + 
  labs(x = "Distance (m)", y = "Mean Correlation", color = "Correlation Type",
       title = "Correlation of Salinity and Hydro with Distance") + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))

# Save as TIFF
ggsave("/Volumes/MaloneLab/Research/Whitezone/Result/Spatial_Modelling/Results/correlation_vs_distance_new.tif", plot = p, width = 8, height = 6, dpi = 300, units = "in", device = "tiff")

lm_salinity <- lm(mean_correlation.salinity ~ distance_m, data = WZ_change_mechanism)
lm_hydro <- lm(mean_correlation.hydro ~ distance_m, data = WZ_change_mechanism)

summary(lm_salinity)
summary(lm_hydro)

thrd <- 0
combined_result <- combined_result %>%
  mutate(mechanism = case_when(
    signed_distance_m > 0 & mean_correlation.salinity > 0 & mean_correlation.hydro > thrd ~ "SLR - WZ up",
    signed_distance_m > 0 & mean_correlation.salinity > 0 & mean_correlation.hydro < thrd ~ "Water Management - WZ up",
    signed_distance_m < 0 & mean_correlation.salinity > 0 & mean_correlation.hydro < thrd ~ "Water Management - WZ down",
    signed_distance_m < 0 & mean_correlation.salinity > 0 & mean_correlation.hydro > thrd ~ "SLR - WZ down",
    TRUE ~ "Unknown or NA"
  ))

################## Make plot for the mechanism map !!!! ####################
combined_result <- combined_result %>%
  mutate(ID = as.factor(ID))

ENP.transects <- ENP.transects %>%
  mutate(ID = as.factor(ID))

ENP.transects_with_corr <- ENP.transects %>%
  left_join(combined_result, by = "ID")

ENP.transects_with_corr_clean <- ENP.transects_with_corr[!is.na(ENP.transects_with_corr$WZ_trend_direction), ]
ENP.transects_with_corr_clean <- ENP.transects_with_corr_clean %>%
  mutate(ID = as.character(ID)) %>%
  left_join(
    line_sf %>%
      st_drop_geometry() %>%
      mutate(Polygon_ID = as.character(Polygon_ID)) %>%
      select(Polygon_ID, annual_rate_m),
    by = c("ID" = "Polygon_ID")
  )
ENP.transects_with_corr_clean %>%
  select(ID, annual_rate_m, geometry) %>%
  st_write("/Volumes/MaloneLab/Research/Whitezone/Result/Spatial_Modelling/Results/ENP_with_rate.shp", delete_dsn = TRUE)

plot_mechanism_map(
  ENP_transects = ENP.transects_with_corr_clean,
  output_dir = '/Volumes/MaloneLab/Research/Whitezone/Result/Spatial_Modelling/mechanism_0.shp'
)


ggplot(data = ENP.transects_with_corr_clean) +
  geom_sf(aes(color = mechanism), size = 1) +
  scale_color_manual(values = c('#d7191c','#fdae61','#abdda4','#2b83ba')) +
  theme_minimal(base_size = 15) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),       # remove grid lines
    axis.text = element_blank(),        # remove axis text
    axis.ticks = element_blank(),       # remove axis ticks
    axis.title = element_blank()        # remove axis titles
  ) +
  labs(color = "Mechanism")

ggsave("/Volumes/MaloneLab/Research/Whitezone/Result/Spatial_Modelling/mechanism_map.tif", width = 10, height = 8, dpi = 300)
############### Making plot for all the sites and transects - Data section ####
plot_station_sources(
  eden_summary = eden_summary,
  data = salinity_summary, 
  db_summary = db_summary, 
  ENP_transects = transects_in_polygon)

############### Plot WZ trend, signs of correlation with salinity and wl ###
plot_transect_signs(ENP.transects_with_corr_clean, "WZ_trend_direction", "WZ trend")

########### Output shapefiles for eden and dbhdyro #####
export_points_shapefile(
  df = eden_summary,
  lon_col = "lon_decimal",
  lat_col = "lat_decimal",
  output_path = "/Volumes/MaloneLab/Research/Whitezone/Result/Spatial_Modelling/Results/eden_summary_points.shp"
)
