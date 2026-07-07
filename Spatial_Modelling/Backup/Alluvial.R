library(terra)
library(dplyr)
library(tidyr)
library(ggalluvial)
library(stringr)

# --- Step 1: Define raster stack directory ---
raster_dir <- "/Volumes/MaloneLab/Research/Whitezone/Result/Emulation/Final/CDthreshold06/CD_closest_GE/Subregion/Region_3"  # UPDATE THIS

# List and sort raster files
raster_files <- list.files(raster_dir, pattern = "\\.tif$", full.names = TRUE)
raster_files <- sort(raster_files)

# Optional: select specific years
selected_years <- c(1985, 1994, 2003, 2013, 2023)
all_years <- 1984:2024
selected_indices <- match(selected_years, all_years)
raster_files <- raster_files[selected_indices]

# --- Step 2: Read raster stack ---
rasters <- rast(raster_files)  # terra SpatRaster stack
rasters[rasters == 0] <- NA   # Treat 0 as NA

# --- Step 3: Extract pixel-wise trajectories ---
mat <- as.matrix(rasters)       # rows = pixels, columns = years
mat <- mat[complete.cases(mat), ]  # remove rows with any NA

# --- Step 4: Count transition patterns ---
df <- as.data.frame(mat)
colnames(df) <- as.character(selected_years)
df$count <- 1

agg_df <- df %>%
  group_by(across(all_of(as.character(selected_years)))) %>%
  summarise(count = n(), .groups = "drop")

# Create a unique trajectory identifier for each row
agg_df$trajectory <- apply(agg_df[, as.character(selected_years)], 1, paste, collapse = "-")

# --- Step 5: Prepare data for ggalluvial ---
long_df <- agg_df %>%
  pivot_longer(cols = all_of(as.character(selected_years)),
               names_to = "Year",
               values_to = "Class") %>%
  mutate(Year = as.integer(Year))

# Map numeric classes to factors with labels and set order
class_map <- c("1" = "VZ", "2" = "WZ", "3" = "TZ")
long_df$Class <- factor(class_map[as.character(long_df$Class)],
                        levels = c("TZ", "WZ", "VZ"))  # top to bottom order

# Add the trajectory column repeated for each year
long_df$trajectory <- rep(agg_df$trajectory, each = length(selected_years))

# Select top 20 largest trajectories to plot (optional, for clarity)
top_trajectories <- agg_df %>%
  arrange(desc(count)) %>%
  slice(1:20) %>%
  pull(trajectory)

# Filter long_df for those trajectories
plot_df <- long_df %>%
  filter(trajectory %in% top_trajectories)

# Now plot only with shaded flows, no boundary lines
ggplot(plot_df,
       aes(x = Year, stratum = Class, alluvium = trajectory,
           y = count, fill = Class)) +
  geom_alluvium(alpha = 0.7, color = NA) +  # shaded curves, no outlines
  scale_fill_manual(values = c("TZ" = "#ffffbf", "WZ" = "#d7191c", "VZ" = "#2b83ba")) +
  theme_minimal(base_family = "Times") +
  labs(title = "Land Cover Transitions (Top 20 Trajectories)",
       y = "Pixel Count", x = "Year") +
  theme(legend.position = "top")
