# ============================================================
# 2. EXACT COLORS FROM REFERENCE FIGURE
# ============================================================
region_colors <- c(
  "West"    = "#F8C532",  # Yellow
  "Central" = "#F68D29",  # Orange
  "East"    = "#D92529"   # Red
)

# ============================================================
# 4. LABEL POSITIONS - MOVE DOWN AND PREVENT OVERLAP
# ============================================================

bbox <- st_bbox(transects_plot)

xmin <- as.numeric(bbox["xmin"])
xmax <- as.numeric(bbox["xmax"])
ymin <- as.numeric(bbox["ymin"])
ymax <- as.numeric(bbox["ymax"])

map_width  <- xmax - xmin
map_height <- ymax - ymin

# Position labels with equal spacing
label_df <- data.frame(
  Region = factor(
    c("West", "Central", "East"),
    levels = c("West", "Central", "East")
  ),
  x = xmin + map_width * c(0.17, 0.50, 0.83),
  y = rep(ymin - 0.24 * map_height, 3),
  Label = c(
    "West\n17.6 PSU",
    "Central\n21.9 PSU",
    "East\n29.8 PSU"
  )
)

# ============================================================
# 5. CREATE UPDATED FIGURE
# ============================================================

p <- ggplot() +
  
  geom_sf(
    data = transects_plot,
    aes(color = Region),
    linewidth = 0.65,
    alpha = 0.85
  ) +
  
  geom_sf(
    data = aoi_plot,
    fill = NA,
    color = "grey35",
    linewidth = 0.65
  ) +
  
  geom_label(
    data = label_df,
    aes(
      x = x,
      y = y,
      label = Label,
      fill = Region
    ),
    color = "black",
    size = 6,
    fontface = "bold",
    family = "Times New Roman",
    label.size = 0,
    label.padding = unit(0.35, "lines"),
    lineheight = 1.15
  ) +
  
  scale_color_manual(values = region_colors) +
  scale_fill_manual(values = region_colors) +
  
  coord_sf(
    xlim = c(
      xmin - 0.05 * map_width,
      xmax + 0.05 * map_width
    ),
    ylim = c(
      ymin - 0.42 * map_height,
      ymax + 0.10 * map_height
    ),
    datum = NA,
    expand = FALSE
  ) +
  
  theme_void(base_family = "Times New Roman") +
  
  theme(
    legend.position = "none",
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    plot.margin = margin(20, 25, 20, 25)
  )

print(p)

# ============================================================
# 6. SAVE UPDATED FIGURE
# ============================================================

output_dir <- "/Users/bz294/Documents/WhiteZone/Figure/Backup"

ggsave(
  file.path(output_dir, "SESE_transects_salinity_regions.tif"),
  plot = p,
  width = 12,
  height = 7,
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

ggsave(
  file.path(output_dir, "SESE_transects_salinity_regions.png"),
  plot = p,
  width = 12,
  height = 7,
  dpi = 300,
  bg = "white"
)


## FIGURE 1 ##

# ============================================================
# EXAMPLE USAGE
# ============================================================

# result <- plot_land_cover_time_series(
#   raster_folder_path = "/path/to/annual_landcover",
#   output_path = "/path/to/Land_Cover_1984_2024.tif",
#   extend_years = NULL,
#   lowlimit = 15,
#   uplimit = 110,
#   startingyear = 1984
# )

##


library(ggplot2)
library(grid)

# ============================================================
# 1. CORRECTED DATA FROM REVISED MANUSCRIPT
#    Transect rates: 1984–2024
# ============================================================

df <- data.frame(
  Region = rep(c("West", "Central", "East"), each = 2),
  Class = rep(c("Degraded zone", "Vegetated zone"), 3),
  Rate = c(
    16.2,  20.3,   # West
    32.0, -14.7,   # Central
    23.4,  -1.4    # East
  )
)

df$Region <- factor(
  df$Region,
  levels = c("West", "Central", "East")
)

df$Class <- factor(
  df$Class,
  levels = c("Degraded zone", "Vegetated zone")
)

# ============================================================
# 2. ORIGINAL COLORS
# ============================================================

class_colors <- c(
  "Degraded zone" = "#B08B5B",
  "Vegetated zone" = "#00441B"
)

region_colors <- c(
  "West" = "#F8C532",
  "Central" = "#F68D29",
  "East" = "#D92529"
)

# ============================================================
# 3. BAR POSITIONS AND VALUE LABELS
# ============================================================

df$x <- as.numeric(df$Region) +
  ifelse(df$Class == "Degraded zone", -0.19, 0.19)

df$label_y <- ifelse(
  df$Rate >= 0,
  df$Rate + 3.5,
  df$Rate - 4.0
)

df$rate_label <- sprintf("%+.1f", df$Rate)

# ============================================================
# 4. CORRECTED SALINITY LABELS
# ============================================================

region_df <- data.frame(
  Region = factor(
    c("West", "Central", "East"),
    levels = c("West", "Central", "East")
  ),
  x = 1:3,
  y = rep(-35, 3),
  label = c(
    "West\n17.9 PSU",
    "Central\n21.9 PSU",
    "East\n29.8 PSU"
  )
)

# ============================================================
# 5. CREATE FIGURE
# ============================================================

p <- ggplot() +
  
  # Horizontal gridlines
  geom_hline(
    yintercept = c(-20, 20, 40),
    color = "#E5E5E5",
    linewidth = 0.5
  ) +
  
  # Zero baseline
  geom_hline(
    yintercept = 0,
    color = "grey45",
    linewidth = 0.9
  ) +
  
  # Bars
  geom_col(
    data = df,
    aes(x = x, y = Rate, fill = Class),
    width = 0.32
  ) +
  
  # Numerical labels
  geom_text(
    data = df,
    aes(x = x, y = label_y, label = rate_label),
    size = 5.5,
    color = "black",
    family = "Times New Roman"
  ) +
  
  # Expansion annotation
  annotate(
    "text",
    x = 0.58,
    y = 39,
    label = "expansion \u2191",
    color = "black",
    fontface = "bold",
    size = 8,
    hjust = 0,
    family = "Times New Roman"
  ) +
  
  # Retraction annotation
  annotate(
    "text",
    x = 0.58,
    y = -15,
    label = "retraction \u2193",
    color = "black",
    fontface = "bold",
    size = 8,
    hjust = 0,
    family = "Times New Roman"
  ) +
  
  # Regional salinity labels
  geom_text(
    data = region_df,
    aes(x = x, y = y, label = label, color = Region),
    size = 6,
    fontface = "bold",
    lineheight = 1.15,
    family = "Times New Roman"
  ) +
  
  scale_fill_manual(values = class_colors) +
  
  scale_color_manual(
    values = region_colors,
    guide = "none"
  ) +
  
  scale_x_continuous(
    limits = c(0.5, 3.5),
    breaks = NULL,
    expand = expansion(mult = 0)
  ) +
  
  scale_y_continuous(
    breaks = c(-20, 0, 20, 40),
    limits = c(-48, 49),
    expand = expansion(mult = 0)
  ) +
  
  labs(
    x = NULL,
    y = expression(Rate~(m~yr^{-1})),
    fill = NULL
  ) +
  
  theme_minimal(
    base_size = 18,
    base_family = "Times New Roman"
  ) +
  
  theme(
    legend.position = "none",
    
    axis.title.y = element_text(
      size = 19,
      color = "black",
      margin = margin(r = 12)
    ),
    
    axis.text.y = element_text(
      size = 20,
      color = "black"
    ),
    
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    panel.grid = element_blank(),
    
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    
    plot.margin = margin(
      t = 15,
      r = 25,
      b = 20,
      l = 20
    )
  )

# ============================================================
# 6. DISPLAY FIGURE
# ============================================================

print(p)

# ============================================================
# 7. SAVE FIGURES
# ============================================================

output_dir <- "/Users/bz294/Documents/WhiteZone/Figure/Backup"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  file.path(
    output_dir,
    "SESE_Salinity_Change_Rates_Corrected.tif"
  ),
  plot = p,
  width = 11,
  height = 6.5,
  units = "in",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

ggsave(
  file.path(
    output_dir,
    "SESE_Salinity_Change_Rates_Corrected.png"
  ),
  plot = p,
  width = 11,
  height = 6.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

cat("Corrected figures saved to:", output_dir, "\n")


### FIGURE #3 ####

### FIGURE #3 ###

library(ggplot2)
library(grid)

# ============================================================
# 1. SCENARIO DATA
# ============================================================

df <- data.frame(
  Scenario = rep(
    c(
      "More freshwater",
      "Current trajectory",
      "Less freshwater"
    ),
    each = 2
  ),
  
  Class = rep(
    c("Vegetated zone", "Degraded zone"),
    times = 3
  ),
  
  Area = c(
    97.2, 54.0,    # More freshwater
    76.2, 67.0,    # Current trajectory
    63.0, 74.0     # Less freshwater
  )
)

# IMPORTANT:
# 97.2 and 76.2 are values reported in the manuscript.
# The remaining values (54.0, 67.0, 63.0, 74.0)
# are based on the original rounded chart values.
# Replace them with exact scenario outputs when available.


# ============================================================
# 2. BAR POSITIONS
# ============================================================

df$x <- c(
  0.8, 1.2,   # More freshwater
  2.1, 2.5,   # Current trajectory
  3.4, 3.8    # Less freshwater
)


# ============================================================
# 3. 2024 REFERENCE VALUES
# ============================================================

# Original approximate baseline values
# Replace with measured 2024 values when available

vz_2024 <- 62
dz_2024 <- 51


# ============================================================
# 4. COLORS
# ============================================================

class_colors <- c(
  "Vegetated zone" = "#00441B",
  "Degraded zone" = "#B08B5B"
)


# ============================================================
# 5. CREATE FIGURE
# ============================================================

p <- ggplot(df) +
  
  # ----------------------------------------------------------
# Bars
# ----------------------------------------------------------

geom_col(
  aes(
    x = x,
    y = Area,
    fill = Class
  ),
  width = 0.36
) +
  
  # ----------------------------------------------------------
# 2024 VZ reference line
# ----------------------------------------------------------

geom_hline(
  yintercept = vz_2024,
  color = "#00441B",
  linetype = "dashed",
  linewidth = 0.65
) +
  
  # ----------------------------------------------------------
# 2024 DZ reference line
# ----------------------------------------------------------

geom_hline(
  yintercept = dz_2024,
  color = "#B08B5B",
  linetype = "dashed",
  linewidth = 0.65
) +
  
  # ----------------------------------------------------------
# Bar values: ONE DECIMAL PLACE
# ----------------------------------------------------------

geom_text(
  aes(
    x = x,
    y = Area + 2,
    label = sprintf("%.1f", Area)
  ),
  size = 4,
  fontface = "bold",
  family = "Times New Roman",
  color = "black",
  vjust = 0
) +
  
  # ----------------------------------------------------------
# 2024 VZ label
# ----------------------------------------------------------

annotate(
  "text",
  x = 4.18,
  y = vz_2024 + 5,
  label = "2024 VZ",
  hjust = 0,
  color = "#00441B",
  size = 4,
  fontface = "bold",
  family = "Times New Roman"
) +
  
  # ----------------------------------------------------------
# 2024 DZ label
# ----------------------------------------------------------

annotate(
  "text",
  x = 4.18,
  y = dz_2024 - 5,
  label = "2024 DZ",
  hjust = 0,
  color = "#B08B5B",
  size = 4,
  fontface = "bold",
  family = "Times New Roman"
) +
  
  # ----------------------------------------------------------
# Class colors
# ----------------------------------------------------------

scale_fill_manual(
  values = class_colors,
  breaks = c(
    "Vegetated zone",
    "Degraded zone"
  )
) +
  
  # ----------------------------------------------------------
# X axis
# ----------------------------------------------------------

scale_x_continuous(
  breaks = c(1.0, 2.3, 3.6),
  labels = c(
    "More freshwater\n(low-salinity scenario)",
    "Current trajectory",
    "Less freshwater\n(high-salinity scenario)"
  ),
  expand = expansion(mult = 0)
) +
  
  # ----------------------------------------------------------
# Y axis
# ----------------------------------------------------------

scale_y_continuous(
  breaks = c(0, 25, 50, 75, 100),
  expand = expansion(mult = 0)
) +
  
  coord_cartesian(
    xlim = c(0.4, 5.0),
    ylim = c(0, 112),
    clip = "off"
  ) +
  
  # ----------------------------------------------------------
# Labels
# ----------------------------------------------------------

labs(
  title = "50-year land-cover scenarios (2024–2074)",
  x = NULL,
  y = expression("Area in 2074 (km"^2*")"),
  fill = NULL
) +
  
  # ----------------------------------------------------------
# TIMES NEW ROMAN FONT
# ----------------------------------------------------------

theme_minimal(
  base_size = 12,
  base_family = "Times New Roman"
) +
  
  theme(
    
    text = element_text(
      family = "Times New Roman"
    ),
    
    plot.title = element_text(
      size = 17,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    
    axis.text.x = element_text(
      size = 9,
      color = "black",
      lineheight = 1.15,
      margin = margin(t = 10)
    ),
    
    axis.text.y = element_text(
      size = 11,
      color = "black"
    ),
    
    axis.title.y = element_text(
      size = 12,
      color = "black",
      margin = margin(r = 10)
    ),
    
    axis.ticks = element_blank(),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    
    panel.grid.major.y = element_line(
      color = "#E8E8E8",
      linewidth = 0.4
    ),
    
    legend.position = "bottom",
    
    legend.text = element_text(
      size = 11,
      color = "black"
    ),
    
    legend.key.size = grid::unit(
      0.45, "cm"
    ),
    
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    
    plot.margin = margin(
      15, 15, 10, 15
    )
  )


# ============================================================
# 6. DISPLAY FIGURE
# ============================================================

print(p)


# ============================================================
# 7. SAVE FIGURES
# ============================================================

output_dir <- "/Users/bz294/Documents/WhiteZone/Figure/Backup"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# PNG
ggsave(
  filename = file.path(
    output_dir,
    "SESE_Future_Scenarios_2074_Vertical.png"
  ),
  plot = p,
  width = 7.5,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)


# TIFF
ggsave(
  filename = file.path(
    output_dir,
    "SESE_Future_Scenarios_2074_Vertical.tif"
  ),
  plot = p,
  width = 7.5,
  height = 6,
  units = "in",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)


cat(
  "\nFigures saved to:\n",
  output_dir,
  "\n"
)
