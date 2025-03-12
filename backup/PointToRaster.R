library(sf)
library(raster)
library(terra)
##library(gdalUtils)

shapefile_path <- "/Users/bz294/Documents/WhiteZone/Code/output_points_2013.shp"
points <- st_read(shapefile_path)

resolution <- 30  

# Create a blank raster with the same extent as the shapefile
r <- raster(extent(st_bbox(points)), resolution = resolution, crs = st_crs(points)$proj4string)

# Rasterize the points based on 'y_region'
r <- rasterize(as(points, "Spatial"), r, field = "y_region", fun = mean)  # Adjust 'fun' as needed
#r_terra <- rast(r)
#r_filled <- terra::focal(r_terra, w = matrix(1, 3, 3), fun = "mean", na.policy = "omit", na.rm = TRUE)

#polygon_mask <- st_read("/Users/bz294/Documents/WhiteZone/Points/Shapefiles/SEregion.shp")

# Mask the raster to exclude areas outside the polygon
#r_filled <- mask(r_filled, vect(polygon_mask))

output_path <- "/Users/bz294/Documents/WhiteZone/Code/raster_2013.tif"
writeRaster(r, output_path, overwrite = TRUE)
