##' .. This function is designed to make a transect layer to measure landscape change. The function requires the following parms: ..
##'
##' @title makeTransects, this function will produce transects that are perpendicular to two reference lines.
##' @param line1 is the first lie
##' @param line2 is the second line
##' @param transect.n is the number of transects you want
##' @return sfc a simple feature of the transects and their lengths (length.org; Meters)
##' @author Sparkle L. Malone

library(sf)
library(sfheaders)
library(tidyverse)
library(ggplot2)

coast.sf <- read_sf('/Volumes/Malonelab/Research/Whitezone/Spatial_Data/Coastal.kml')
inland.sf <- read_sf('/Volumes/Malonelab/Research/Whitezone/Spatial_Data/Inland.kml')
aoi <- read_sf('/Volumes/MaloneLab/Research/ENP/The Whitezone/WhiteZoneMap/SamplingGrid/AOI.kml')

plot(coast.sf$geometry)
plot(inland.sf$geometry, add=T)
plot(aoi$geometry, add=T)

makeTransect <- function(line1, line2, transect.n){

line1.points <- st_sample(line1, size=transect.n, type='regular') %>% st_cast( "POINT") %>% st_as_sf()
line2.points <- st_sample(line2, size=transect.n, type='regular') %>% st_cast( "POINT") %>% st_as_sf()

line1.points <- line1.points[1:transect.n,]
line2.points <- line2.points[1:transect.n,]

line1.points$ID <- seq(1, transect.n,1)
line1.points$Type <- "line1"
line2.points$ID <- seq(1, transect.n,1)
line2.points$Type <- "line2"

## Put all points into a single file
aoi.points <- rbind( line1.points,line2.points )
coordinates <- st_coordinates(aoi.points)  
aoi.points$latitude <- coordinates[, 1] 
aoi.points$longitude <- coordinates[, 2] 

# reformat the point file: 
#m <- aoi.points %>% as.data.frame %>% select("latitude", "longitude", "ID")
m <- as.data.frame(aoi.points)[, c("latitude", "longitude", "ID")]
m2 <-data.matrix(m) # Make the file a matrix of lat, long, and ID

## set the order by the 'ID' column
m2 <- m2[ order( m2[,3] ), ]

## then use `sfheaders` to create your linestrings
sfc <- sfheaders::sfc_linestring(obj = m2, linestring_id = 3) %>% st_as_sf() ## 3rd column)
st_crs(sfc) <- st_crs(line1)
sfc$length.org <- st_length(sfc) # add the length of each transect here:

# Show a plot of the transects
sfc$ID <- seq(1, transect.n, 1)

final.plot <- ggplot( ) + geom_sf( data = line1 ) + geom_sf( data = line2 ) + geom_sf( data = sfc)
print( final.plot)
 # once you get the length I cant plot anymore ...

return(sfc)

}

ENP.transects <- makeTransect(line1 = coast.sf, line2 =  inland.sf, transect.n= 100)

