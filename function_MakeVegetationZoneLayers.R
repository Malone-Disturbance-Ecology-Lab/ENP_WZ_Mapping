
# 
##' .. his function is designed to turn a rast stack into polygons. ..
##'
##' @title polygonize_SpatRast, this function will produce a simplie feature.
##' @param stack raster stack of vegetation zones 1, 2, or 3
##' @return sfc a simple feature of the polygons
##' @author Sparkle L. Malone
##' 
library(terra)
library(tidyverse)
library(sf)

polygonize_SpatRast <- function( stack){
  for ( i in 1: nlyr(stack)){
    print(i)
    rast <- stack[[i]]
    names(rast) <- 'layers'
    wz.vector <- as.polygons(rast, dissolve=T) 
    
    wz.vector %>% plot
    
    wz.vector$time <- time(rast)
    
    if(  i == 1 ){
      simpleF <- wz.vector
    } else{
      simpleF <- rbind( simpleF, wz.vector)
    }
    
  }  
  return(simpleF %>% st_as_sf)
}

##' @title Trans_Intersection, this function will produce a simplie feature.
##' @param transects Transect layer produces with makeTransects
##' @param polygons  Polygon layer produces with polygonize_SpatRast
##' @return sfc a simple feature of the transects
##' @author Sparkle L. Malone
Trans_Intersection <- function( transects, polygons ){
  
  Intersect <- st_intersection(transects, polygons)
  Intersect$length.new <- st_length(Intersect)
  Intersect$length.Delta <- (Intersect$length.new - Intersect$length.org)/Intersect$length.org *100
  
  print(Intersect$length.Delta %>% hist )
  
  return( Intersect)
}