library(dplyr)
library(sf)  
library(sp)
library(rgdal)
library(raster)
library(ggplot2)
library(tmap)    # for static and interactive maps
library(viridis)
library(spData) # contains datasets used in this example
library(spDataLarge)
library(tmap)    # for static and interactive maps
library(grid)
library(USAboundaries)
# test, should work
# loading library order matters
# might need to install.packages("tmap") 
# then load libraries again
# no idea why this happens
us_states_map = tm_shape(us_states, projection = 2163) + tm_polygons() + 
  tm_layout(frame = FALSE)
us_states_map

setwd("~/Documents/UCDavis/Whitehead/map/")
shapefile_dir <- "~/Documents/UCDavis/Whitehead/map/10m_physical/"

# tutorial from Ryan
# https://ryanpeek.github.io/mapping-in-R-workshop/vig_workflow_in_R_snowdata.html#
ogrInfo(dsn="./data", layer="ne_10m_rivers_lake_centerlines")
hucs_sp<- readOGR(dsn = "./data", layer = "ne_10m_rivers_lake_centerlines") # takes a few seconds!
proj4string(hucs_sp) # check projection, should be WGS84
hucs_sf <- st_read("data/ne_10m_rivers_lake_centerlines.shp") # fast!
st_crs(hucs_sp)
extent(hucs_sp)
plot(st_geometry(hucs_sf), col="darkblue", axes=TRUE)

# import fish sampling points
# tutorial:
# https://www.neonscience.org/dc-csv-to-shapefile-r
locations <- read.csv("FishSamplingLOCATIONS.csv",stringsAsFactors = FALSE)
str(locations)
head(locations$Lat)
head(locations$Long)

# create custom crs
utm18nCRS <- crs(hucs_sp)
utm18nCRS
class(utm18nCRS)
plotLocations <- SpatialPointsDataFrame(locations[,c(5,4)],
                                                locations,    #the R object to convert
                                                proj4string = utm18nCRS)   # assign a CRS 
crs(plotLocations)
plot(plotLocations, 
     main="Fish Sampling Locations", axes = TRUE)
extent(hucs_sf)
plotLocations.extent<-extent(plotLocations)
xmin <- plotLocations.extent@xmin
xmax <- plotLocations.extent@xmax
ymin <- plotLocations.extent@ymin
ymax <- plotLocations.extent@ymax
# can't get sampling locations to plot on top?
plot(hucs_sp,
     main="Fish Sampling Locations",
     xlim=c(xmin,xmax),
     ylim=c(ymin,ymax))
tmap_mode("plot")
colors = c("red","blue")
shapes = c(19,17)
tm_shape(us_states, projection = 2163) +
  tm_polygons() +
  tm_shape(plotLocations) +
  tm_symbols(col="Native.Physiology",palette=colors,size=2,title.size=4) +
  tm_layout(frame = FALSE) +
  tm_text("Species",ymod=1,auto.placement=TRUE) +
  tm_scale_bar(position = c("LEFT","BOTTOM"),text.size = 1)


