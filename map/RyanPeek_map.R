# Libraries ---------------------------------------------------------------

library(dplyr)
library(sf)  
library(ggplot2)
library(ggthemes)
library(tmap)    # for static and interactive maps
library(viridis)
#library(spData) # contains datasets used in this example
#library(spDataLarge)
library(grid)
library(USAboundaries)
library(ggspatial)


# download package to interact with natural earth data
#devtools::install_github("ropenscilabs/rnaturalearth")
#devtools::install_github("ropenscilabs/rnaturalearthdata")
#install.packages("rnaturalearthhires",
#                 repos = "http://packages.ropensci.org",
#                 type = "source")
library(rnaturalearth)


# Base US Map -------------------------------------------------------------

us_states<- us_states() %>% filter(!stusps=="HI", !stusps=="AK", !stusps=="PR") #%>% st_transform(2163)
us_states_map = tm_shape(us_states, projection = 2163) + tm_polygons() + 
  tm_layout(frame = FALSE)
us_states_map

# GET DATA ----------------------------------------------------------------

#rivers
hucs_sf <- ne_download(scale = 10, type = 'rivers_lake_centerlines', category = 'physical', returnclass = "sf")
#plot(hucs_sf$geometry, col="blue2")
st_crs(hucs_sf)
us_states_map + tm_shape(hucs_sf) + tm_lines(col="blue")

# import fish sampling points
locations <- read.csv("https://raw.githubusercontent.com/johnsolk/RNAseq_15killifish/master/map/FishSamplingLOCATIONS.csv",stringsAsFactors = FALSE) %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326, remove=F) # make spatial as SF

st_crs(locations)

# make a map ------------------------------------------------------------
cols <- c("FW" = "#FF3300", "M" = "#0000FF")
ggplot() + 
  geom_sf(data= us_states, fill="gray50", alpha=0.5) + 
  geom_point(data=locations, aes(x=Long, y=Lat, 
                                 fill=Native.Physiology), pch=21, size=6) + 
  theme_bw(base_size = 12) +
  annotation_scale(location = "bl",style = "bar",pad_y=unit(0.2, "cm")) +
  coord_sf() + 
  theme(plot.background = element_blank(),
        legend.position = c(0.1, 0.2),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        legend.key = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(c(0, 0, 0 ,0), "mm")) +
  scale_colour_manual(values = cols,aesthetics = c("colour", "fill")) +
  ggrepel::geom_label_repel(data=locations, size=5,aes(x = Long, y=Lat, label=Species), segment.color = "black",
                            cex=3, segment.alpha = 1,
                            force = 3,  box.padding=1, nudge_y = 0.3) 

ggsave(filename = "figs/map_example_killifish.png", dpi=300, width = 8, height = 6, units="in")
