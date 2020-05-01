###Session 2 - Rilquer Mascarenhas

##Loading packages
library(raster)
library(rasterVis)
library(rgeos)
library(rgdal)
library(RColorBrewer)
library(viridis)
library("ggplot2")
library("sf")
library("ggspatial")
library("rnaturalearth")
library("rnaturalearthdata")

##Reading my data
data <- read.table('data/plants_occur.txt', sep = '\t', header = T)

##Reading climatic variables
#Creating study area extent
enclaves_extent <- extent(-41, -34,-8,-2)
#Loading tif files into stack and cropping to extent
bio<-crop(stack(list.files(path = 'data/wc0-5/', pattern = '\\.tif$', full.names = T)),enclaves_extent)
##Dividing by ten to get real Celsius values
bio <- bio/10
#Cropping the shapefile for Brazil states
br <- crop(readOGR(dsn = 'data/br/', layer = 'br'),enclaves_extent)

#Creating a vector with species name from the data loaded
spp <- levels(data[,1])

#Creating a color ramp
mycol=colorRampPalette(c('peachpuff2','darkgreen'))

##Plotting South America
class(world)
theme_set(theme_bw())
ggplot(data = world) +
  geom_sf(fill= "antiquewhite")+
  coord_sf(xlim = c(-91.531219, -28.23167), ylim = c(-56.791423,17.911542), expand = FALSE)+
  theme(panel.grid.major = element_line(color = gray(.9), linetype = "dashed", size = 0),
        panel.background = element_rect(fill = "aliceblue"))+
  ggsave('output/south_america.tiff')

##Using a loop to plot the points for each species on bio1 variable background
for (i in spp) {
  sp.data <- data[which(data[,1]==i),]
  gplot(bio[[1]]) +
    geom_tile(aes(fill=factor(value)))+
    scale_fill_manual(values = mycol(85))+
    theme(legend.position = c(0.2, 0.8))+
    annotation_scale(location = "br", width_hint = 0.5,pad_x = unit(0.5, "cm"), pad_y = unit(13, "cm"))+
    annotation_north_arrow(location = "br", which_north = "true", 
                           pad_x = unit(0.75, "in"), pad_y = unit(4, "in"),
                           style = north_arrow_fancy_orienteering)+
    coord_sf(xlim = c(-41,-34), ylim = c(-8,-2), expand = FALSE, crs = "+proj=longlat +datum=WGS84 +no_defs")+
    geom_polygon(data = br, aes(x = long, y = lat, group = group), colour = "black",
                 fill = NA, size = 0.1)+
    geom_point(data = sp.data, aes(x = Longitude,
                                   y = Latitude, color = Plant.species), size = 2)+
    scale_color_manual(name = "Species", values = 'red')+
    scale_x_discrete(name = "Longitude")+
    scale_y_discrete(name = "Latitude")+
    theme(panel.grid.major = element_line(color = gray(.9), linetype = "dashed", size = 0),
          panel.background = element_rect(fill = "aliceblue"), legend.position = 'none')+
    ggtitle(paste0('Occurences points for ',i))+
    ggsave(paste0('output/maps/',i,'.tiff'))
}
