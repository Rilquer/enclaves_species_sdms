##Lab 6 - Rilquer Mascarenhas

library(raster)
library(rgeos)
library(rgdal)
library(sf)
library(ggplot2)

##Loading data
treelop.points <- read.csv('data/Treelop_data/Treelop_sightings.csv')
usa <- st_read('data/Treelop_data/USA_states/USA.shp')
landcover <- st_read('data/Treelop_data/Landcover/lulc250k.shp')
st_crs(landcover) = "+proj=longlat +datum=WGS84 +no_defs" #Converting crs of landcover to match other layers
roads <- st_read('data/Treelop_data/Roads/ne_10m_roads.shp')

#Subsetting usa states to NY state only, and savin
NewYorkState <- usa[which(usa$name=='New York'),]
dir.create('data/Treelop_data/NY_state')
st_write(NewYorkState,dsn = 'data/Treelop_data/NY_state/',layer = 'NewYorkState', driver = 'ESRI Shapefile')

##Intersecting other layers to reduce to NY state area (st_intersection for vector)
##and mask for raster
NY_landcover <- st_intersection(landcover,NewYorkState)
NY_roads <- st_intersection(roads,NewYorkState)
bio1 <- raster('data/Treelop_data/NY_regional_bio_1.tif')
NY_bio_1 <- mask(bio1,NewYorkState)
writeRaster(NY_bio_1,'data/Treelop_data/NY_bio_1.tif')

##Creating suitable raster from observed occurences
NY_bio_1_suitable <- NY_bio_1
low <- round(min(extract(NY_bio_1,SpatialPoints(treelop.points[,2:3]))),0)
high <- round(max(extract(NY_bio_1,SpatialPoints(treelop.points[,2:3]))),0)
NY_bio_1_suitable <- (NY_bio_1_suitable > low)*NY_bio_1_suitable & (NY_bio_1_suitable < high)*NY_bio_1_suitable 

##Transforming suitable raster into polygon for further analyses
NY_bio_1_suitable_polygon <- rasterToPolygons(NY_bio_1_suitable, dissolve = TRUE)
NY_bio_1_suitable_polygon <- NY_bio_1_suitable_polygon[which(NY_bio_1_suitable_polygon$layer==1),]

##Creating vector of preferred habitats and subsetting landcover shapefile through
##LUTEXT based on that vector.
prefered_habitats <- c("Deciduous forest land","Evergreen forest land","Forested wetland","Mixed forest land","Shrub and brush rangeland")
Treelop_preferable_landcover <- landcover[which(landcover$LUTEXT %in% prefered_habitats),]

##Creating buffer around roads with st_buffer, using appr. 3.33 km (0.03 in arc degress)
NY_roads_buffered <- st_buffer(NY_roads, dist = 0.03)

##Transforming sf objects into spatial objects.
NY_roads_buffered_sp <- as_Spatial(NY_roads_buffered)
Treelop_preferable_landcover_sp <- as_Spatial(Treelop_preferable_landcover)
NY_bio_1_suitable_polygon_sp <- as_Spatial(NY_bio_1_suitable_polygon)

##Using gIntersection and gDifference to overlay the three datasets and look
##for suitable areas with appropriate landcover and far enough from major roads
LC_temp_intersection <- gIntersection(Treelop_preferable_landcover_sp,NY_bio_1_suitable_polygon_sp)
LC_temp_roads_difference <- gDifference(LC_temp_intersection@polyobj,NY_roads_buffered_sp)
LC_temp_roads_difference <- st_as_sf(LC_temp_roads_difference)

#Plotting final map
library("ggspatial")
library("rnaturalearth")
library("rnaturalearthdata")

ggplot(data = usa) +
  geom_sf(fill= "antiquewhite")+
  coord_sf(xlim = c(-127.734030, -63.685048), ylim = c(21.730400,48.820799), expand = FALSE)+
  theme(panel.grid.major = element_line(color = gray(.9), linetype = "dashed", size = 0),
        panel.background = element_rect(fill = "aliceblue"))

ggplot(data = LC_temp_roads_difference) +
  geom_sf(fill= "red", size = 0.1)+
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-80.178559,-73.010285), ylim = c(40.486597,45.207088), expand = FALSE)+
  geom_sf(data = NewYorkState, colour = "black", fill = NA, size = 0.1)+
  geom_point(data = treelop.points, aes(x = Longitude,
                              y = Latitude), size = 1, color = 'green')+
  scale_x_discrete(name = "Longitude")+
  scale_y_discrete(name = "Latitude")+
  theme(panel.grid.major = element_line(color = gray(.9), linetype = "dashed", size = 0),
        panel.background = element_rect(fill = "aliceblue"), legend.position = "bottom")+
  ggsave('output/lab6.pdf')
