require(dismo)
require(ENMeval)
require(phyloclim)
require(sp)
require(rgdal)
require(rgeos)
#################################################################################################################
## ENMTools package
# https://github.com/danlwarren/ENMTools
install.packages("devtools")
library(devtools)
devtools::install_github("danlwarren/ENMTools")
library(ENMTools)

#################################################################################################################
#### Niche comparisons using Schoener's D. 
# Load two rasters to compare
r1<-raster('Pathway/to/first/raster.tif') # can be any raster type (i.e. .asci, .tif, .bil)
r2<-raster('Pathway/to/second/raster.tif')
flav<-raster("E_m_flavifrons.asc")
mac<-raster("E_m_macaco.asc")

# Calculate niche overlap between the rasters. Change stat to 'I' for Hellinger's I test, and 
# to 'D' for Schoener's D. Remember to use ?nicheOverlap to see the arguments for this function. 
Overlap <- nicheOverlap(flav, mac, stat='D', mask=TRUE, checkNegatives=TRUE)
# This function is used by ENMeval to calculate the niche overlap between candidate models during the
# model tuning process. Set overlap to TRUE to get pairwise Schoener's D values for all candidate models. 

## ENMTools
ENMOverlap<-raster.overlap(flav, mac)
#################################################################################################################
#### Maxent model tuning using ENMeval
# Load environmental rasters as stack. Change pattern to .asc for ascii.
Env<-stack(list.files(path = "Pathway/to/environmental/data", pattern = '\\.tif$', full.names = T))
# Load locality information with 3 columns of "Species", "Longitude", "Latitude"
locs<-read.csv('Pathway/to/locality/data.csv')

# Here, method refers to data partitioning method. Categoricals refers to the names of any environmental data
# that are categorical. fc refers to feature classes to use. bg.coords refers to the background points used every time. You
# can set this as an object (i.e. a .csv file). RMvalues refers to the regularization multipliers to use. 
# Remember to check out ?ENMevaluate
res <- ENMevaluate(occ = locs, env = Env, method='block', categoricals=NULL, 
                   fc = c("L", 'LQ', "H", "LQH"), 
                   bg.coords=NULL, RMvalues=seq(1, 5, 0.5), overlap = T)

###################################################################################################################
#### Identity Test (same as Warren 2008, ENMtools)
setwd('data/Lemurdata/')
# Load environmental variables
env<-stack(list.files(path='Lemur_layers', pattern = '\\.asc$', full.names=T))
# Load occurrence records for both species
setwd('C:/Users/pgalante/Documents/Projects/QGIS_tutorial/RGGS_GIS/Session6_data') # These should be csv files of records where columns are: "Species, X, Y".
flav <- read.csv("e_m_flav.csv")
maca <- read.csv("e_m_maca.csv")
# Change the species columns to just the species' names
flav[1] <- as.factor('flavifrons')
maca[1] <- as.factor('macaco')
# row bind them so all occurrences are in 3 rows of Species, X, Y
sites<-rbind(maca, flav)
species <- c('flavifrons','macaco')
# Change the column names of sites
colnames(sites)<-c("species","longitude","latitude")
samples <- sites[grep(paste(species, collapse = "|"), sites$species), ]
# Tell R where maxent is (the copy that in with dismo).
maxent.exe <- paste(system.file(package="dismo"),"/java/maxent.jar", sep = "")
### ?niche.equivalency.test
nicheEquivalency<-niche.equivalency.test(p = samples, env = env, app=maxent.exe, dir = 'NicheEquivalence')
# Note that you can also perform the background test using bg.similarity.test. For more see ?bg.similarity.test

####Using ENMTools
# Load environmental variables
env<-stack(list.files(path='C:/Users/pgalante/Documents/Projects/QGIS_tutorial/RGGS_GIS/Session6_data/Lemur_layers', pattern = '\\.asc$', full.names=T))
setwd('C:/Users/pgalante/Documents/Projects/QGIS_tutorial/RGGS_GIS/Session6_data')
## We need all occurrence data and environmental layers to be in decimal degrees, which we can do
# First read in csv data
maca<-read.csv('e_m_maca.csv')[,2:3]
flavi<-read.csv('e_m_flav.csv')[,2:3]
# Now, convert into a SpatialPoints object with the proj4string for the area of interest- here for Madagascar
macasp<-SpatialPoints(maca, proj4string = CRS('+proj=utm +zone=38 +south +ellps=intl +towgs84=-189,-242,-91,0,0,0,0 +units=m +no_defs'))
flavisp<-SpatialPoints(flavi, proj4string = CRS('+proj=utm +zone=38 +south +ellps=intl +towgs84=-189,-242,-91,0,0,0,0 +units=m +no_defs'))
# Convert to lat long for decimal degrees
macadd<-spTransform(macasp, crs('+proj=longlat +datum=WGS84 +no_defs'))
flavidd<-spTransform(flavisp, crs('+proj=longlat +datum=WGS84 +no_defs'))
# Convert back to a dataframe and set the column names
maca.df<-as.data.frame(macadd)
flavi.df<-as.data.frame(flavidd)
colnames(maca.df)<-c('long','lat')
colnames(flavi.df)<-c('long','lat')
# Now we need to change the coordinate reference system of the environmental layers. Set the CRS of this, as before
crs(env)<-'+proj=utm +zone=38 +south +ellps=intl +towgs84=-189,-242,-91,0,0,0,0 +units=m +no_defs'
# Then reproject it to the desired CRS
env2<-projectRaster(env, crs='+proj=longlat +datum=WGS84 +no_defs')
env2<-stack(env2)

env<-stack(list.files(path='Lemur_layers', pattern = '\\.asc$', full.names = T))
flav<-enmtools.species()
flav$species.name <- "flav"
flav$presence.points <- flavi.df
flav$range <- background.raster.buffer(flav$presence.points, 50000, mask = env2)
flav$background.points <- background.points.buffer(points = flav$presence.points,
                                                   radius = 20000, n = 1000, mask = env2[[1]])
mac<-enmtools.species()
mac$species.name <- "mac"
mac$presence.points <- maca.df
mac$range <- background.raster.buffer(mac$presence.points, 50000, mask = env2)
mac$background.points <- background.points.buffer(points = mac$presence.points,
                                                  radius = 20000, n = 1000, mask = env2[[1]])
colnames(flav$presence.points)<- colnames(flav$background.points) <- c("Longitude","Latitude")
colnames(mac$presence.points)<- colnames(mac$background.points) <- c("Longitude","Latitude")

id.glm <- identity.test(species.1 = flav, species.2 = mac, env = env2, type = "glm", nreps = 100)
library(ggplot2)
rep.values <- data.frame(id.glm$reps.overlap[2:nrow(id.glm$reps.overlap),1:2])
rep.values <- data.frame(nicheEquivalency$null.distribution)
ggplot(data = rep.values, aes(x = D)) + geom_histogram(color = 'black',
                                                       fill = 'white', binwidth = 0.11)+
  geom_vline(xintercept = 0.1867669, color="red", linetype="dashed", size=0.5)+
  ggtitle("Null distributation and observed value of Schoner's D (Identity test)")+
  ggsave('output/identity_test_lemur.pdf')
