
library(rJava)
options(java.parameters = "-Xmx1g")
library(raster)
library(rgdal)
library(rgeos)
library(spocc)
library(spThin)
library(dismo)
library(ENMeval)
library(dplyr)
library(wallace)
library(maps)

#Also loading some necessary functions from Wallace package:
source(system.file('shiny/funcs', 'functions.R', package = 'wallace'))

# Now loading occurence data:

spp <- c('Cassia ferruginea','Myrcia silvatica', 'Scinax nebulosus','Sclerurus cearensis')
colors <- c('green','orange','purple','blue')
occs.xy <- list()
occs.xy[[1]] <- read.csv("data/occs/cassia_ferruginea.csv")
occs.xy[[2]] <- read.csv("data/occs/myrcia_silvatica.csv")
occs.xy[[3]] <- read.csv("data/occs/s_nebulosus.csv")
occs.xy[[4]] <- read.csv("data/occs/s_cearensis.csv")
for (i in 1:length(occs.xy)) {
  colnames(occs.xy[[i]]) <- c('longitude','latitude')
  sp::coordinates(occs.xy[[i]]) <- ~ longitude + latitude
  occs.xy[[i]] <- crop(occs.xy[[i]],enclaves_extent)
}

# Reading environmental data, BR shapefiles and setting other geographical boundaries:

enclaves_extent <- extent(-41, -34,-8,-2)
br <- crop(readOGR(dsn = 'data/br/', layer = 'br'),enclaves_extent)
env <- crop(stack(list.files(path = 'data/wc0-5/', pattern = '\\.tif$', full.names = T)),
            enclaves_extent)

## Getting background points through minimum convex polygon:

# Creating empty lists
bgExt <- list()
envsBgCrop <- list()
envsBgMsk <- list()
bg.xy <- list()

for (i in 1:length(occs.xy)) {
  bgExt[[i]] <- mcp(occs.xy[[i]])
  bgExt[[i]] <- rgeos::gBuffer(bgExt[[i]], width = 0.5)
  
  # crop the environmental rasters by the background extent shape
  envsBgCrop[[i]] <- raster::crop(env, bgExt[[i]])
  # mask the background extent shape from the cropped raster
  envsBgMsk[[i]] <- raster::mask(envsBgCrop[[i]], bgExt[[i]])
  # sample random background points
  bg.xy[[i]] <- dismo::randomPoints(envsBgMsk[[i]], 1000)
  # convert matrix output to data frame
  bg.xy[[i]] <- as.data.frame(bg.xy[[i]])
  plot(br)
  points(bg.xy[[i]], col = 'black', pch = 20)
  points(occs.xy[[i]], col = colors[i], pch = 20)
  title(paste0(spp[i],' - Presence and Background'))
  dev.copy(tiff,filename=paste0('output/',spp[i],' - Presence and Background.tif'),
           width = 6, height = 5, units = "in", res = 500)
  dev.off()
}

#Running bioclim:
presvals <- list()
bc <- list()
bc.projections <- list()
for (i in 1:length(occs.xy)) {
  message('Extracting environmental values for presence points of ',spp[i])
  presvals[[i]] <- extract(env, occs.xy[[i]])
  message('Running bioclim for ',spp[i])
  bc[[i]] <- bioclim(presvals[[i]])
  bc.projections[[i]] <- predict(env, bc[[i]])
}

# Plotting and saving bioclim models:
dir.create('output/bioclim_models')
dir.create('output/bioclim_models/current')
for (i in 1:length(bc.projections)) {
  par(bty='n')
  r.range <- c(min(na.omit(values(bc.projections[[i]]))),max(na.omit(values(bc.projections[[i]]))))
  plot(bc.projections[[i]], axes=FALSE,legend.width=1,
       axis.args=list(at=round(seq(r.range[1], r.range[2], 0.1),2),
                      labels=round(seq(r.range[1], r.range[2], 0.1),2)),
       legend.args=list(text='Suitability', side=4, font=2, line=-2, cex=0.8))
  map.scale(x=-37, y=-3, ratio=FALSE, relwidth=0.25)
  title(paste0(spp[i],' - Bioclim Current'))
  dev.copy(tiff,filename=paste0('output/bioclim_models/current/',spp[i],'.tif'),
           width = 6, height = 5, units = "in", res = 500)
  dev.off()
}

## Running and evaluating Maxent models:
# Selecting spatial partition scheme
group.data <- list()
occs.grp <- list()
bg.grp <- list()
for (i in 1:length(occs.xy)) {
  group.data[[i]] <- ENMeval::get.checkerboard1(occ=occs.xy[[i]], env=envsBgMsk[[i]],
                                                bg.coords=bg.xy[[i]],aggregation.factor=2)
  # pull out the occurrence and background partition group numbers from the list
  occs.grp[[i]] <- group.data[[i]][[1]]
  bg.grp[[i]] <- group.data[[i]][[2]]
}

# define the vector of regularization multipliers to test
rms <- seq(0.5, 3.5, 1)
# iterate model building over all chosen parameter settings
e <- list()
for (i in 1:length(occs.xy)) {
  e[[i]] <- ENMeval::ENMevaluate(occs.xy[[i]], envsBgMsk[[i]], bg.coords = bg.xy[[i]],
                                 RMvalues = rms, fc = c('L', 'LQ', 'H', 'LQH', 'LQHP'),
                                 method = 'user', occs.grp[[i]], bg.grp[[i]], clamp = TRUE,
                                 algorithm = "maxnet")
  # unpack the results data frame, the list of models, and the RasterStack of raw predictions
  evalTbl[[i]] <- e[[i]]@results
  evalMods[[i]] <- e[[i]]@models
  names(evalMods[[i]]) <- e[[i]]@results$settings
  evalPreds[[i]] <- e[[i]]@predictions
}
