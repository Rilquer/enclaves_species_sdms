install.packages('dismo')
install.packages('rgdal')
install.packages('rJava')
library(rJava)
library(dismo)
library(rgdal)

## We can fairly easily create a Maxent model in R. However, to do this you must:
# 1. Download Maxent from http://www.cs.princeton.edu/~schapire/maxent/
# 2. Copy the file Maxent.jar and paste a copy of it in ".../R(version)/library/dismo/java". 
#    There should already be "dismo.jar" in that folder. Add Maxent..ar to it

# Create a rasterStack of your bioclim data
enclaves_extent <- extent(-41, -34,-8,-2)
samer.extent <- extent(-89.755388,-29.110858,-55.978471,13.231814)
env <- crop(stack(list.files(path = 'data/wc0-5/', pattern = '\\.tif$', full.names = T)),enclaves_extent)
env.lgm <- crop(stack(list.files(path = 'data/lgm30sec', pattern = '\\.tif$', full.names = T)),enclaves_extent)
# Load locality information as two columns of longitude, latitude.
spp <- c('Cassia ferruginea','Myrcia silvatica', 'Scinax nebulosus','Sclerurus cearensis')
locs <- list()
locs[[4]] <- read.csv("data/occs/s_cearensis.csv")
locs[[2]] <- read.csv("data/occs/cassia_ferruginea.csv")
locs[[3]] <- read.csv("data/occs/myrcia_silvatica.csv")
locs[[1]] <- read.csv("data/occs/s_nebulosus.csv")


###Create that optimal model using all localities (no withheld data)
###x = rasterstack, p = locs object as data.frame, a = background coords, factors = categorical variables, 
### make true only the arguments wanted. E.g. using autofeatures lets the algorithm choose features based on the
# number of occurrences. Otherwise, features are used by default. You can use a feature by changing it to 'true'.
dir.create('output/maxent_models/')
mod <- list()
for (i in 1:length(spp)) {
  dir.create(paste0('output/maxent_models/',spp[i]))
  mod[[i]] <- maxent(
    x=env, # bio stack
    p=locs[[i]], # locality csv
    factors = NULL,
    path = paste0('output/maxent_models/',spp[i]),
    args=c(
      'betamultiplier=1',
      #'linear=false',
      #'quadratic=false',
      #'product=false',
      #'threshold=false',
      #'hinge=false',
      'threads=2',
      'responsecurves=true',
      'jackknife=true',
      'askoverwrite=false',
      'autofeature=true'
    )
  )
}
dir.create('output/maxent_models/current/')
pred.mod <- list()
for (i in 1:length(spp)) {
  pred.mod[[i]] <- predict(
    object = mod[[i]],
    x = env,
    filename = paste0('output/maxent_models/current/',spp[i],'.asc'),
    na.rm = T,
    format = 'ascii',#or GTiff
    overwrite = F,
    args = "logistic"
  )
}

dir.create('output/maxent_models/LGM')
pred.lgm <- list()
for (i in 1:length(spp)) {
  pred.lgm[[i]] <- predict(
    object = mod[[i]],
    x = env.lgm,
    filename = paste0('output/maxent_models/LGM/',spp[i],'.asc'),
    na.rm = T,
    format = 'ascii',#or GTiff
    overwrite = F,
    args = "logistic"
  )
}

library(RColorBrewer)
library(maps)
##Reading results and making maps
maps.current <- list()
files <- list.files(path = 'output/maxent_models/current/', pattern = '.asc', full.names = T)
for (i in files) {
  maps.current[[i]] <- raster(i)
}

maps.lgm <- list()
files <- list.files(path = 'output/maxent_models/LGM/', pattern = '.asc', full.names = T)
for (i in files) {
  maps.lgm[[i]] <- raster(i)
}

for (i in 1:length(maps.current)) {
  par(bty='n')
  plot(maps.current[[i]], col=brewer.pal(9,'GnBu'), legend=FALSE, axes=FALSE)
  r.range <- c(min(na.omit(values(maps.current[[i]]))), max(na.omit(values(maps.current[[i]]))))
  plot(maps.current[[i]],legend.only=TRUE,col=brewer.pal(9,'GnBu'),legend.width=1,
       axis.args=list(at=c(as.integer(seq(r.range[1], r.range[2], 0.1))),
                      labels=round(seq(r.range[1], r.range[2], 0.1),2)),
       legend.args=list(text='Suitability', side=4, font=2, line=-2, cex=0.8))
  map.scale(x=-37, y=-3, ratio=FALSE, relwidth=0.25)
  title(paste0(spp[i],' - Current'))
  dev.copy(tiff,filename=paste0('output/maxent_models/current/',spp[i],'.tif'),
           width = 6, height = 5, units = "in", res = 500)
  dev.off()
}

for (i in 1:length(maps.lgm)) {
  par(bty='n')
  plot(maps.lgm[[i]], col=brewer.pal(9,'GnBu'), legend=FALSE, axes=FALSE)
  r.range <- c(min(na.omit(values(maps.lgm[[i]]))), max(na.omit(values(maps.lgm[[i]]))))
  plot(maps.lgm[[i]],legend.only=TRUE,col=brewer.pal(9,'GnBu'),legend.width=1,
       axis.args=list(at=c(as.integer(seq(r.range[1], r.range[2], 0.1))),
                      labels=round(seq(r.range[1], r.range[2], 0.1),2)),
       legend.args=list(text='Suitability', side=4, font=2, line=-2, cex=0.8))
  map.scale(x=-37, y=-3, ratio=FALSE, relwidth=0.25)
  title(paste0(spp[i],' - LGM'))
  dev.copy(tiff,filename=paste0('output/maxent_models/LGM/',spp[i],'.tif'),
           width = 6, height = 5, units = "in", res = 500)
  dev.off()
}





