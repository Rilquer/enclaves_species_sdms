library(ENMeval)

#####Set Working Directory#####
setwd("C:/Users/renat/OneDrive/Desktop/me")
      
#####Loading Data#####
#ALL OCCURRENCES
locs<-read.csv("euploea.csv") #load occurence points
locs<-locs[,2:3] #keep columns with lon/lat
locs[locs==""]<-NA #turn blanks into NA
colnames(locs) <- c("x","y") #rename columns
locs<- subset(locs, !is.na(x)) #remove rows with NAs
locs<- subset(locs, !is.na(y)) #remove rows with NAs

env<-stack((list.files(path="C:/Users/renat/OneDrive/Desktop/me/wc2.1_10m_bio",pattern = '.tif', full.names = T)[c(17,10,11)]))
env<-crop(env, extent(locs)) #Crop Layers to extent of occurence points

ev<-ENMevaluate(occ = locs, env = env, bg.coords = NULL, RMvalues = seq(1,3,1), fc = c("L","LQ","H"),
                method = "block", n.bg = 5000, rasterPreds = T, parallel = T, numCores = 8, algorithm = "maxent.jar")

#Function to find best model
# optimize
optimize <- function(res) {
  ###Remove any candidate model which has an AUC less than 0.51= models with no discrimination
  opt.auc <- res[res[,4] >= 0.5,]
  ###Remove any candidates which have no parameters
  no.param <- opt.auc[opt.auc[,13] > 1,]
  ###Remove any candidates where the AIC score was NA (too many parameters)
  noAICNA<- no.param[which(!is.na(no.param$AICc)),]
  ###Remove any models which have an OR of zero
  noOR0 <- noAICNA[noAICNA[,9] != 0,]
  ###Order the remaining list by lowest OR then highest AUC, sequentially
  ordered<-noOR0[with(noOR0, order(avg.test.or10pct, -avg.test.AUC)), ]
  ###Grab the settings of that first model (the optimal model)
  ordered[1,]
}
best<-optimize(ev@results)
best

plot(ev1@predictions)

plot(ev@predictions$L_3) #insert best model here, instead of L_3
points(ev@occ.pts, pch=21, bg=ev@occ.grp)