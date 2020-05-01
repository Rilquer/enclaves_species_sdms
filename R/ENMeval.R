library(ENMeval)
setwd('/home/pgalante/Projects/PracticeData')
locs<-read.csv('locs.csv')
colnames(locs) <- c("x","y")
env<-stack(list.files(path = "/home/pgalante/layers/wc2-5/", pattern = '\\.tif$', full.names = T))
env<-crop(env, extent(locs))

env <- crop(stack(list.files(path = 'data/wc0-5/', pattern = '\\.tif$', full.names = T)),enclaves_extent)
ev<-ENMevaluate(occ = locs[[1]], env = env, bg.coords = NULL, RMvalues = seq(1,3,1), fc = c("L","LQ","H"),
                method = "block", n.bg = 5000, rasterPreds = T, parallel = T, numCores = 8, algorithm = "maxnet")

ev@results

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
  ordered<-noOR0[with(noOR0, order(Mean.OR10, -Mean.AUC)), ]
  ###Grab the settings of that first model (the optimal model)
  ordered[1,]
}
best<-optimize(ev@results)

#####Modified from Adam B Smith @ Missouri Botanical Garden (earthskysea.org)
###Create that optimal model using all localities (no withheld data)
###x = rasterstack, p = locs object as data.frame, a = background coords, factors = categoricals,
### make true only the arguments wanted
mod <- maxent(
  x=env, # bio stack
  p=locs, # locality csv
  a= NULL, # background coords
  #path= "C:/Users/pgalante/Documents/Projects/Madagascar/endemism_diversity/tests", # path to save to
  args=c(
    'betamultiplier=3',
    'linear=true',
    'quadratic=false',
    'product=false',
    'threshold=false',
    'hinge=false',
    'threads=2',
    'responsecurves=true',
    'jackknife=true',
    'askoverwrite=false'
  )
)


##Predict that optimal model: select output format wanted, output file type wantes, and pathway to save it
Pred.Mod<- predict(
  object = mod,
  x = env,
  #filename = "C:/Users/pgalante/Documents/Projects/Madagascar/endemism_diversity/tests", #where do you want the prediction saved?
  na.rm = T,
  format = 'GTiff',#or GTiff
  overwrite = F,
  args = "cloglog"
)

