# install.packages(c("sf","terra", "stars", "lwgeom", "raster", "igraph"))
library(sf)
library(terra)
library(stars)
library(lwgeom)
library(raster)
library(gdistance)


###
# Linear interpolation function
#   rescale value c from range a-b to range y-z
###
lerp <- function(c, a, b, y, z) {
  (c - a) * (z - y) / (b - a) + y
}


###
# Logistic function used to link matrix patch area to distance decay of the edge effect
#   https://stackoverflow.com/questions/55725139/fit-sigmoid-function-s-shape-curve-to-data-using-python
# 
# Params:
#   `x` is the value for which you want to return the y
#   `L` is responsible for scaling the output range from [0,1] to [0,L]
#   `b` adds bias to the output and changes its range from [0,L] to [b,L+b]
#   `k` is responsible for scaling the input, which remains in (-inf,inf)
#   `x0` is the point in the middle of the Sigmoid, i.e. the point where Sigmoid should originally output the value 1/2 [since if x=x0, we get 1/(1+exp(0)) = 1/2].
# 
###
logistic <- function(x, L=1.0, b=0.0, k, x0) {
  L / (1.0 + exp(-k*(x-x0))) + b
}
colnePatches

plot(colnePatches$Area,colnePatches$VCI_AOP)

###
# This makes a call to the logistic function and then scales the result to enforce a value of 1 for the 
#   stated maximum value ([max]), 0 for 0 and 0.5 where x = x0 (maximum val/2). This is achieved using a 
#   linear  interpolation function to scale either from 0-0.5 (where x<0.5) or 0.5-[max] (otherwise). This 
#   is necessary as scaling from 0-max can displace the case so that y!=0.5 where x=x0 if the adjustment 
#   required at x=0 and x=[max] are not equal.
# 
# Params:
#   `x` is the value for which you want to return the y
#   `fullEdgeEffectArea` is the x value that should be coerced to 1
# 
###
coercedLogistic <- function(x, fullEdgeEffectArea) {
  
  # work out order of magnitude and convert to scale factor
  # TODO: this is an arbitrary parameter
  k <- 1 / (10^(floor(log(fullEdgeEffectArea, 10))-1))
  
  # x0 is the midpoint of the scale
  x0 <- fullEdgeEffectArea*0.5
  
  # get the 'raw' logistic value for the current x
  vx = logistic(x, k=k, x0=x0)
  
  # stretch it either up towards 1 or down towards 0 to coerce the scale to (x=0)=0, (x-x0)=0.5, (x=fullEdgeEffectArea)=1
  if (vx < 0.5) lerp(vx, logistic(0, k=k, x0=x0), 0.5, 0, 0.5) else lerp(vx, 0.5, logistic(fullEdgeEffectArea, k=k, x0=x0), 0.5, 1)
}


###
# Returns the connectivity Value for a given landscape and set of parameters (i.e., species)
# 
# Params: 
#   - `patches` = habitat patches (shapefile)
#   - `matrix` + matrix patches with edge and cost in attribute tables (shapefile)
#   - `specialism` = one of "interior", "edge" or "generalist"
#   - `dispersalDist` = mean dispersal distance of the species being modeled
#   - `habComponent` = component to determine the sensitivity of the species to patch size 
#       (i.e. strength of relationship between extinction and patch size) 
#   - `colonizationRate` = component setting colonization/survival success
#   - `edgeSensitivity` = sets how sensitive the species is to "specialism". 
#       Setting this to very small values for either edge or interior basically renders it a generalist
#   - `minPatchSize` = optional minimum viable patch size
#   - `perim` = if TRUE then connectivity matrix is uplifted by edge density (perimeter-edge ratio) of 
#       habitat. If FALSE then does nothing
#   - `fullEdgeEffectArea` = the are after which a full edge effect is realised
###


sl<-st_read("Data/slHeterog.shp")
plot(sl)
ss<-st_read("Data/ssHeterog.shp")

mixed<-st_read("Data/mixedHeterog.shp")

#only need columns 6-8 (edge, cost and geometry)
ss<-ss[,6:8]
sl<-sl[,6:8]
mixed<-mixed[,6:8]

#set habitat polygons to 1 (or skip this step to set to 0)
for(i in 1:nrow(ss)){

  if(ss$Cost[i]==0){
  ss$Cost[i]=1
}else{ss$Cost[i]=ss$Cost[i]}
}


for(i in 1:nrow(sl)){
  if(sl$Cost[i]==0){
    sl$Cost[i]=1
  }else{sl$Cost[i]=sl$Cost[i]}
}

for(i in 1:nrow(mixed)){
  if(mixed$Cost[i]==0){
    mixed$Cost[i]=1
  }else{mixed$Cost[i]=mixed$Cost[i]}
}


ggplot(sl) +
  geom_sf(aes(fill = Cost))+theme_pubr(base_size = 20,legend="right")+coord_sf(label_axes = "")


ggplot(ss) +
  geom_sf(aes(fill = Cost))+theme_pubr(base_size = 20,legend="right")+coord_sf(label_axes = "")


ggplot(mixed) +
  geom_sf(aes(fill = Cost))+theme_pubr(base_size = 20,legend="right")+coord_sf(label_axes = "")#+font("legend.text",size=29)



slPatches<-sl[sl$Cost==0,]
ssPatches<-ss[ss$Cost==0,]
mixedPatches<-mixed[mixed$Cost==0,]
patches<-ssPatches
unique(sl$Cost)

ssMatrix<-ss[ss$Cost>1,]
plot(ssMatrix)


slMatrix<-sl[sl$Cost>1,]
plot(slMatrix)

mixedMatrix<-mixed[mixed$Cost>1,]
(mixedMatrix)



ss<-ss[,6:8]
sl<-sl[,6:8]
mixed<-mixed[,6:8]
mixed
patches<-mixedPatches
matrix<-mixedMatrix
sl
ss
castorConnect<-function(patches, matrix, specialism, dispersalDist, colonizationRate, 
                        minPatchSize, habComponent, edgeSensitivity, perim, fullEdgeEffectArea) { 
  
  ## begin with a little quality control...
  
  # a landscape with no patches gets an 'easy' 0
  if (nrow(patches) == 0){
    return(0)
  }
  
  # catch the edge case (hah!) if the entire area is covered by a single patch - this
  #  causes a fail in the algorithm (patches doesn't return anything) so needs to be 
  #  caught beforehand
  if (nrow(patches) == 1){
    return(1)
  }
  
  
  ## rasterise the patches
  
  # first rasterize the matrix patches by edge attribute (column 1) and set NAs to 0 to keep everything same extent
  patchR <- st_rasterize(matrix[,1])
  patchR[is.na(patchR)] <- 0
  #res(rast(patchR))
  
  ### next set up for loop to collate all edge rasters coming OUT from matrix patches
  
  # initialise empty list
  edgeList <- list()
  
  # loop through the matrix of patches
  for(i in 1:nrow(matrix)) {
    
    # get patch area
    area.i <- st_area(matrix[i,])
    
    # Add area.i as argument to logistic function and get place on X axis by subtracting half range of values (10Ha as max)
    #  resulting logN is passed to exponential curve (distEdge) further down
    logN <- coercedLogistic(as.numeric(area.i), fullEdgeEffectArea=fullEdgeEffectArea)
    
    # rasterize the patch
    rMat <- rasterize(matrix[i,], field=matrix$Edge[i], rast(patchR))
    
    # get distance raster 
    rMatDist <- gridDistance(rMat, target=NA)
    
    # set zero (i.e. the patch) to 1 so that it has max edge
    rMatDist[rMatDist==0] <- 1
    
    # get edge for each patch from the attribute table (this needs to be provided beforehand)
    edge <- matrix$Edge[i]
    
    # isolate areas within cost threshold 
    rEdge <- rMatDist <= edge
    
    # get distances within edge threshold
    rEdge <- rEdge * rMatDist
    
    # set zeros to NA otherwise all areas outside the patch are given 1
    rEdge[rEdge==0] <- NA
    
    # this is the bit that sets shape of neg. exp. curve
    distEdge <- -log(logN)/ edge 
    
    # this creates the raster kernel
    edgeRast <- exp(-distEdge * rEdge)
    
    # add zeros back in so we can sum rasters later
    edgeRast[is.na(edgeRast)] <- 0
    #plot(edgeRast)
    # collect all rasters
    edgeList[[i]] <- edgeRast
    
  } # --close for loop
  
  # add up all edge rasters
  sumRast <- do.call("sum", edgeList)
  #plot(sumRast)
  # set max as one otherwise subtracting values from habitat patches later gives 
  #   negative area which doesn't make sense
  sumRast[sumRast > 1] <- 1 
  #plot(sumRast)
  
  
  # need to rasterize the patches to get at edge kernel and area
  r<-st_rasterize(patches) 
  r<-rast(r)  #convert to terra library.
  
  # increase resolution 5x
  # TODO: need to revisit this for using actual patch polygons
  r<-disagg(r, fact=5)
  
  # get patches (this will account for adjacent cells forming one patch, not in 
  #   hypothetical landscape but later on when modelling real ones)
  slClump<-patches(r, directions=8) 
  
  # polygonize the extent to get AOI later
  ext <- as.polygons(ext(vect(matrix)),crs=crs(patches)) 
  
  # set NA values to 0 to use when making the the distance raster later on
  slClump[is.na(slClump)] <- 0 
  
  #get "inverse edge" for edge species habitat
  invEdge <- 1-sumRast
  
  
  ## select by species type
  
  # interior specialist
  if(specialism == "interior") {
    
    # for an interior specialist we remove all edge cell values from patch cells 
    #   (i.e. remove a proportion of the area)
    edgeClump <- 1-(sumRast * edgeSensitivity)
    
    # 
    extClump <- extract(edgeClump, patches, fun=sum, method="simple", bind=T)
    
    # get cell area
    cellArea <- res(edgeClump)[1]^2 
    
    # create a vector for the new area
    extClump$areaMod <- extClump$layer * cellArea
    
    # edge specialist
  } else if(specialism=="edge") {
    
    # for an edge specialist just leave as the edge habitat within the patch (i.e. inverse of above)
    edgeClump <- 1-(invEdge * edgeSensitivity)
    
    # 
    extClump <- extract(edgeClump, patches, fun=sum, method="simple", bind=T)
    
    # get cell area
    cellArea <- res(edgeClump)[1]^2 
    
    # create a vector for the new area
    extClump$areaMod <- extClump$layer * cellArea
    
    # generalist
  } else if (specialism=="generalist") {
    
    # for generalist all patch cells are one i.e. area not changed so set extClump 
    #   to just be the habitat patches  
    extClump <- patches
    
    #get area of patches
    extClump$areaMod <- expanse(vect(patches))
  }
  
  # for loop to do any final modifications according to minimum patch size settings
  extClump$areaFin <- NA
  for(i in 1:nrow(extClump)) {
    #print(extClump$areaMod[i])
    
    # if the modified area is too small it should be penalised
    if(extClump$areaMod[i] < minPatchSize) {
      
      # apply penalisation for patches below minPatchSize
      extClump$areaFin[i] <- extClump$areaMod[i] - (minPatchSize-extClump$areaMod[i])^habComponent}
    
    else {
      
      # just accept current value
      extClump$areaFin[i] <- extClump$areaMod[i]
      
    }
    
    # enforce 0 as minimum area
    if(extClump$areaFin[i] < 0) { extClump$areaFin[i] <- 0 }
  }
  
  ###########################################Least Cost Path calculations##############################
  #head(ssMatrix)
  # convert matrix patches to cost raster (costs in col 2)
  #head(matrix)
  costRast<-st_rasterize(matrix[,2])
  
  # gdistance takes raster package objects so convert 
  costRast<-rast(costRast)
  costRast[is.na(costRast)]<-1
  #plot(costRast)
  
  cost1<-costRast
  cost1[costRast>1]=10
  #plot(cost1)
  plot(sumRast)
  costPatch<-cost1*sumRast
  #plot(costPatch)
  costPatch[costPatch>1]<-1
  finCost<-costPatch*costRast
  plot(finCost)
  #hist(finCost)
  # set patches (these are NA in the matrix layer) to zero cost so that patches are 
  #   "free movement" (also same as edge-to-edge distance)
  #costRast[costRast==0] <- sumRast
  
  #costRastIntra<-costRast*sumRast
  #plot(costRastIntra)
  
  #costRast1<-(costRastIntra>1)*costRast
  #plot(costRast1)
  #hist(costRastIntra)
  # calculate transition matrix using inline function
  land_cost <- transition(raster(costRast), transitionFunction=function(x) 1 / mean(x), 8)
  
  # set destination points as centroids of the patches
  # TODO: centroid might not always be inside the polygon - we need a representative point instead
  sites <- SpatialPoints(st_coordinates(st_point_on_surface(patches)))
  
  # init cost matrix and loop through each row
  costMat <- matrix(0, nrow=nrow(patches), ncol=nrow(patches))
  n.patch<-nrow(costMat)  
  for (i in 1:n.patch) {
    c <- gdistance::costDistance(land_cost, sites[i], sites) 
    costMat[i,] <- c
    # print(i / n.patch * 100)
  }
  
  # this is a matrix of least cost distances between patches 
  distMat <- costMat
  
  #create basis for probability matrix
  distMat<-as.matrix(distMat, nrow=nrow(patches), nrow=nrow(patches)) 
  
  # in case modelling on single patch then need if else to side step matrix error here
  #if(nrow(distMat) > 1) { 
  distMat <- apply(distMat, MARGIN=1, FUN=as.numeric) # annoyingly need to make sure all elements are numeric here
  
  #} else {
  #distMat <- distMat
  #} 
  
  # set alpha which determines colonization probability of the species 
  alpha= -log(colonizationRate) / dispersalDist
  
  # 
  #if(nrow(distMat) == 1){
  
  # for single patches just set diagonal to 1
  # A.prob <- matrix(1, nrow=nrow(distMat), ncol=ncol(distMat))
  
  # TODO:is this not already a matrix?
  #A.prob <- as.matrix(A.prob)
  
  # store probabilities in adjacency table
  #graph.Aprob <- graph.adjacency(A.prob, mode="undirected", weighted=T)
  
  #} else {
  
  # 
  A.prob <- matrix(0, nrow=nrow(distMat), ncol=ncol(distMat))
  
  # negative exponential of colonization kernel x distance get probability
  A.prob <- exp(-alpha * distMat) 
  
  # set diag to zero
  # TODO: but this is another parameter that could be tweaked?
  diag(A.prob) <- 0 
  
  # final matrix for connectivity graph
  A.prob <- as.matrix(A.prob)
  
  # final matrix for connectivity graph
  graph.Aprob <- graph.adjacency(A.prob, mode="undirected", weighted=T)
  #}
  
  
  ## calculate landscape-scale connectivity
  
  # calculate all shortest paths between nodes
  pstar.mat <- shortest.paths(graph.Aprob, weights= -log(E(graph.Aprob)$weight))
  
  # back-transform to probabilities of connectedness
  pstar.mat <- exp(-pstar.mat)                                                   
  
  # sum probabilities of connectedness
  pMatSum <- sum(pstar.mat)
  
  # get study area in m2
  AL <- expanse(ext)
  
  # get area vector for PC metric
  area <- extClump$areaFin
  
  # sum areas from the vector
  areaSum <- sum(area)
  
  # get product of all patch areas ij and multiply by probabilities above
  PCmatNew <- outer(area,area) * pstar.mat 
  
  # calculate sum of edge ratio (perimeter / area) for each patch
  edgeRatio <- sum(st_perimeter(patches) / (st_area(patches)))
  
  # uplift by 1+ edge ratio if perim is set to TRUE
  # JJH: I added brackets here to correct the calculation
  pcMatSum <- sum(PCmatNew)
  
  # divide by total area of the study squared to get the PC metric  
  pcMod <- pcMatSum / as.numeric(AL^2) 
  #pcMod <- sqrt(pcMatSum) / as.numeric(AL)
  # compile results into list and return
  return(list(areaSum, edgeRatio, pMatSum, pcMatSum, pcMod))
  #print(x)
}


######################## Do some testing  

# load data
#SL<-st_read("data/SL_SmallGridDissolve.shp")
#slMat<-st_read("data/SL_MatrixSmallDissolve.shp")

head(slPatches)
slPatches<-slPatches[,c(6,7,8)]

head(slMatrix)
slMatrix<-slMatrix[,c(6,7,8)]

# run function
testConnectIntSL<-castorConnect(
  patches=slPatches,
  matrix=slMatrix,
  specialism="interior",
  dispersalDist=5000,
  minPatchSize=0, 
  habComponent=0.9,
  colonizationRate=0.1,
  edgeSensitivity=0.09,
  perim=FALSE,
  fullEdgeEffectArea=1000000
)

######################## Do some testing  

# load data
ssPatches

#plot(mixedMatrix)

# run function
testConnect<-castorConnect(
  patches=mixedPatches,
  matrix=mixedMatrix,
  specialism="edge",
  dispersalDist=5000,
  minPatchSize=0, 
  habComponent=0.9,
  colonizationRate=0.05,
  edgeSensitivity=0.5,
  perim=FALSE,
  fullEdgeEffectArea=10000
)

testConnect

plot(ssMatrix)

edgeEffectFun<-function(x){
  testConnect<-castorConnect(
    patches=ssPatches,
    matrix=ssMatrix,
    specialism="interior",
    dispersalDist=1000,
    minPatchSize=00, 
    habComponent=0.9,
    colonizationRate=x,
    edgeSensitivity=0.99,
    perim=FALSE,
    fullEdgeEffectArea=1000000,
    return(testConnect)
  )
  
}

sampN<-seq(0.05,1,by=0.05)
sampN

###################################################
###############################################
############################################LO SENS

edgeEffectFun<-function(x){
  testConnect<-castorConnect(
    patches=mixedPatches,
    matrix=mixedMatrix,
    specialism="edge",
    dispersalDist=500,
    minPatchSize=00, 
    habComponent=0.9,
    colonizationRate=x,
    edgeSensitivity=0.1,
    perim=FALSE,
    fullEdgeEffectArea=10000000
    
  )
  print(x)
  return(testConnect)
}



edgeMixedLoSenseDisp500<-lapply(sampN,edgeEffectFun)

edgeMixedLoSenseDisp500

edgeMixedLoSenseDisp500<-do.call("rbind",edgeMixedLoSenseDisp500) # have this

edgeMixedLoSenseDisp500DF<-data.frame(edgeMixedLoSenseDisp500)






###########################################################

edgeEffectFun<-function(x){
  testConnect<-castorConnect(
    patches=mixedPatches,
    matrix=mixedMatrix,
    specialism="interior",
    dispersalDist=500,
    minPatchSize=00, 
    habComponent=0.9,
    colonizationRate=x,
    edgeSensitivity=0.1,
    perim=FALSE,
    fullEdgeEffectArea=10000000
    
  )
  print(x)
  return(testConnect)
}

intMixedLoSenseDisp500<-lapply(sampN,edgeEffectFun)

intMixedLoSenseDisp500

intMixedLoSenseDisp500<-do.call("rbind",intMixedLoSenseDisp500) # have this

intMixedLoSenseDisp500DF<-data.frame(intMixedLoSenseDisp500)

intMixedLoSenseDisp500DF




###################################################

edgeEffectFun<-function(x){
  testConnect<-castorConnect(
    patches=ssPatches,
    matrix=ssMatrix,
    specialism="interior",
    dispersalDist=500,
    minPatchSize=00, 
    habComponent=0.9,
    colonizationRate=x,
    edgeSensitivity=0.1,
    perim=FALSE,
    fullEdgeEffectArea=10000000
    
  )
  print(x)
  return(testConnect)
}


intSSLoSenseDisp500<-lapply(sampN,edgeEffectFun)

intSSLoSenseDisp500

intSSLoSenseDisp500<-do.call("rbind",intSSLoSenseDisp500) # have this

intSSLoSenseDisp500DF<-data.frame(intSSLoSenseDisp500)

intSSLoSenseDisp500DF

#plot(sampN,intSSLoSenseDisp500DF[,5])
#plot(sampN,intSLLoSenseDisp500DF[,5])





############

edgeEffectFun<-function(x){
  testConnect<-castorConnect(
    patches=ssPatches,
    matrix=ssMatrix,
    specialism="edge",
    dispersalDist=500,
    minPatchSize=00, 
    habComponent=0.9,
    colonizationRate=x,
    edgeSensitivity=0.1,
    perim=FALSE,
    fullEdgeEffectArea=10000000
    
  )
  print(x)
  return(testConnect)
}


edgeSSLoSenseDisp500<-lapply(sampN,edgeEffectFun)

edgeSSLoSenseDisp500

edgeSSLoSenseDisp500<-do.call("rbind",edgeSSLoSenseDisp500) # have this

edgeSSLoSenseDisp500DF<-data.frame(edgeSSLoSenseDisp500)

edgeSSLoSenseDisp500DF




###############################################################

edgeEffectFun<-function(x){
  testConnect<-castorConnect(
    patches=slPatches,
    matrix=slMatrix,
    specialism="edge",
    dispersalDist=500,
    minPatchSize=00, 
    habComponent=0.9,
    colonizationRate=x,
    edgeSensitivity=0.1,
    perim=FALSE,
    fullEdgeEffectArea=10000000
    
  )
  print(x)
  return(testConnect)
}

edgeSLLoSenseDisp500<-lapply(sampN,edgeEffectFun)

edgeSLLoSenseDisp500

edgeSLLoSenseDisp500<-do.call("rbind",edgeSLLoSenseDisp500) # have this

edgeSLLoSenseDisp500DF<-data.frame(edgeSLLoSenseDisp500)





#################################################

edgeEffectFun<-function(x){
  testConnect<-castorConnect(
    patches=slPatches,
    matrix=slMatrix,
    specialism="interior",
    dispersalDist=500,
    minPatchSize=00, 
    habComponent=0.9,
    colonizationRate=x,
    edgeSensitivity=0.1,
    perim=FALSE,
    fullEdgeEffectArea=10000000
    
  )
  print(x)
  return(testConnect)
}

intSLLoSenseDisp500<-lapply(sampN,edgeEffectFun)


intSLLoSenseDisp500<-do.call("rbind",intSLLoSenseDisp500) # have this

intSLLoSenseDisp500

intSLLoSenseDisp500DF<-data.frame(intSLLoSenseDisp500)

intSLLoSenseDisp500DF

###########################################
##############################################
#################HI SENS######################

edgeEffectFun<-function(x){
  testConnect<-castorConnect(
    patches=mixedPatches,
    matrix=mixedMatrix,
    specialism="edge",
    dispersalDist=500,
    minPatchSize=00, 
    habComponent=0.9,
    colonizationRate=x,
    edgeSensitivity=0.9,
    perim=FALSE,
    fullEdgeEffectArea=10000000
    
  )
  print(x)
  return(testConnect)
}


edgeMixedHiSenseDisp500<-lapply(sampN,edgeEffectFun)

edgeMixedHiSenseDisp500

edgeMixedHiSenseDisp500<-do.call("rbind",edgeMixedHiSenseDisp500) # have this

edgeMixedHiSenseDisp500DF<-data.frame(edgeMixedHiSenseDisp500)

head(edgeMixedHiSenseDisp500DF)



###########################################################

edgeEffectFun<-function(x){
  testConnect<-castorConnect(
    patches=mixedPatches,
    matrix=mixedMatrix,
    specialism="interior",
    dispersalDist=500,
    minPatchSize=00, 
    habComponent=0.9,
    colonizationRate=x,
    edgeSensitivity=0.9,
    perim=FALSE,
    fullEdgeEffectArea=10000000
    
  )
  print(x)
  return(testConnect)
}


intMixedHiSenseDisp500<-lapply(sampN,edgeEffectFun)

intMixedHiSenseDisp500

intMixedHiSenseDisp500<-do.call("rbind",intMixedHiSenseDisp500) # have this

intMixedHiSenseDisp500DF<-data.frame(intMixedHiSenseDisp500)







###################################################
#################################################### Do NEXT!

edgeEffectFun<-function(x){
  testConnect<-castorConnect(
    patches=ssPatches,
    matrix=ssMatrix,
    specialism="interior",
    dispersalDist=500,
    minPatchSize=00, 
    habComponent=0.9,
    colonizationRate=x,
    edgeSensitivity=0.9,
    perim=FALSE,
    fullEdgeEffectArea=10000000
    
  )
  print(x)
  return(testConnect)
}


intSSHiSenseDisp500<-lapply(sampN,edgeEffectFun)

intSSHiSenseDisp500

intSSHiSenseDisp500<-do.call("rbind",intSSHiSenseDisp500) # have this

intSSHiSenseDisp500DF<-data.frame(intSSHiSenseDisp500)

intSSHiSenseDisp500DF



##################################################################

edgeEffectFun<-function(x){
  testConnect<-castorConnect(
    patches=ssPatches,
    matrix=ssMatrix,
    specialism="edge",
    dispersalDist=500,
    minPatchSize=00, 
    habComponent=0.9,
    colonizationRate=x,
    edgeSensitivity=0.9,
    perim=FALSE,
    fullEdgeEffectArea=10000000
    
  )
  print(x)
  return(testConnect)
}

edgeSSHiSenseDisp500<-lapply(sampN,edgeEffectFun)

edgeSSHiSenseDisp500

edgeSSHiSenseDisp500<-do.call("rbind",edgeSSHiSenseDisp500) # have this

edgeSSHiSenseDisp500DF<-data.frame(edgeSSHiSenseDisp500)




###############################################################

edgeEffectFun<-function(x){
  testConnect<-castorConnect(
    patches=slPatches,
    matrix=slMatrix,
    specialism="interior",
    dispersalDist=500,
    minPatchSize=00, 
    habComponent=0.9,
    colonizationRate=x,
    edgeSensitivity=0.9,
    perim=FALSE,
    fullEdgeEffectArea=10000000
    
  )
  print(x)
  return(testConnect)
}

intSLHiSenseDisp500<-lapply(sampN,edgeEffectFun)

intSLHiSenseDisp500

intSLHiSenseDisp500<-do.call("rbind",intSLHiSenseDisp500) # have this

intSLHiSenseDisp500DF<-data.frame(intSLHiSenseDisp500)
intSLHiSenseDisp500DF

########################################################

edgeEffectFun<-function(x){
  testConnect<-castorConnect(
    patches=slPatches,
    matrix=slMatrix,
    specialism="edge",
    dispersalDist=500,
    minPatchSize=00, 
    habComponent=0.9,
    colonizationRate=x,
    edgeSensitivity=0.9,
    perim=FALSE,
    fullEdgeEffectArea=10000000
    
  )
  print(x)
  return(testConnect)
}

edgeSLHiSenseDisp500<-lapply(sampN,edgeEffectFun)

edgeSLHiSenseDisp500

edgeSLHiSenseDisp500<-do.call("rbind",edgeSLHiSenseDisp500) # have this

edgeSLHiSenseDisp500DF<-data.frame(edgeSLHiSenseDisp500)

edgeSLHiSenseDisp500DF


#################################################Generalists


edgeEffectFun<-function(x){
  testConnect<-castorConnect(
    patches=slPatches,
    matrix=slMatrix,
    specialism="generalist",
    dispersalDist=500,
    minPatchSize=00, 
    habComponent=0.9,
    colonizationRate=x,
    edgeSensitivity=0.95,
    perim=FALSE,
    fullEdgeEffectArea=10000000
    
  )
  print(x)
  return(testConnect)
}

genSLDisp500<-lapply(sampN,edgeEffectFun)

genSLDisp500

genSLDisp500<-do.call("rbind",genSLDisp500) # have this

genSLDisp500DF<-data.frame(genSLDisp500)

################################################



edgeEffectFun<-function(x){
  testConnect<-castorConnect(
    patches=ssPatches,
    matrix=ssMatrix,
    specialism="generalist",
    dispersalDist=500,
    minPatchSize=00, 
    habComponent=0.9,
    colonizationRate=x,
    edgeSensitivity=0.95,
    perim=FALSE,
    fullEdgeEffectArea=10000000
    
  )
  print(x)
  return(testConnect)
}

genSSDisp500<-lapply(sampN,edgeEffectFun)

genSSDisp500

genSSDisp500<-do.call("rbind",genSSDisp500) # have this

genSSDisp500DF<-data.frame(genSSDisp500)

#################################################################



edgeEffectFun<-function(x){
  testConnect<-castorConnect(
    patches=mixedPatches,
    matrix=mixedMatrix,
    specialism="generalist",
    dispersalDist=500,
    minPatchSize=00, 
    habComponent=0.9,
    colonizationRate=x,
    edgeSensitivity=0.95,
    perim=FALSE,
    fullEdgeEffectArea=10000000
    
  )
  print(x)
  return(testConnect)
}

genMixedDisp500<-lapply(sampN,edgeEffectFun)

genMixedDisp500

genMixedDisp500<-do.call("rbind",genMixedDisp500) # have this

genMixedDisp500DF<-data.frame(genMixedDisp500)




###########################################################

figData<-cbind(sampN,edgeMixedLoSenseDisp500DF[,5],intMixedLoSenseDisp500DF[,5],
               intSSLoSenseDisp500DF[,5],edgeSSLoSenseDisp500DF[,5],
               edgeSLLoSenseDisp500DF[,5],edgeMixedHiSenseDisp500DF[,5],
               intMixedHiSenseDisp500DF[,5],intSSHiSenseDisp500DF[,5],edgeSSHiSenseDisp500DF[,5],
               edgeSLHiSenseDisp500DF[,5],intSLHiSenseDisp500DF[,5],intSLLoSenseDisp500DF[,5])
               #genSLDisp500DF[,5],genSSDisp500DF[,5],genMixedDisp500DF[,5])

intData<-cbind(sampN,intMixedLoSenseDisp500DF[,5],
               intSSLoSenseDisp500DF[,5],
               intSLLoSenseDisp500DF[,5])


intData<-as.data.frame(intData)

intData<-apply(intData,MARGIN=2,FUN=unlist)

intData<-as.data.frame(intData)

colnames(intData)<-c("dispersal_Rate","intMixedLo","intSSLo","intSLLo")

loIntPlot<- ggplot(intData)+geom_line(aes(dispersal_Rate,intSSLo,colour="SS"),linewidth=2)+
  
  geom_line(aes(dispersal_Rate,intMixedLo,colour="Mixed"),linewidth=2)+
  
  geom_line(aes(dispersal_Rate,intSLLo,colour="SL"),linewidth=2)+
  
  xlab("Dispersal Probability at maxDist")+ylab("RHI")+scale_color_manual(name = "Landscape", 
                                                                          
                                                                          values = c("SL" = "darkblue", "SS" = "red","Mixed"="black"))#+theme(text = element_text(size = 20))   



loIntPlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)

###########################################################


edgeData<-cbind(sampN,edgeMixedLoSenseDisp500DF[,5],
               edgeSSLoSenseDisp500DF[,5],
               edgeSLLoSenseDisp500DF[,5])


edgeData<-as.data.frame(edgeData)

edgeData<-apply(edgeData,MARGIN=2,FUN=unlist)

edgeData<-as.data.frame(edgeData)

colnames(edgeData)<-c("dispersal_Rate","edgeMixedLo","edgeSSLo","edgeSLLo")

loEdgePlot<- ggplot(edgeData)+geom_line(aes(dispersal_Rate,edgeSSLo,colour="SS"),linewidth=2)+
  
  geom_line(aes(dispersal_Rate,edgeMixedLo,colour="Mixed"),linewidth=2)+
  
  geom_line(aes(dispersal_Rate,edgeSLLo,colour="SL"),linewidth=2)+
  
  xlab("Dispersal Probability at maxDist")+ylab("RHI")+scale_color_manual(name = "Landscape", 
                                                                          
                                                                          values = c("SL" = "darkblue", "SS" = "red","Mixed"="black"))#+theme(text = element_text(size = 20))   



loEdgePlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)






####################

figData<-as.data.frame(figData)

figData<-apply(figData,MARGIN=2,FUN=unlist)

figData<-as.data.frame(figData)

colnames(figData)<-c("dispersal_Rate","edgeMixedLo","intMixedLo","intSSLo","edgeSSLo",
                     "edgeSLLo","edgeMixedHi","intMixedHi","intSSHi","edgeSSHi","edgeSLHi",
                     "intSLHi","intSLLo")#"genSL","genSS","genMixed")

(figData)

library(ggplot2)
library(ggpubr)

loIntPlot<- ggplot(figData)+geom_line(aes(dispersal_Rate,intSSLo,colour="SS"),linewidth=2)+
  
  geom_line(aes(dispersal_Rate,intMixedLo,colour="Mixed"),linewidth=2)+
  
  geom_line(aes(dispersal_Rate,intSLLo,colour="SL"),linewidth=2)+
  
  xlab("Dispersal Probability at maxDist")+ylab("RHI")+scale_color_manual(name = "Landscape", 
                                                        
                                                        values = c("SL" = "darkblue", "SS" = "red","Mixed"="black"))#+theme(text = element_text(size = 20))   



loIntPlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)



######################################################################################################

hiIntPlot<- ggplot(figData)+geom_line(aes(dispersal_Rate,intSSHi,colour="SS"),linewidth=2)+
  
  geom_line(aes(dispersal_Rate,intMixedHi,colour="Mixed"),linewidth=2)+
  
  geom_line(aes(dispersal_Rate,intSLHi,colour="SL"),linewidth=2)+
  
  xlab("Dispersal Probability at maxDist")+ylab("RHI")+scale_color_manual(name = "Landscape", 
                                                                          
                                                                          values = c("SL" = "darkblue", "SS" = "red","Mixed"="black"))#+theme(text = element_text(size = 20))   



hiIntPlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)

######################################################################################################

hiEdgePlot<- ggplot(figData)+geom_line(aes(dispersal_Rate,edgeSSHi,colour="SS"),linewidth=2)+
  
  geom_line(aes(dispersal_Rate,edgeMixedHi,colour="Mixed"),linewidth=2)+
  
  geom_line(aes(dispersal_Rate,edgeSLHi,colour="SL"),linewidth=2)+
  
  xlab("Dispersal Probability at maxDist")+ylab("RHI")+scale_color_manual(name = "Landscape", 
                                                                          
                                                                          values = c("SL" = "darkblue", "SS" = "red","Mixed"="black"))#+theme(text = element_text(size = 20))   



hiEdgePlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)


############################################################################################################

loEdgePlot<- ggplot(figData)+geom_line(aes(dispersal_Rate,edgeSSLo,colour="SS"),linewidth=2)+
  
  geom_line(aes(dispersal_Rate,edgeMixedLo,colour="Mixed"),linewidth=2)+
  
  geom_line(aes(dispersal_Rate,edgeSLLo,colour="SL"),linewidth=2)+
  
  xlab("Dispersal Probability at maxDist")+ylab("RHI")+scale_color_manual(name = "Landscape", 
                                                                          
                                                                          values = c("SL" = "darkblue", "SS" = "red","Mixed"="black"))#+theme(text = element_text(size = 20))   



loEdgePlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)


#############################################################

genPlot<- ggplot(figData)+geom_line(aes(dispersal_Rate,genSS,colour="SS"),linewidth=2)+
  
  geom_line(aes(dispersal_Rate,genMixed,colour="Mixed"),linewidth=2)+
  
  geom_line(aes(dispersal_Rate,genSL,colour="SL"),linewidth=2)+
  
  xlab("Dispersal Probability at maxDist")+ylab("RHI")+scale_color_manual(name = "Landscape", 
                                                                          
                                                                          values = c("SL" = "darkblue", "SS" = "red","Mixed"="black"))#+theme(text = element_text(size = 20))   



genPlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)





