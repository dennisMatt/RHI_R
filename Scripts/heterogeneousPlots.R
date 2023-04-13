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

sl<-st_read("Data/slHeterog.shp")

ss<-st_read("Data/ssHeterog.shp")

mixed<-st_read("Data/mixedHeterog.shp")

#only need columns 6-8 (edge, cost and geometry)
ss<-ss[,6:8]
sl<-sl[,6:8]
mixed<-mixed[,6:8]

#set within-habitat patch cost value to 1 (this is just to reproduce Fig. 4 in Dennis et al., 2023, we otherwise deal with habitat movement 
#cost values inside the RHI connectivity functions)

# for SS landscapes
for(i in 1:nrow(ss)){

  if(ss$Cost[i]==0){
  ss$Cost[i]=1
}else{ss$Cost[i]=ss$Cost[i]}
}

#for SL landscapes
for(i in 1:nrow(sl)){
  if(sl$Cost[i]==0){
    sl$Cost[i]=1
  }else{sl$Cost[i]=sl$Cost[i]}
}

#for Mixed landscapes
for(i in 1:nrow(mixed)){
  if(mixed$Cost[i]==0){
    mixed$Cost[i]=1
  }else{mixed$Cost[i]=mixed$Cost[i]}
}


#########create plots for heterogeneous landscapes
ggplot(sl) +
  geom_sf(aes(fill = Cost))+theme_pubr(base_size = 20,legend="right")+coord_sf(label_axes = "")


ggplot(ss) +
  geom_sf(aes(fill = Cost))+theme_pubr(base_size = 20,legend="right")+coord_sf(label_axes = "")


ggplot(mixed) +
  geom_sf(aes(fill = Cost))+theme_pubr(base_size = 20,legend="right")+coord_sf(label_axes = "")

#isolate patches into new object (if setting habitat cost to 0 then replace 1 here with 0)

slPatches<-sl[sl$Cost==1,]
ssPatches<-ss[ss$Cost==1,]
mixedPatches<-mixed[mixed$Cost==1,]

#create objects for matrix
ssMatrix<-ss[ss$Cost>1,]

slMatrix<-sl[sl$Cost>1,]


mixedMatrix<-mixed[mixed$Cost>1,]

#set up the RHI function to reproduce analyses in Dennis et al. 2023 Fig 4

#######
#######This function has nine arguments - patches = habitat patches: shapefile 
######################################## matrix = edge source polygons 
#####################################specialism = one of "interior", "edge" or "generalist"
################################### edge = the size of the edge effect (for a neutral landscape, will have to update this to handle "edge rasters" in real landscapes)
#################################### edgeIntensity = a component (suggest between 0.01 -1) that determines the shape of the kernel (distance decay of edge effect)
#################################### maxDist = maximum dispersal distance of the species being modeled
#####################################dispersalRate = component setting rate of dispersal success
#####################################edgeSensitivity = sets how sensitive the species is to "specialism". Setting this to increasingly small values for either edge or interior moves the species further towards generalist
###################################### fullEdgeEffectArea = sets the area above which edge source patches exert maximal (full) edge effect
rhiConnectHet<-function(patches,matrix,specialism,edge,edgeIntensity,maxDist,
                     dispersalRate,edgeSensitivity,fullEdgeEffectArea){ 
  
  ## begin with a little quality control...
  
  # a landscape with no patches gets an 'easy' 0
  if (nrow(patches) == 0){
    return(0)
  }
  
  #  if the entire area is covered by a single patch - this
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
    
    # extract raster values from within patches
    extClump <- extract(edgeClump, patches, fun=sum, method="simple", bind=T)
    
    # get cell area
    cellArea <- res(edgeClump)[1]^2 
    
    # create a vector for the new area
    extClump$areaMod <- extClump$layer * cellArea
    
    # edge specialist
  } else if(specialism=="edge") {
    
    # for an edge specialist just leave as the edge habitat within the patch (i.e. inverse of above)
    edgeClump <- 1-(invEdge * edgeSensitivity)
    
    # extract raster values from within patches
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
  

  
  
  ###########################################Least Cost Path calculations##############################
  
  # convert matrix patches to cost raster (costs in col 2)
  
  costRast<-st_rasterize(matrix[,2])
  
   
  costRast<-rast(costRast)
  
  #set habitat patches (NA in the data to 1)
  costRast[is.na(costRast)]<-1
  
  
  # gdistance takes raster package objects so convert 
  land_cost <- transition(raster(costRast), transitionFunction=function(x) 1 / mean(x), 8)
  
  # set destination points as centroids of the patches
  
  sites <- SpatialPoints(st_coordinates(st_point_on_surface(patches)))
  
  # init cost matrix and loop through each row
  costMat <- matrix(0, nrow=nrow(patches), ncol=nrow(patches))
  n.patch<-nrow(costMat)  
  for (i in 1:n.patch) {
    c <- gdistance::costDistance(land_cost, sites[i], sites) 
    costMat[i,] <- c
    
  }
  
  # this is a matrix of least cost distances between patches 
  distMat <- costMat
  
  #create basis for probability matrix
  distMat<-as.matrix(distMat, nrow=nrow(patches), nrow=nrow(patches)) 
  
  
  
  distMat <- apply(distMat, MARGIN=1, FUN=as.numeric) # make sure all elements are numeric here
  
  # set alpha which determines colonization probability of the species 
  alpha= -log(dispersalRate) / maxDist
  
  # init empty matrix for adjacency matrix 
  A.prob <- matrix(0, nrow=nrow(distMat), ncol=ncol(distMat))
  
  # negative exponential of colonization kernel x distance get probability
  A.prob <- exp(-alpha * distMat) 
  
  # set diag to zero
  diag(A.prob) <- 0 
  
  # final matrix for connectivity graph
  A.prob <- as.matrix(A.prob)
  
  # final matrix for connectivity graph
  graph.Aprob <- graph.adjacency(A.prob, mode="undirected", weighted=T)
  #}
  
  
  ## calculate RHI
  
  # calculate all shortest paths between nodes
  pstar.mat <- shortest.paths(graph.Aprob, weights= -log(E(graph.Aprob)$weight))
  
  # back-transform to probabilities of connectedness
  pstar.mat <- exp(-pstar.mat)                                                   
  
 
  
  # get study area in m2
  AL <- expanse(ext)
  
  # get area vector 
  area <- extClump$areaMod
  
  # sum areas from the vector
  areaSum <- sum(area)
  
  # get product of all patch areas ij and multiply by probabilities above
  PCmat <- outer(area,area) * pstar.mat 
  

  #sum above
  pcMatSum <- sum(PCmat)
  
  # divide by total area of the study squared to get the PC metric  
  pcMod <- pcMatSum / as.numeric(AL^2) 
  
  
  # compile results into list and return
  return(list(areaSum, pcMatSum, pcMod))
  #print(x)
}




######################## Do some testing  

#plot(mixedMatrix)

# run function
testConnect<-rhiConnect(
  patches=mixedPatches,
  matrix=mixedMatrix,
  specialism="edge",
  maxDist=5000,
  dispersalRate=0.05,
  edgeSensitivity=0.5,
  fullEdgeEffectArea=10000
)

#####create function that takes disperalRate values and passes to RHI function
###########################Dispersal Fun

dispFun<-function(x,configPatches,configMatrix,speciesGroup,sens,dist,FEEA){
  
  intRHI<-rhiConnectHet(patches = configPatches,matrix = configMatrix,specialism = speciesGroup,maxDist = dist,dispersalRate = x,edgeSensitivity = sens,fullEdgeEffectArea = FEEA)
  print(x)
  
  return(intRHI)
  
}


############################################### Run a sequence of models iterating over values for dispersalRate
###############################################

#create vector of dispersal rate values from 0.01-0.99

dispN<-c(0.01,seq(0.05,0.95,by=0.05),0.99)

############################################LOW SENSITIVITY TESTS 

##edge/mixed
edgeMixedLoSenseDisp500<-lapply(dispN,dispFun,configPatches=mixedPatches,configMatrix=mixedMatrix,
                                speciesGroup="edge",sens=0.1,dist=500,FEEA=100000)
                                

edgeMixedLoSenseDisp500

edgeMixedLoSenseDisp500<-do.call("rbind",edgeMixedLoSenseDisp500) # have this

edgeMixedLoSenseDisp500DF<-data.frame(edgeMixedLoSenseDisp500)





###########################################################

#int/mixed
intMixedLoSenseDisp500<-lapply(dispN,dispFun,configPatches=mixedPatches,configMatrix=mixedMatrix,
                                speciesGroup="interior",sens=0.1,dist=500,FEEA=100000)


intMixedLoSenseDisp500

intMixedLoSenseDisp500<-do.call("rbind",intMixedLoSenseDisp500) # have this

intMixedLoSenseDisp500DF<-data.frame(intMixedLoSenseDisp500)

intMixedLoSenseDisp500DF




###################################################


#int/SS
intSSLoSenseDisp500<-lapply(dispN,dispFun,configPatches=ssPatches,configMatrix=ssMatrix,
                              speciesGroup="interior",sens=0.1,dist=500,FEEA=100000)

intSSLoSenseDisp500

intSSLoSenseDisp500<-do.call("rbind",intSSLoSenseDisp500) # have this

intSSLoSenseDisp500DF<-data.frame(intSSLoSenseDisp500)

intSSLoSenseDisp500DF






############

#edge/SS

edgeSSLoSenseDisp500<-lapply(dispN,dispFun,configPatches=ssPatches,configMatrix=ssMatrix,
                             speciesGroup="edge",sens=0.1,dist=500,FEEA=100000)

edgeSSLoSenseDisp500

edgeSSLoSenseDisp500<-do.call("rbind",edgeSSLoSenseDisp500) # have this

edgeSSLoSenseDisp500DF<-data.frame(edgeSSLoSenseDisp500)

edgeSSLoSenseDisp500DF




###############################################################

#edge/SL

edgeSLLoSenseDisp500<-lapply(dispN,dispFun,configPatches=slPatches,configMatrix=slMatrix,
                             speciesGroup="edge",sens=0.1,dist=500,FEEA=100000)

edgeSLLoSenseDisp500

edgeSLLoSenseDisp500<-do.call("rbind",edgeSLLoSenseDisp500) # have this

edgeSLLoSenseDisp500DF<-data.frame(edgeSLLoSenseDisp500)





#################################################

#int/SL

intSLLoSenseDisp500<-lapply(dispN,dispFun,configPatches=slPatches,configMatrix=slMatrix,
                            speciesGroup="interior",sens=0.1,dist=500,FEEA=100000)


intSLLoSenseDisp500<-do.call("rbind",intSLLoSenseDisp500) # have this

intSLLoSenseDisp500

intSLLoSenseDisp500DF<-data.frame(intSLLoSenseDisp500)

intSLLoSenseDisp500DF

###########################################
##############################################
#################HIGH SENSITIVITY TESTS######################

#edge/mixed
edgeMixedHiSenseDisp500<-lapply(dispN,dispFun,configPatches=mixedPatches,configMatrix=mixedMatrix,
                                speciesGroup="edge",sens=0.9,dist=500,FEEA=100000)

edgeMixedHiSenseDisp500

edgeMixedHiSenseDisp500<-do.call("rbind",edgeMixedHiSenseDisp500) # have this

edgeMixedHiSenseDisp500DF<-data.frame(edgeMixedHiSenseDisp500)

head(edgeMixedHiSenseDisp500DF)



###########################################################

#int/mixed
intMixedHiSenseDisp500<-lapply(dispN,dispFun,configPatches=mixedPatches,configMatrix=mixedMatrix,
                               speciesGroup="interior",sens=0.9,dist=500,FEEA=100000)

intMixedHiSenseDisp500

intMixedHiSenseDisp500<-do.call("rbind",intMixedHiSenseDisp500) # have this

intMixedHiSenseDisp500DF<-data.frame(intMixedHiSenseDisp500)


#int/SS

intSSHiSenseDisp500<-lapply(dispN,dispFun,configPatches=ssPatches,configMatrix=ssMatrix,
                            speciesGroup="interior",sens=0.9,dist=500,FEEA=100000)

intSSHiSenseDisp500

intSSHiSenseDisp500<-do.call("rbind",intSSHiSenseDisp500) # have this

intSSHiSenseDisp500DF<-data.frame(intSSHiSenseDisp500)

intSSHiSenseDisp500DF



##################################################################

#edge/SS

edgeSSHiSenseDisp500<-lapply(dispN,dispFun,configPatches=ssPatches,configMatrix=ssMatrix,
                             speciesGroup="edge",sens=0.9,dist=500,FEEA=100000)

edgeSSHiSenseDisp500

edgeSSHiSenseDisp500<-do.call("rbind",edgeSSHiSenseDisp500) # have this

edgeSSHiSenseDisp500DF<-data.frame(edgeSSHiSenseDisp500)




###############################################################

#int/SL

intSLHiSenseDisp500<-lapply(dispN,dispFun,configPatches=slPatches,configMatrix=slMatrix,
                            speciesGroup="interior",sens=0.9,dist=500,FEEA=100000)

intSLHiSenseDisp500

intSLHiSenseDisp500<-do.call("rbind",intSLHiSenseDisp500) # have this

intSLHiSenseDisp500DF<-data.frame(intSLHiSenseDisp500)
intSLHiSenseDisp500DF

########################################################

#edge/SL

edgeSLHiSenseDisp500<-lapply(dispN,dispFun,configPatches=slPatches,configMatrix=slMatrix,
                             speciesGroup="edge",sens=0.9,dist=500,FEEA=100000)

edgeSLHiSenseDisp500

edgeSLHiSenseDisp500<-do.call("rbind",edgeSLHiSenseDisp500) # have this

edgeSLHiSenseDisp500DF<-data.frame(edgeSLHiSenseDisp500)

edgeSLHiSenseDisp500DF


#################################################Generalists

#GEN/SL

genSLDisp500<-lapply(dispN,dispFun,configPatches=slPatches,configMatrix=slMatrix,
                     speciesGroup="generalist",sens=0.1,dist=500,FEEA=100000)

genSLDisp500

genSLDisp500<-do.call("rbind",genSLDisp500) # have this

genSLDisp500DF<-data.frame(genSLDisp500)

################################################

#GEN/SS

genSSDisp500<-lapply(dispN,dispFun,configPatches=ssPatches,configMatrix=ssMatrix,
                     speciesGroup="generalist",sens=0.1,dist=500,FEEA=100000)

genSSDisp500

genSSDisp500<-do.call("rbind",genSSDisp500) # have this

genSSDisp500DF<-data.frame(genSSDisp500)

#################################################################

#GEN/mixed

genMixedDisp500<-lapply(dispN,dispFun,configPatches=mixedPatches,configMatrix=mixedMatrix,
                        speciesGroup="generalist",sens=0.1,dist=500,FEEA=100000)

genMixedDisp500

genMixedDisp500<-do.call("rbind",genMixedDisp500) # have this

genMixedDisp500DF<-data.frame(genMixedDisp500)




########################################################### CREATE DATA FRAME OF RESULTS


figData<-as.data.frame(cbind(dispN,edgeMixedLoSenseDisp500DF[,3],intMixedLoSenseDisp500DF[,3],
               intSSLoSenseDisp500DF[,3],edgeSSLoSenseDisp500DF[,3],
               edgeSLLoSenseDisp500DF[,3],edgeMixedHiSenseDisp500DF[,3],
               intMixedHiSenseDisp500DF[,3],intSSHiSenseDisp500DF[,3],edgeSSHiSenseDisp500DF[,3],
               edgeSLHiSenseDisp500DF[,3],intSLHiSenseDisp500DF[,3],intSLLoSenseDisp500DF[,3],
               genSLDisp500DF[,3],genSSDisp500DF[,3],genMixedDisp500DF[,3]))

figData
####################


#ENSURE NUMERIC
figData<-apply(figData,MARGIN=2,FUN=unlist)

figData<-as.data.frame(figData)

colnames(figData)<-c("dispersal_Rate","edgeMixedLo","intMixedLo","intSSLo","edgeSSLo",
                     "edgeSLLo","edgeMixedHi","intMixedHi","intSSHi","edgeSSHi","edgeSLHi",
                     "intSLHi","intSLLo","genSL","genSS","genMixed")



library(ggplot2)
library(ggpubr)

###############Low Sensitivity Interior models

loIntPlot<- ggplot(figData)+geom_line(aes(dispersal_Rate,intSSLo,colour="SS"),linewidth=2)+
  
  geom_line(aes(dispersal_Rate,intMixedLo,colour="Mixed"),linewidth=2)+
  
  geom_line(aes(dispersal_Rate,intSLLo,colour="SL"),linewidth=2)+
  
  xlab("Dispersal Probability at maxDist")+ylab("RHI")+scale_color_manual(name = "Landscape", 
                                                        
                                                        values = c("SL" = "darkblue", "SS" = "red","Mixed"="black"))#+theme(text = element_text(size = 20))   



loIntPlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)



######################################################################################################
###################################High Sensitivity Interior Models


hiIntPlot<- ggplot(figData)+geom_line(aes(dispersal_Rate,intSSHi,colour="SS"),linewidth=2)+
  
  geom_line(aes(dispersal_Rate,intMixedHi,colour="Mixed"),linewidth=2)+
  
  geom_line(aes(dispersal_Rate,intSLHi,colour="SL"),linewidth=2)+
  
  xlab("Dispersal Probability at maxDist")+ylab("RHI")+scale_color_manual(name = "Landscape", 
                                                                          
                                                                          values = c("SL" = "darkblue", "SS" = "red","Mixed"="black"))#+theme(text = element_text(size = 20))   



hiIntPlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)

######################################################################################################
##############################################High Sensitivity Edge Models

hiEdgePlot<- ggplot(figData)+geom_line(aes(dispersal_Rate,edgeSSHi,colour="SS"),linewidth=2)+
  
  geom_line(aes(dispersal_Rate,edgeMixedHi,colour="Mixed"),linewidth=2)+
  
  geom_line(aes(dispersal_Rate,edgeSLHi,colour="SL"),linewidth=2)+
  
  xlab("Dispersal Probability at maxDist")+ylab("RHI")+scale_color_manual(name = "Landscape", 
                                                                          
                                                                          values = c("SL" = "darkblue", "SS" = "red","Mixed"="black"))#+theme(text = element_text(size = 20))   



hiEdgePlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)


############################################################################################################
##################################Low Sensitivity Edge Models
loEdgePlot<- ggplot(figData)+geom_line(aes(dispersal_Rate,edgeSSLo,colour="SS"),linewidth=2)+
  
  geom_line(aes(dispersal_Rate,edgeMixedLo,colour="Mixed"),linewidth=2)+
  
  geom_line(aes(dispersal_Rate,edgeSLLo,colour="SL"),linewidth=2)+
  
  xlab("Dispersal Probability at maxDist")+ylab("RHI")+scale_color_manual(name = "Landscape", 
                                                                          
                                                                          values = c("SL" = "darkblue", "SS" = "red","Mixed"="black"))#+theme(text = element_text(size = 20))   



loEdgePlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)


#############################################################
############################################GENERALISTS
genPlot<- ggplot(figData)+geom_line(aes(dispersal_Rate,genSS,colour="SS"),linewidth=2)+
  
  geom_line(aes(dispersal_Rate,genMixed,colour="Mixed"),linewidth=2)+
  
  geom_line(aes(dispersal_Rate,genSL,colour="SL"),linewidth=2)+
  
  xlab("Dispersal Probability at maxDist")+ylab("RHI")+scale_color_manual(name = "Landscape", 
                                                                          
                                                                          values = c("SL" = "darkblue", "SS" = "red","Mixed"="black"))#+theme(text = element_text(size = 20))   



genPlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)





