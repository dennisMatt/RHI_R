setwd("C:/colneValley/colneBig")

load("ws280223.RData")

############################################SL####################################################
install.packages(c("dplyr","sf","terra","igraph","stars","ggplot2","ggpubr"))
library(sf)
library(lwgeom)
library(stars)
library(starsExtra)

setwd("C:/colneValley/colneBig")

SL<-st_read("gridCon_10Ha_SelFin.shp")
SS<-st_read("gridCon_5Ha_SelFin.shp")

st_area(SL)
st_area(ss)

plot(ss$geometry)
mixed<-st_read("gridCon_mixed_10_5.shp")
mixed2<-st_read("gridCon_mixed_10_5v2.shp")
plot(mixed$geometry)
plot(SL$geometry)
plot(SS$geometry)
sum(st_perimeter(SL))
sum(st_perimeter(SS))
sum(st_perimeter(mixed))
sum(st_perimeter(mixed2))

(st_area(SL))#check area of SL patches
(st_area(SS))#check area of SS patches
(st_area(mixed))

mean(st_distance(SL))
mean(st_distance(SS))
mean(st_distance(mixed))
mean(st_distance(mixed2))

library(dplyr)
library(terra)
library(igraph)
library(stars)

patches<-SS

#######
#######This function has six arguments - patches = habitat patches: shapefile 
######################################## specialism = one of "interior", "edge" or "generalist"
################################### edge = the size of the edge effect (for a neutral landscape, will have to update this to handle "edge rasters" in real landscapes)
#################################### edgeComponent = a component (suggest between 0.01 -1) that determines the shape of the kernel (distance decay of edge effect)
#################################### dispersalDist = mean dispersal distance of the species being modeled
#################################### habComponent = component to determine the sensitivity of the species to patch size (i.e. strength of relationship between extinction and patch size) 
#####################################colonizationRate = component setting colonization/survival success
#####################################edgeSensitivity = sets how sensitive the species is to "specialism". Setting this to very small values for either edge or interior basically renders it a generalist
rhiCost<-function(patches,specialism,edge,edgeIntensity,maxDist,dispersalRate,edgeSensitivity,patchCost){ 
  
  
  r<-st_rasterize(patches) # need to rasterize the patches to get at edge kernel and area
  r<-rast(r)#convert to terra library, it's better.
  
  r<-disagg(r,fact=2)# get data to some resolution to reliably establish relative area later on will need to revisit this for using actual patch polygons
  
  
  Clump<-patches(r,directions=4) # get patches (this will account for adjacent cells forming one patch, not in hypothetical landscape but later on when modelling real ones)
  
  
  #####################this next bit is to just buffer the whole study area slightly so patches on the edge still get some edge effect
  
  extentNew<-ext(patches) #get extent
  
  
  ext<-as.polygons(extentNew,crs=crs(patches)) #polygonize the extent
  
  
  extBuff<-terra::buffer(ext,width=100) #use some nominal buffer distance (the size doesn't affect the results as is)
  
  
  
  Clump<-extend(Clump,extBuff) # extend the data area
  
  
  Clump[is.na(Clump)] <-0 #set NA value to use when making the the distance raster later on
  
  distClump<-gridDist(Clump$id_lyr.1,target=0)
  
  distEdge=log(edgeIntensity)/edge# this sets the edge distance and the shape of neg. exp. curve
  
  edgeRast<-exp(distEdge*distClump) # this creates an edge raster where 1 = full edge effect i.e. remove that cell(this happens below) and 0 = no edge effect
  plot(edgeRast)
  invEdge<-1-edgeRast
  
  
  #Need to select species type. 
  if(specialism=="interior"){
    
  edgeClump<-1-(edgeRast*edgeSensitivity)} #Interior removes all edge cell values from patch cells (i.e. removes a proportion of the area)
  
  else if(specialism=="edge"){edgeClump<- 1-(invEdge*edgeSensitivity)} #For edge just leave as the edge habitat within the patch i.e. inverse of above
  
  else if(specialism=="generalist"){edgeClump=distClump>0}# for generalist all patch cells are one i.e. area not changed
  
  #plot(allRast)
  
  extClump<-extract(edgeClump, patches, fun=sum, method="simple", bind=T)
  plot(edgeClump)
  head(extClump)
  # get cell area
  cellArea<-res(edgeClump)[1]^2 
  
  extClump$areaMod<-extClump$id_lyr.1*cellArea #create a vector for the new area
  
  
  
  
  if(patchCost=="FULL"){
    costRaster<-edgeRast
    costRaster[costRaster==1]=2
  costRaster[costRaster<2]=1}
  
  else if(patchCost=="EDGE"){
    costRaster<-edgeRast
    costRaster[costRaster==1]=2}
  
  
  land_cost <- transition(raster(costRaster), transitionFunction=function(x) 1 / mean(x), 8)
  
  
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
  
  if(nrow(distMat)>1){ # in case modelling on single patch then need if else to side step matrix error here
    
  distMat<-apply(distMat,MARGIN=1,FUN=as.numeric)}else{distMat=distMat} # annoyingly need to make sure all elements are numeric here
  
  alpha= -log(dispersalRate)/maxDist# set alpha which determines colonization probability of the species 
  
  if(nrow(distMat)==1){A.prob<- matrix(1, nrow=nrow(distMat), ncol=ncol(distMat)) # for single patches just set diagonal to 1
  
  A.prob<-as.matrix(A.prob)
  #A.prob
  graph.Aprob <- graph.adjacency(A.prob, mode="undirected", weighted=T)}
  
  else{
    
    
    A.prob <- matrix(0, nrow=nrow(distMat), ncol=ncol(distMat))
    A.prob <- exp(-alpha*distMat) #negative exponential of colonization kernel x distance get probability
    
    diag(A.prob) <- 0 # set diag to zero - but this is another parameter that could be tweaked?
    
    A.prob<-as.matrix(A.prob) #final matrix for connectivity graph
    
    
    graph.Aprob <- graph.adjacency(A.prob, mode="undirected", weighted=T) #create adjacency graph
    
  }
  
  
  #calculate landscape-scale connectivity
  pstar.mat <- shortest.paths(graph.Aprob, weights= -log(E(graph.Aprob)$weight)) #calculate all shortest paths between nodes
  pstar.mat <- exp(-pstar.mat)                                                   #back-transform to probabilities of connectedness
  pMatSum<-sum(pstar.mat)
  AL<-expanse(extBuff)# get study area in m2
  
  
  area<-extClump$areaMod # get area vector for PC metric
  
  
  
  
  
  PCmat <- outer(area,area)*pstar.mat #get product of all patch areas ij (to the power of something less that 1) and multiply by probabilities above.
  
  pcMatSum<-sum(PCmat)
  
  
    
    
    pcMod <- pcMatSum/as.numeric(AL^2) #divide by total area of the study squared to get the PC metric  
    
  
  
  
  
  return(pcMod)
}

#set vector of maxDist values

dispVec<-seq(500,10000,by=500)


#Function

dispFun<-function(x,cost,config,speciesGroup){
  
  costRHI<-rhiCost(patches = config,specialism = speciesGroup,edge = 100,patchCost = cost,edgeIntensity = 0.5,maxDist = x,dispersalRate = 0.05,edgeSensitivity = 0.5)
  print(x)
  
  return(costRHI)

}

################################Edge species Cost Tests

#species=EDGE config= SS patchCost=edgeCost

dispCostEdgeSS_EDGE<-lapply(dispVec,dispFun,config=SS,speciesGroup="edge",cost="EDGE")


dispCostEdgeSS_EDGE<-do.call("rbind",dispCostEdgeSS_EDGE)
dispCostEdgeSS_EDGEDF<-data.frame(dispCostEdgeSS_EDGE)



####species=EDGE config= SS patchCost=FULL
dispCostEdgeSS_FULL<-lapply(dispVec,dispFun,config=SS,speciesGroup="edge",cost="FULL")


dispCostEdgeSS_FULL<-do.call("rbind",dispCostEdgeSS_FULL)
dispCostEdgeSS_FULLDF<-data.frame(dispCostEdgeSS_FULL)



#species=EDGE config= SL patchCost=edgeCost

dispCostEdgeSL_EDGE<-lapply(dispVec,dispFun,config=SL,speciesGroup="edge",cost="EDGE")

dispCostEdgeSL_EDGE<-do.call("rbind",dispCostEdgeSL_EDGE)
dispCostEdgeSL_EDGEDF<-data.frame(dispCostEdgeSL_EDGE)



####species=EDGE config= SS patchCost=FULL
dispCostEdgeSL_FULL<-lapply(dispVec,dispFun,config=SL,speciesGroup="edge",cost="FULL")


dispCostEdgeSL_FULL<-do.call("rbind",dispCostEdgeSL_FULL)
dispCostEdgeSL_FULLDF<-data.frame(dispCostEdgeSL_FULL)


#####################################Interior species

#species=INT config= SS patchCost=edgeCost

dispCostINTSS_EDGE<-lapply(dispVec,dispFun,config=SS,speciesGroup="interior",cost="EDGE")


dispCostINTSS_EDGE<-do.call("rbind",dispCostINTSS_EDGE)
dispCostINTSS_EDGEDF<-data.frame(dispCostINTSS_EDGE)



####species=INT config= SS patchCost=FULL
dispCostINTSS_FULL<-lapply(dispVec,dispFun,config=SS,speciesGroup="interior",cost="FULL")


dispCostINTSS_FULL<-do.call("rbind",dispCostINTSS_FULL)
dispCostINTSS_FULLDF<-data.frame(dispCostINTSS_FULL)

dispCostINTSS_FULLDF


#species=INT config= SL patchCost=edgeCost

dispCostINTSL_EDGE<-lapply(dispVec,dispFun,config=SL,speciesGroup="interior",cost="EDGE")

dispCostINTSL_EDGE<-do.call("rbind",dispCostINTSL_EDGE)
dispCostINTSL_EDGEDF<-data.frame(dispCostINTSL_EDGE)



####species=INT config= SS patchCost=FULL
dispCostINTSL_FULL<-lapply(dispVec,dispFun,config=SL,speciesGroup="interior",cost="FULL")


dispCostINTSL_FULL<-do.call("rbind",dispCostINTSL_FULL)
dispCostINTSL_FULLDF<-data.frame(dispCostINTSL_FULL)



dispCostINTSL_FULLDF

#########

costDataAll<-data.frame(cbind(dispVec,dispCostINTSS_FULLDF,dispCostINTSL_FULLDF,dispCostEdgeSS_EDGEDF,
                              dispCostEdgeSL_EDGEDF,dispCostEdgeSL_FULLDF,dispCostEdgeSS_FULLDF,
                              dispCostINTSS_EDGEDF,dispCostINTSL_EDGEDF))

colnames(costDataAll)<-c("Dispersal Distance","costSS_Int_full","costSL_Int_full",
                         "costSS_Edge_edge","costSL_Edge_edge","costSL_Edge_full","costSS_Edge_full",
                         "costSS_Int_edge","costSL_Int_edge")

head(costDataAll)

costDataAll<-apply(costDataAll,MARGIN = 2,FUN=unlist)

costDataAll<-data.frame(costDataAll)

costDataAll$ssDiffInt<-(costDataAll$costSS_Int_edge -costDataAll$costSS_Int_full)/costDataAll$costSS_Int_full*100

costDataAll$slDiffInt<-(costDataAll$costSL_Int_edge -costDataAll$costSL_Int_full)/costDataAll$costSL_Int_full*100

head(costDataAll)

costDataAll$ssDiffEdge<-(costDataAll$costSS_Edge_edge -costDataAll$costSS_Edge_full)/costDataAll$costSS_Edge_full*100

costDataAll$slDiffEdge<-(costDataAll$costSL_Edge_edge -costDataAll$costSL_Edge_full)/costDataAll$costSL_Edge_full*100



costDataAll



#allCostData<-cbind(costDataEdge,costDataInt[,2:7])


####################################################

intPlot<- ggplot(costDataAll)+geom_line(aes(dispVec,slDiffInt,colour="SL"),linewidth=2)+
  
  geom_line(aes(dispVec,ssDiffInt,colour="SS"),linewidth=2)+
  
  
  xlab("Dispersal distance")+ylab("% Change in RHI")+scale_color_manual(name = "Landscape", 
                                                            
                                                            values = c("SL" = "darkblue", "SS" = "darkred"))#+theme(text = element_text(size = 20))   



intPlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)


edgePlot<- ggplot(costDataAll)+geom_line(aes(dispVec,slDiffEdge,colour="SL"),linewidth=2)+
  
  geom_line(aes(dispVec,ssDiffEdge,colour="SS"),linewidth=2)+
  
  
  xlab("Dispersal distance")+ylab("% Change in RHI")+scale_color_manual(name = "Landscape", 
                                                                        
                                                                        values = c("SL" = "darkblue", "SS" = "darkred"))#+theme(text = element_text(size = 20))   


  

edgePlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)



###################################

head(allCostData)

edgePlot<- ggplot(allCostData)+geom_line(aes(Dispersal_Distance,diffSL_edge,colour="SL"),linewidth=1.5)+
  
  geom_line(aes(Dispersal_Distance,diffSS_edge,colour="SS"),linewidth=1.5)+
  
  
  xlab("Dispersal distance")+ylab("RHI")+scale_color_manual(name = "Species group", 
                                                            
                                                            values = c("SL" = "darkblue", "SS" = "darkred"))#+theme(text = element_text(size = 20))   



edgePlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)



###############################
dispEdgeCostSS<-lapply(dispVec,dispFun)

dispAllEdge

dispEdgeCostSS<-do.call("rbind",dispEdgeCostSS)

dispEdgeCostSS<-data.frame(dispEdgeCostSS)

dispEdgeCostSS


########################

dispAllEdge<-lapply(dispVec,dispFun)

dispAllEdge

dispAllEdge<-do.call("rbind",dispAllEdge)

dispAllEdgeSS<-data.frame(dispAllEdge)

(dispAllEdgeSS)

##################

dispAllEdgeSL<-lapply(dispVec,dispFun)

dispAllEdgeSL

dispAllEdgeSL<-do.call("rbind",dispAllEdgeSL)

dispAllEdgeSL<-data.frame(dispAllEdgeSL)


####################################

dispEdgeSS<-lapply(dispVec,dispFun)

dispEdgeSS

dispEdgeSS<-do.call("rbind",dispEdgeSS)

dispEdgeSS<-data.frame(dispEdgeSS)

dispEdgeSS

dispAllEdgeSS


head(dispAllEdge)

plot(dispVec,dispAllEdge$X3)

res(allRast)

dispE2E_edgeSL

dispData<-cbind(dispVec,dispE2E_edgeSLDF[,3],dispE2E_edgeSL_CentroidDF[,3],dispE2E_edgeSS_CentroidDF[,3],dispE2E_edgeSSDF[,3],
                dispE2E_intSLDF[,3],dispE2E_intSL_CentroidDF[,3],dispE2E_intSSDF[,3],dispE2E_intSS_CentroidDF[,3])


colnames(dispData)<-c("Dispersal_distance","edgeSL","edgeSL_C","edgeSS_C","edgeSS","intSL","intSL_C","intSS","intSS_C")


dispData<-apply(dispData,MARGIN = 2,unlist)


dispData<-data.frame(dispData)

dispData

intPlot<- ggplot(dispData)+geom_line(aes(Dispersal_distance,intSlDiff,colour="SL"),linewidth=2)+
  
  geom_line(aes(Dispersal_distance,intssDiff,colour="SS"),linewidth=2)+

  
  xlab("Dispersal distance")+ylab("RHI")+scale_color_manual(name = "Species group", 
                                                        
                                                        values = c("SL" = "darkblue", "SS" = "darkred"))#+theme(text = element_text(size = 20))   



intPlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)

################################################################################

dispData

edgePlot<- ggplot(dispData)+geom_line(aes(Dispersal_distance,edgeSlDiff,colour="SL"),linewidth=2)+
  
  geom_line(aes(Dispersal_distance,edgessDiff,colour="SS"),linewidth=2)+
  
  xlab("Dispersal distance")+ylab("% change in RHI")+scale_color_manual(name = "Species group", 
                                                            
                                                            values = c("SL" = "darkblue", "SS" = "darkred"))#+theme(text = element_text(size = 20))   



edgePlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)



dispData$edgeSlDiff<-(dispData$edgeSL-dispData$edgeSL_C)/dispData$edgeSL_C*100
dispData$edgessDiff<-(dispData$edgeSS-dispData$edgeSS_C)/dispData$edgeSS_C*100


dispData$intSlDiff<-(dispData$intSL-dispData$intSL_C)/dispData$intSL_C*100
dispData$intssDiff<-(dispData$intSS-dispData$intSS_C)/dispData$intSS_C*100



dispData

###################################################################################################################


dispE2E_edgeSS_Centroid<-lapply(dispVec,dispFun)

dispE2E_edgeSS_Centroid<-do.call("rbind",dispE2E_edgeSS_Centroid)

dispE2E_edgeSS_CentroidDF<-data.frame(dispE2E_edgeSS_Centroid)

dispE2E_edgeSS_CentroidDF
################################################################

dispE2E_edgeSL_Centroid<-lapply(dispVec,dispFun)

dispE2E_edgeSL_Centroid<-do.call("rbind",dispE2E_edgeSL_Centroid)

dispE2E_edgeSL_CentroidDF<-data.frame(dispE2E_edgeSL_Centroid)

dispE2E_edgeSL_CentroidDF

######################

dispE2E_intSL_Centroid<-lapply(dispVec,dispFun)

dispE2E_intSL_Centroid<-do.call("rbind",dispE2E_intSL_Centroid)

dispE2E_intSL_CentroidDF<-data.frame(dispE2E_intSL_Centroid)

dispE2E_intSL_CentroidDF
################################################

dispE2E_intSS_Centroid<-lapply(dispVec,dispFun)

dispE2E_intSS_Centroid<-do.call("rbind",dispE2E_intSS_Centroid)

dispE2E_intSS_CentroidDF<-data.frame(dispE2E_intSS_Centroid)

dispE2E_intSS_CentroidDF


####################################
dispE2E_edgeSL<-lapply(dispVec,dispFun)

dispE2E_edgeSL<-do.call("rbind",dispE2E_edgeSL)

dispE2E_edgeSLDF<-data.frame(dispE2E_edgeSL)

dispE2E_edgeSLDF

############################################################

dispE2E_intSL<-lapply(dispVec,dispFun)

dispE2E_intSL<-do.call("rbind",dispE2E_intSL)

dispE2E_intSLDF<-data.frame(dispE2E_intSL)

################################################

dispE2E_intSS<-lapply(dispVec,dispFun)

dispE2E_intSS<-do.call("rbind",dispE2E_intSS)

dispE2E_intSSDF<-data.frame(dispE2E_intSS)





##########################################################################

dispE2E_edgeSS<-lapply(dispVec,dispFun)

dispE2E_edgeSS<-do.call("rbind",dispE2E_edgeSS)

dispE2E_edgeSSDF<-data.frame(dispE2E_edgeSS)

dispE2E_edgeSSDF

######################

dispE2E_edgeSL<-lapply(dispVec,dispFun)

dispE2E_edgeSL<-do.call("rbind",dispE2E_edgeSL)

dispE2E_edgeSLDF<-data.frame(dispE2E_edgeSL)

dispE2E_edgeSLDF

############################################################

dispE2E_intSL<-lapply(dispVec,dispFun)

dispE2E_intSL<-do.call("rbind",dispE2E_intSL)

dispE2E_intSLDF<-data.frame(dispE2E_intSL)

################################################

dispE2E_intSS<-lapply(dispVec,dispFun)

dispE2E_intSS<-do.call("rbind",dispE2E_intSS)

dispE2E_intSSDF<-data.frame(dispE2E_intSS)




#######################################################


intensityIntSL_CentroidDF<-cbind(sampN,intensityIntSL_CentroidDF)

intensityIntSL_CentroidDF

intIntSLDF<-intensityIntSL_CentroidDF[,c(1,6)]

plot(intIntSLDF)

##################################################
intesnsityEdgeSL_Centroid<-lapply(sampN,dispFun)



intensityEdgeSL_Centroid<-do.call("rbind", intesnsityEdgeSL_Centroid)

intensityEdgeSL_CentroidDF<-data.frame(intesnsityEdgeSL_Centroid)

intensityEdgeSL_CentroidDF<-cbind(sampN,intensityEdgeSL_CentroidDF)

(intensityEdgeSL_CentroidDF)

edgeIntSLDF<-intensityEdgeSL_CentroidDF[,c(1,8)]

#intesnsityEdgeSL_CentroidDF<-apply(dispDF,MARGIN = 2,unlist)  


slIntData<-cbind(intIntSLDF,edgeIntSLDF)

slIntData<-slIntData[,c(1,2,4)]
slIntData

colnames(slIntData)<-c("edgeIntensity","RHI_intSL","RHI_edgeSL")

plot(slIntData)

############################################################################SS



intensityIntSS<-lapply(sampN,dispFun)



intensityIntSS<-do.call("rbind", intensityIntSS)

intensityIntSS<-data.frame(intensityIntSS)

intensityIntSSDF<-cbind(sampN,intensityIntSS)

intensityIntSSDF



intIntSSDF<-intensityIntSSDF[,c(1,6)]

plot(intIntSLDF)

##################################################

intensityEdgeSS<-lapply(sampN,dispFun)



intensityEdgeSS<-do.call("rbind", intensityEdgeSS)

intensityEdgeSS<-data.frame(intesnsityEdgeSS)

intensityEdgeSS<-cbind(sampN,intensityEdgeSS)

(intensityEdgeSS[,c(1,6)])

edgeIntSSDF<-intensityEdgeSS[,c(1,6)]

#intesnsityEdgeSL_CentroidDF<-apply(dispDF,MARGIN = 2,unlist)  


ssIntData<-cbind(intIntSSDF,edgeIntSSDF)

ssIntData<-ssIntData[,c(1,2,4)]
ssIntData

colnames(ssIntData)<-c("edgeIntensity","RHI_intSS","RHI_edgeSS")

plot(ssIntData$edgeIntensity,ssIntData$RHI_intSS)

############################################################################Mixed

plot(mixed)

dispFun<-function(x){
  
  testConnect<-castorConnect(patches = mixed,specialism = "edge",edge = 100,edgeComponent = x,dispersalDist = 5000,minPatchSize = 1, habComponent = 0.5,colonizationRate = 0.05,perim=FALSE,edgeSensitivity = 0.5)
  print(x)
  
  return(testConnect)
  
}

intensityIntMixed<-lapply(sampN,dispFun)



intensityIntMixed<-do.call("rbind", intensityIntMixed)

intensityIntMixed<-data.frame(intensityIntMixed)

intensityIntMixed<-cbind(sampN,intensityIntMixed)

intensityIntMixed



intIntMixedDF<-intensityIntMixed[,c(1,4)]

plot(intIntMixedDF)

##################################################

intensityEdgeMixed<-lapply(sampN,dispFun)



intensityEdgeMixed<-do.call("rbind", intensityEdgeMixed)

intensityEdgeMixed<-data.frame(intensityEdgeMixed)

intensityEdgeMixed<-cbind(sampN,intensityEdgeMixed)

intensityEdgeMixed

(intensityEdgeMixed[,c(1,4)])

edgeIntMixedDF<-intensityEdgeMixed[,c(1,4)]

#intesnsityEdgeSL_CentroidDF<-apply(dispDF,MARGIN = 2,unlist)  

edgeIntMixedDF

mixedIntData<-cbind(intIntMixedDF,edgeIntMixedDF)

mixedIntData<-mixedIntData[,c(1,2,4)]
plot(mixedIntData)

colnames(mixedIntData)<-c("edgeIntensity","RHI_intMixed","RHI_edgeMixed")



allIntData<-cbind(slIntData,ssIntData,mixedIntData)


############################################################################
head(allIntData)

allIntData<-allIntData[,c(1,2,3,5,6,8,9)]

class(allIntData[,3])

allIntData<-apply(allIntData,MARGIN = 2,unlist)

allIntData<-data.frame(allIntData)

library(ggplot2)
library(ggpubr)
###################################################################################

slPlot<- ggplot(allIntData)+geom_line(aes(edgeIntensity,RHI_intSL,colour="interior"),linewidth=2)+
  
  geom_line(aes(edgeIntensity,RHI_edgeSL,colour="edge"),linewidth=2)+
  
  xlab("edge Intensity")+ylab("RHI")+scale_color_manual(name = "Species group", values = c("interior" = "darkblue", "edge" = "red"))


slPlot+ theme_pubr()

ssPlot<- ggplot(allIntData)+geom_line(aes(edgeIntensity,RHI_intSS,colour="interior"),linewidth=2)+
  
  geom_line(aes(edgeIntensity,RHI_edgeSS,colour="edge"),linewidth=2)+
  
  xlab("edge Intensity")+ylab("RHI")+scale_color_manual(name = "Species group", 
                                                        
  values = c("interior" = "darkblue", "edge" = "red"))#+theme(text = element_text(size = 20))   

head(allIntData)

mixedPlot<- ggplot(allIntData)+geom_line(aes(edgeIntensity,RHI_intMixed,colour="interior"),linewidth=2)+
  
  geom_line(aes(edgeIntensity,RHI_edgeMixed,colour="edge"),linewidth=2)+
  
  xlab("edge Intensity")+ylab("RHI")+scale_color_manual(name = "Species group", 
                                                        
                                                        values = c("interior" = "darkblue", "edge" = "red"))#+theme(text = element_text(size = 20))   



ssPlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)
 

slPlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)

mixedPlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)


 
ggpar(ssPlot,font.legend = c(20,"bold","black"))
  
ggpar(ssPlot,legend = "right", legend.title = "Dose (mg)",
      font.legend = c(20, "bold", "red"))+theme_pubr(base_size = 15)+labs_pubr(base_size = 15)




head(allIntData)

slDat<-allIntData[,c(1,2,3)]

ssDat<-allIntData[,c(1,4,5)]


colnames(slDat)<-c("edgeIntensity","interior","edge")

colnames(ssDat)<-c("edgeIntensity","interior","edge")


eiSL<-ggline(slDat,x="edgeIntensity",y=c("interior","edge"),merge=TRUE,
       
       plot_type = "l",size=2,ylab = "RHI", ggtheme = theme_pubr())
       

eiSL+font("legend.text",size=15)+font("xlab",size=15)+font("ylab",size = 15)


eiSS<-ggline(ssDat,x="edgeIntensity",y=c("interior","edge"),merge=TRUE,
             
             plot_type = "l",size=2,ylab = "RHI", ggtheme = theme_pubr())


eiSS+font("legend.text",size=15)




   #add = "loess",  # Add regression line
          
add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 1, label.sep = "\n"))


########################################


ggplot(allIntData)+geom_line(aes(edgeIntensity,RHI_intSS,colour="interior"),size=1)+
  
  geom_line(aes(edgeIntensity,RHI_edgeSS,colour="edge"),size=1)+
  
  xlab("edge Intensity")+ylab("RHI")+scale_color_manual(name = "Species group", values = c("interior" = "darkblue", "edge" = "red"))

###########################################################

ggplot(allIntData)+geom_line(aes(edgeIntensity,RHI_intMixed,colour="interior"),size=1)+
  
  geom_line(aes(edgeIntensity,RHI_edgeMixed,colour="edge"),size=1)+
  
  xlab("edge Intensity")+ylab("RHI")+scale_color_manual(name = "Species group", values = c("interior" = "darkblue", "edge" = "red"))


devtools::install_github("tidyverse/ggplot2")
library(ggplot2)
# Install
#if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")


library(ggpubr)




#########################################################################################


####################################################################################
dispVec<-seq(500,20000,by=500)

dispEdgeSL<-lapply(dispVec,dispFun)

dispDFSLEdge_Centroid<-do.call("rbind", dispDataEdgeSL_Centroid)

dispDFSLEdge_Centroid<-data.frame(dispDFSLEdge_Centroid)

dispDFSLEdge_Centroid<-cbind(dispVec,dispDFSLEdge_Centroid)

dispDF<-apply(dispDF,MARGIN = 2,unlist)

head(dispDFEdgeCentroid)

dispDFSL_Edge<-data.frame(dispDF)

colnames(dispDFSLEdge_Centroid)<-c("dispDistSLEdge_Centroid","areaSumSLEdge_Centroid","edgeRatioSLEdge_Centroid","pMatSumSLEdge_Centroid","pcMatSumSLEdge_Centroid","pcModSLEdge_Centroid")

dispDFSLEdge_E2E
dispDFSLInt_E2E
dispDFSSInt_E2E
dispDFSSEdge_E2E
dispDFSSEdge_Centroid
dispDFSSInt_Centroid
dispDFSLInt_Centroid
dispDFSLEdge_Centroid


allDispData<-cbind(dispDFSLEdge_E2E,dispDFSLInt_E2E,dispDFSSInt_E2E,dispDFSSEdge_E2E,dispDFSSEdge_Centroid,dispDFSSInt_Centroid,dispDFSLEdge_Centroid,dispDFSLInt_Centroid)

allDispData<-apply(allDispData,MARGIN = 2,unlist)

allDispData<-data.frame(allDispData)

allDispData$edgeSSDiff=allDispData$pcModSSEdge_e2e-allDispData$pcModSSEdge_Centroid

allDispData$intSSDiff=allDispData$pcModSSInt_e2e-allDispData$pcModSSInt_Centroid

allDispData$edgeSLDiff=allDispData$pcModSLEdge_e2e-allDispData$pcModSLEdge_Centroid

allDispData$intSLDiff=allDispData$pcModSLInt_e2e-allDispData$pcModSLInt_Centroid

allDispData$intSSDiffPC<-allDispData$intSSDiff/allDispData$pcModSSInt_Centroid*100
allDispData$intSLDiffPC<-allDispData$intSLDiff/allDispData$pcModSLInt_Centroid*100
allDispData$edgeSLDiffPC<-allDispData$edgeSLDiff/allDispData$pcModSLEdge_Centroid*100
allDispData$edgeSSDiffPC<-allDispData$edgeSSDiff/allDispData$pcModSSEdge_Centroid*100


head(allDispData)

ggplot(allDispData)+geom_line(aes(dispDistSLEdge_e2e,edgeSSDiffPC,colour="edgeSS"),size=1)+
  geom_line(aes(dispDistSLEdge_e2e,edgeSLDiffPC,colour="edgeSL"),size=1)+xlab("Dispersal Distance")+ylab("% increase in connectivity")+scale_color_manual(name = "Species group", values = c("edgeSS" = "darkblue", "edgeSL" = "red"))


ggplot(allDispData)+geom_line(aes(dispDistSLInt_e2e,intSSDiffPC,colour="interiorSS"),size=1)+
  geom_line(aes(dispDistSLInt_e2e,intSLDiffPC,colour="interiorSL"),size=1)+xlab("Dispersal Distance")+ylab("% increase in connectivity")+scale_color_manual(name = "Species group", values = c("interiorSS" = "darkblue", "interiorSL" = "red"))



###############
head(allDispData)
dispDFIntCentroid
dispDFEdgeCentroid
dispDFSSEdgeCentroid
dispDFSSIntCentroid

dispDFSLInt
dispDFSLEdge

dispSLData<-cbind(dispDFSLInt,dispDFSLEdge)



dispSLData<-apply(dispSLData,MARGIN = 2,unlist)

dispSLData<-data.frame(dispSLData)

allSS<-cbind(dispDFSSIntCentroid,dispSSIntData)

head(allSS)


head(dispSLData)

dispSSIntData
dispSLData

ggplot(dispSLData)+geom_line(aes(dispDistanceSLInt,pcModSLInt,colour="interior"),size=1)+
  geom_line(aes(dispDistanceSLInt,pcModSLEdge,colour="edge"),size=1)+xlab("Dispersal Distance")+ylab("connectivity")+scale_color_manual(name = "Species group", values = c("interior" = "darkblue", "edge" = "red"))




write.csv(dispDFSS_Int,file="dispSS_Int_0_01edgeComp.csv")
write.csv(dispDFSS_Edge,file="dispMixed_edge_0_01edgeComp.csv")


#############################################edge component


sampN<-seq(0.05,1,by=0.05)

edgeFun<-function(x){
  
  testConnect<-castorConnect(patches = SS,specialism = "edge",edge = 100,edgeComponent = x,dispersalDist = 5000,minPatchSize = 0, habComponent = 0.9,colonizationRate = 0.05,perim=TRUE,edgeSensitivity = 0.5)
  print(x)
  #print(area)
  
  return(testConnect)
  
}

#dispVec<-seq(500,20000,by=500)

edgeDataSSIntensityMP<-lapply(sampN,edgeFun)



edgeDFSSIntenseMP<-do.call("rbind",edgeDataSSIntensityMP)
intDFSIntenseMP
edgeDFSSIntenseMP<-data.frame(edgeDFSSIntenseMP)

edgeDFSSIntenseMP<-cbind(sampN,edgeDFSSIntenseMP)
intDFMixedIntenseMP


edgeDF<-apply(edgeDF,MARGIN = 2,unlist)
#intDataSLSensMP



(intDFSLSenseMP)

edgeDFSLSenseMP


##########################################
intDFSSSenseMP

edgeDFSLSenseMP


#sensDFSS_Edge<-data.frame(edgeDF)

edgeDFSLSenseMP<-edgeDFSLSenseMP[,2:7]

colnames(edgeDFSSIntenseMP)<-c("intenseEdgeMP","areaSumEdgeMP","edgeRatioEdgeMP","pMatSumEdgeMP","pcMatSumEdgeMP","pcModEdgeMP")
colnames(intDFSSIntenseMP)<-c("intenseIntMP","areaSumIntMP","edgeRatioIntMP","pMatSumIntMP","pcMatSumIntMP","pcModIntMP")

edgeDFMixedIntenseMP

head(edgeDFSLIntenseMP)
(intDFMixedSenseMP)

IntenseSSData<-cbind(edgeDFSSIntenseMP,intDFSSIntenseMP)

head(mpIntenseData)

IntenseSSData<-apply(IntenseSSData,MARGIN = 2,unlist)

IntenseSSData<-data.frame(IntenseSSData)

IntenseMixedData

ggplot(IntenseSSData)+geom_line(aes(intenseIntMP,pcModIntMP,colour="interior"),size=1)+
  geom_line(aes(intenseEdgeMP,pcModEdgeMP,colour="edge"),size=1)+xlab("Edge intensity")+ylab("connectivity")+scale_color_manual(name = "Species group", values = c("interior" = "darkblue", "edge" = "red"))



####################################################

intDFSSSenseMP
intDFSLSenseMP
intDFMixedSenseMP






ggplot(slSensData)+geom_line(aes(sensitivityInt,pcModIntSL,colour="interior"),size=1)+
  geom_line(aes(sensitivityEdge,pcModEdgeSL,colour="edge"),size=1)+xlab("Edge Sensitivity")+ylab("connectivity")+scale_color_manual(name = "Species group", values = c("interior" = "darkblue", "edge" = "red"))



#############################################
intDataMixedSens<-lapply(sampN,edgeFun)

intDFmixedSense<-do.call("rbind", intDataMixedSens)

intDFmixedSense<-data.frame(intDFmixedSense)

intDFmixedSense<-cbind(sampN,intDFmixedSense)



edgeDF<-apply(edgeDF,MARGIN = 2,unlist)

head(intDFmixedSense)

#sensDFSS_Edge<-data.frame(edgeDF)

colnames(intDFmixedSense)<-c("sensitivity","areaSum","edgeRatio","pMatSum","pcMatSum","pcMod")
#################################################

intDataSLSens<-lapply(sampN,edgeFun)

intDFSLSense<-do.call("rbind", intDataSLSens)

intDFSLSense<-data.frame(intDFSLSense)

intDFSLSense<-cbind(sampN,intDFSLSense)





head(intDFSLSense)

#sensDFSS_Edge<-data.frame(edgeDF)

colnames(intDFSLSense)<-c("sensitivity","areaSum","edgeRatio","pMatSum","pcMatSum","pcMod")


###################################################

edgeDataSSSens<-lapply(sampN,edgeFun)

edgeDFSSSense<-do.call("rbind", edgeDataSSSens)

edgeDFSSSense<-data.frame(edgeDFSSSense)

edgeDFSSSense<-cbind(sampN,edgeDFSSSense)



edgeDF<-apply(edgeDF,MARGIN = 2,unlist)

head(edgeDFSSSense)

#sensDFSS_Edge<-data.frame(edgeDF)

colnames(edgeDFSSSense)<-c("sensitivity","areaSum","edgeRatio","pMatSum","pcMatSum","pcMod")

(edgeDFSSSense)
(edgeDFSLSense)
(edgeDFmixedSense)
(edgeDFmixed2Sense)

(intDFSSSense)
(intDFSLSense)
(intDFmixedSense)



library(ggplot2)
class(intDFSLSense$pcMod)

intDFSLSense<-apply(intDFSLSense,MARGIN = 2,unlist)
intDFSLSense<-data.frame(intDFSLSense)
intDFSSSense<-apply(intDFSSSense,MARGIN = 2,unlist)
intDFSSSense<-data.frame(intDFSSSense)




colnames(intDFSLSense)<-c("sensitivityInt","areaSumIntSL","edgeRatioIntSL","pMatSumIntSL","pcMatSumIntSL","pcModIntSL")
colnames(edgeDFSLSense)<-c("sensitivityEdge","areaSumEdgeSL","edgeRatioEdgeSL","pMatSumEdgeSL","pcMatSumEdgeSL","pcModEdgeSL")

slSensData<-cbind(intDFSLSense,edgeDFSLSense)
slSensData<-apply(slSensData,MARGIN = 2,unlist)
slSensData<-data.frame(slSensData)


head(slSensData)
ggplot(slSensData)+geom_line(aes(sensitivityInt,pcModIntSL,colour="interior"),size=1)+
  geom_line(aes(sensitivityEdge,pcModEdgeSL,colour="edge"),size=1)+xlab("Edge Sensitivity")+ylab("connectivity")+scale_color_manual(name = "Species group", values = c("interior" = "darkblue", "edge" = "red"))


colnames(intDFSSSense)<-c("sensitivityInt","areaSumIntSS","edgeRatioIntSS","pMatSumIntSS","pcMatSumIntSS","pcModIntSS")
colnames(edgeDFSSSense)<-c("sensitivityEdge","areaSumEdgeSS","edgeRatioEdgeSS","pMatSumEdgeSS","pcMatSumEdgeSS","pcModEdgeSS")

ssSensData<-cbind(intDFSSSense,edgeDFSSSense)
ssSensData<-apply(ssSensData,MARGIN = 2,unlist)
ssSensData<-data.frame(ssSensData)


head(ssSensData)
ggplot(ssSensData)+geom_line(aes(sensitivityInt,pcModIntSS,col="interior"),size=1)+
  geom_line(aes(sensitivityEdge,pcModEdgeSS,col="edge"),size=1)+xlab("Edge Sensitivity")+ylab("connectivity")+xlab("Edge Sensitivity")+ylab("connectivity")+scale_color_manual(name = "Species group", values = c("interior" = "darkblue", "edge" = "red"))

##########################################

colnames(intDFmixedSense)<-c("sensitivityInt","areaSumIntMixed","edgeRatioIntMixed","pMatSumIntMixed","pcMatSumIntMixed","pcModIntMixed")
colnames(edgeDFmixedSense)<-c("sensitivityEdge","areaSumEdgeMixed","edgeRatioEdgeMixed","pMatSumEdgeMixed","pcMatSumEdgeMixed","pcModEdgeMixed")

mixedSensData<-cbind(intDFmixedSense,edgeDFmixedSense)
mixedSensData<-apply(mixedSensData,MARGIN = 2,unlist)
mixedSensData<-data.frame(mixedSensData)


head(mixedSensData)
ggplot(mixedSensData)+geom_line(aes(sensitivityInt,pcModIntMixed,colour="interior"),size=1)+
  geom_line(aes(sensitivityEdge,pcModEdgeMixed,colour="edge"),size=1)+
  xlab("Edge Sensitivity")+ylab("connectivity")+xlab("Edge Sensitivity")+ylab("connectivity")+xlab("Edge Sensitivity")+ylab("connectivity")+scale_color_manual(name = "Species group", values = c("interior" = "darkblue", "edge" = "red"))





class(edgeDFSSSense)
getwd()

write.csv(edgeDFSLSense,file="edgeSensSS_5000.csv")


plot(edgeDFSLSense$sensitivity, edgeDFSLSense$pcMod)
plot(edgeDFSSSense$sensitivity, edgeDFSSSense$pcMod)


plot(SS$geometry)
plot(SL$geometry)
plot(mixed$geometry)
plot(mixed2$geometry)

###########################################################



pairs(sensDFSS_Int)
pairs(sensDFSS_Edge)

edgeDFSS_Edge
dispDFSS_Int
edgeDFSL_Int

write.csv(sensDFSS_Edge,file="sensSS_Edge_disp5000.csv")
write.csv(edgeDFSL_Edge,file="edgeCompSLedge_disp5000.csv")


##############################################################

dispFun<-function(x){
  
  testConnect<-castorConnect(patches = SS,specialism = "edge",edge = 100,edgeComponent = 0.5,dispersalDist = 5000,minPatchSize = 1, habComponent = 0.5,colonizationRate = 0.05,perim=FALSE,edgeSensitivity = 0.5)
  print(x)
  
  return(testConnect)
  
}




###################################################################################
dispFun<-function(x){
  
  testConnect<-castorConnect(patches = SS,specialism = "edge",edge = 100,edgeComponent = 0.5,dispersalDist = 5000,minPatchSize = 1, habComponent = 0.5,colonizationRate = 0.05,perim=FALSE,edgeSensitivity = 0.5)
  print(x)
  
  return(testConnect)
  
}

dispVec<-seq(1000,20000,by=1000)

dispEdgeMixed<-lapply(sampN,dispFun)

plot(SL)

sensEdgeMixed<-do.call("rbind", sensIntMixed)

sensIntMixedDF<-data.frame(sensIntMixed)

sensIntMixedDF





########################################################################

sensIntSS<-lapply(sampN,dispFun)

plot(SL)

sensIntSS<-do.call("rbind", sensIntSS)

sensIntSSDF<-data.frame(sensIntSS)

sensIntSSDF

plot(sensIntSSDF)

#sensitivityIntSL<-cbind(sampN,sensIntSLDF)




####################################################################

sensEdgeSS<-lapply(sampN,dispFun)

plot(SL)

sensEdgeSS<-do.call("rbind", sensEdgeSS)

sensEdgeSSDF<-data.frame(sensEdgeSS)

sensEdgeSSDF
sensIntSSDF
sampN

sensitivityData<-cbind(sampN,sensIntSSDF[,3],sensEdgeSSDF[,3],sensIntSLDF[,3],
                       sensEdgeSLDF[,3],sensIntMixedDF[,3],sensEdgeMixedDF[,3])

colnames(sensitivityData)<-c("edgeSensitivity","intSS","edgeSS","intSL","edgeSL","intMixed","edgeMixed")

plot(sensitivitySS)

sensitivityData<-apply(sensitivityData,MARGIN = 2,unlist)

sensitivityData<-data.frame(sensitivityData)

class(sensitivityData)

slSensPlot<- ggplot(sensitivityData)+geom_line(aes(edgeSensitivity,intSL,colour="interior"),linewidth=2)+
  
  geom_line(aes(edgeSensitivity,edgeSL,colour="edge"),linewidth=2)+
  
  xlab("edgeSensitivity")+ylab("RHI")+scale_color_manual(name = "Species group", values = c("interior" = "darkblue", "edge" = "red"))




slSensPlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)

###########################################################################################


ssSensPlot<- ggplot(sensitivityData)+geom_line(aes(edgeSensitivity,intSS,colour="interior"),linewidth=2)+
  
  geom_line(aes(edgeSensitivity,edgeSS,colour="edge"),linewidth=2)+
  
  xlab("edgeSensitivity")+ylab("RHI")+scale_color_manual(name = "Species group", values = c("interior" = "darkblue", "edge" = "red"))




ssSensPlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)


#####################################################################################

mixedSensPlot<- ggplot(sensitivityData)+geom_line(aes(edgeSensitivity,intMixed,colour="interior"),linewidth=2)+
  
  geom_line(aes(edgeSensitivity,edgeMixed,colour="edge"),linewidth=2)+
  
  xlab("edgeSensitivity")+ylab("RHI")+scale_color_manual(name = "Species group", values = c("interior" = "darkblue", "edge" = "red"))




mixedSensPlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)





#############################################################################
sensEdgeSL<-lapply(sampN,dispFun)

plot(SL)

sensEdgeSL<-do.call("rbind", sensEdgeSL)

sensEdgeSLDF<-data.frame(sensEdgeSL)

sensEdgeSLDF

sensitivityIntSL<-cbind(sampN,sensIntSLDF)




##################################################################

sensEdgeSL<-lapply(sampN,dispFun)

plot(SL)

sensEdgeSL<-do.call("rbind", sensEdgeSL)

sensEdgeSLDF<-data.frame(sensEdgeSL)

sensEdgeSLDF

sensitivityIntSL<-cbind(sampN,sensIntSLDF)



###############################################
sesnIntSL<-lapply(sampN,dispFun)

plot(SL)

sensIntSL<-do.call("rbind", sesnIntSL)

sensIntSLDF<-data.frame(sensIntSL)

sensIntSLDF

sensitivityIntSL<-cbind(sampN,sensIntSLDF)


sensitivityIntSL<-sensitivityIntSL[,c(1,4)]

plot(sensitivityIntSL)





##################################################

sensEdgeMixed<-lapply(sampN,dispFun)



sensEdgeMixed<-do.call("rbind", sensEdgeMixed)

sensEdgeMixedDF<-data.frame(sensEdgeMixed)

#sensEdgeMixedDF<-cbind(sampN,sensEdgeMixedDF)

sensEdgeMixedDF



allIntData<-cbind(slIntData,ssIntData,mixedIntData)


############################################################################
head(allIntData)

allIntData<-allIntData[,c(1,2,3,5,6,8,9)]

class(allIntData[,3])

allIntData<-apply(allIntData,MARGIN = 2,unlist)

allIntData<-data.frame(allIntData)

library(ggplot2)
library(ggpubr)
###################################################################################

slPlot<- ggplot(allIntData)+geom_line(aes(edgeIntensity,RHI_intSL,colour="interior"),linewidth=2)+
  
  geom_line(aes(edgeIntensity,RHI_edgeSL,colour="edge"),linewidth=2)+
  
  xlab("edge Intensity")+ylab("RHI")+scale_color_manual(name = "Species group", values = c("interior" = "darkblue", "edge" = "red"))


slPlot+ theme_pubr()

ssPlot<- ggplot(allIntData)+geom_line(aes(edgeIntensity,RHI_intSS,colour="interior"),linewidth=2)+
  
  geom_line(aes(edgeIntensity,RHI_edgeSS,colour="edge"),linewidth=2)+
  
  xlab("edge Intensity")+ylab("RHI")+scale_color_manual(name = "Species group", 
                                                        
                                                        values = c("interior" = "darkblue", "edge" = "red"))#+theme(text = element_text(size = 20))   

head(allIntData)

mixedPlot<- ggplot(allIntData)+geom_line(aes(edgeIntensity,RHI_intMixed,colour="interior"),linewidth=2)+
  
  geom_line(aes(edgeIntensity,RHI_edgeMixed,colour="edge"),linewidth=2)+
  
  xlab("edge Intensity")+ylab("RHI")+scale_color_manual(name = "Species group", 
                                                        
                                                        values = c("interior" = "darkblue", "edge" = "red"))#+theme(text = element_text(size = 20))   



ssPlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)


slPlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)

mixedPlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)



ggpar(ssPlot,font.legend = c(20,"bold","black"))

ggpar(ssPlot,legend = "right", legend.title = "Dose (mg)",
      font.legend = c(20, "bold", "red"))+theme_pubr(base_size = 15)+labs_pubr(base_size = 15)




head(allIntData)

slDat<-allIntData[,c(1,2,3)]

ssDat<-allIntData[,c(1,4,5)]


colnames(slDat)<-c("edgeIntensity","interior","edge")

colnames(ssDat)<-c("edgeIntensity","interior","edge")


eiSL<-ggline(slDat,x="edgeIntensity",y=c("interior","edge"),merge=TRUE,
             
             plot_type = "l",size=2,ylab = "RHI", ggtheme = theme_pubr())


eiSL+font("legend.text",size=15)+font("xlab",size=15)+font("ylab",size = 15)


eiSS<-ggline(ssDat,x="edgeIntensity",y=c("interior","edge"),merge=TRUE,
             
             plot_type = "l",size=2,ylab = "RHI", ggtheme = theme_pubr())


eiSS+font("legend.text",size=15)



write.csv(dataFin,file="C:/colnevalley/colneBig/dataFinFin.csv")

head(dataFin)


