

############################################SL####################################################
install.packages(c("sf","terra","igraph","stars","ggplot2","ggpubr"))
library(sf)
library(lwgeom)
library(stars)
library(starsExtra)
library(igraph)
library(ggplot2)
library(ggpubr)

#Read in the data
SL<-st_read("Data/slHomog.shp")
SS<-st_read("Data/ssHomog.shp")


#######
#######This function is a modified version of the Realized Habitat Index base function to test with-patch movement cost
######################################It has eight arguments:
###################################### patches = habitat patches: shapefile 
######################################## specialism = one of "interior", "edge" or "generalist"
################################### edge = the size of the edge effect (for a homogenous landscape)
#################################### edgeIntensity = a component (suggest between 0.01 -1) that determines the shape of the kernel (distance decay of edge effect)
#################################### maxDist = mean dispersal distance of the species being modeled
 
#####################################dispersalRate = component setting colonization/survival success
#####################################edgeSensitivity = sets how sensitive the species is to "specialism". Setting this to very small values for either edge or interior basically renders it a generalist
#######################################patchCost = one of "EDGE" or "FULL" to set within-patch cost to the pixel edge value (range 0-1) or to 1 



rhiCost<-function(patches,specialism,edge,edgeIntensity,maxDist,dispersalRate,edgeSensitivity,patchCost){ 
  
  
  r<-st_rasterize(patches) # need to rasterize the patches to get at edge kernel and area
  r<-rast(r)#convert to terra library.
  
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
  
  #for edge specialists need the inverse of edgeRast to represent habitat
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
  
  
  
  # set within-patch movement cost to either "FULL" ~ Euclidean distance or "EDGE" = edge raster
  if(patchCost=="FULL"){
    costRaster<-edgeRast
    costRaster[costRaster==1]=2
  costRaster[costRaster<2]=1}
  
  else if(patchCost=="EDGE"){
    costRaster<-edgeRast
    costRaster[costRaster==1]=2}
  
  
  #create transition matrix for least cost analysis
  land_cost <- transition(raster(costRaster), transitionFunction=function(x) 1 / mean(x), 8)
  
  # create within-patch centroids (st_point_on_surface function ensure centroids are within the patch 
  # - important for irregular shaped polygons)
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
    
  distMat<-apply(distMat,MARGIN=1,FUN=as.numeric)}else{distMat=distMat} # need to make sure all elements are numeric here
  
  alpha= -log(dispersalRate)/maxDist# set alpha which determines colonization probability of the species 
  
  if(nrow(distMat)==1){A.prob<- matrix(1, nrow=nrow(distMat), ncol=ncol(distMat)) # for single patches just set diagonal to 1
  
  A.prob<-as.matrix(A.prob)
  
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
  
  
  area<-extClump$areaMod # get area vector 
  
  
  
  
  
  PCmat <- outer(area,area)*pstar.mat #get product of all patch areas ij and multiply by probabilities above.
  
  # sum PCmat
  pcMatSum<-sum(PCmat)
  
  
    
    
  RHI <- pcMatSum/as.numeric(AL^2) #divide by total area of the study squared to get the RHI metric  
    
  
  
  
  
  return(RHI)
}

#set vector of maxDist values

dispVec<-seq(500,10000,by=500)


#Function for running RHI based on different within-patch cost-estimation method 

dispFun<-function(x,cost,config,speciesGroup){
  
  costRHI<-rhiCost(patches = config,specialism = speciesGroup,edge = 50,patchCost = cost,edgeIntensity = 0.5,maxDist = x,dispersalRate = 0.05,edgeSensitivity = 0.5)
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


#species=INT config= SL patchCost=edgeCost

dispCostINTSL_EDGE<-lapply(dispVec,dispFun,config=SL,speciesGroup="interior",cost="EDGE")

dispCostINTSL_EDGE<-do.call("rbind",dispCostINTSL_EDGE)
dispCostINTSL_EDGEDF<-data.frame(dispCostINTSL_EDGE)



####species=INT config= SS patchCost=FULL

dispCostINTSL_FULL<-lapply(dispVec,dispFun,config=SL,speciesGroup="interior",cost="FULL")


dispCostINTSL_FULL<-do.call("rbind",dispCostINTSL_FULL)
dispCostINTSL_FULLDF<-data.frame(dispCostINTSL_FULL)


######### Create Data Frame

costDataAll<-data.frame(cbind(dispVec,dispCostINTSS_FULLDF,dispCostINTSL_FULLDF,dispCostEdgeSS_EDGEDF,
                              dispCostEdgeSL_EDGEDF,dispCostEdgeSL_FULLDF,dispCostEdgeSS_FULLDF,
                              dispCostINTSS_EDGEDF,dispCostINTSL_EDGEDF))

colnames(costDataAll)<-c("Dispersal Distance","costSS_Int_full","costSL_Int_full",
                         "costSS_Edge_edge","costSL_Edge_edge","costSL_Edge_full","costSS_Edge_full",
                         "costSS_Int_edge","costSL_Int_edge")

head(costDataAll)

#Unlist all data
costDataAll<-apply(costDataAll,MARGIN = 2,FUN=unlist)

#revert to data frame
costDataAll<-data.frame(costDataAll)

#Calculate differences in RHI between within-patch cost estimation methods

costDataAll$ssDiffInt<-(costDataAll$costSS_Int_edge -costDataAll$costSS_Int_full)/costDataAll$costSS_Int_full*100

costDataAll$slDiffInt<-(costDataAll$costSL_Int_edge -costDataAll$costSL_Int_full)/costDataAll$costSL_Int_full*100

costDataAll$ssDiffEdge<-(costDataAll$costSS_Edge_edge -costDataAll$costSS_Edge_full)/costDataAll$costSS_Edge_full*100

costDataAll$slDiffEdge<-(costDataAll$costSL_Edge_edge -costDataAll$costSL_Edge_full)/costDataAll$costSL_Edge_full*100


###################################Create plots for comparison
##############################################################

intPlot<- ggplot(costDataAll)+geom_line(aes(dispVec,slDiffInt,colour="SL"),linewidth=2)+
  
  geom_line(aes(dispVec,ssDiffInt,colour="SS"),linewidth=2)+
  
  
  xlab("Dispersal distance")+ylab("% Change in RHI")+scale_color_manual(name = "Landscape", 
                                                            
                                                            values = c("SL" = "darkblue", "SS" = "darkred"))#+theme(text = element_text(size = 20))   



intPlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)


edgePlot<- ggplot(costDataAll)+geom_line(aes(dispVec,slDiffEdge,colour="SL"),linewidth=2)+
  
  geom_line(aes(dispVec,ssDiffEdge,colour="SS"),linewidth=2)+
  
  
  xlab("Dispersal distance")+ylab("% Change in RHI")+scale_color_manual(name = "Landscape", 
                                                                        
                                                                        values = c("SL" = "darkblue", "SS" = "darkred"))#+theme(text = element_text(size = 20))   


  

edgePlot+theme_pubr(base_size = 22)+font("legend.text",size=29)



###################################END


