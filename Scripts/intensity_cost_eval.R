
############################################SL####################################################
install.packages(c("dplyr","sf","terra","igraph","stars","ggplot2","ggpubr"))
library(sf)
library(lwgeom)
library(stars)
library(starsExtra)
library(gdistance)
library(dplyr)
library(terra)
library(igraph)
library(stars)
library(terra)

######Load in patch configurations
SL<-st_read("Data/slHomog.shp")
SS<-st_read("Data/ssHomog.shp")
mixed<-st_read("Data/mixedHomog.shp")




#######
#######This function has six arguments - patches = habitat patches: shapefile 
######################################## specialism = one of "interior", "edge" or "generalist"
################################### edge = the size of the edge effect (for a neutral landscape, will have to update this to handle "edge rasters" in real landscapes)
#################################### edgeIntensity = a component (suggest between 0.01 -1) that determines the shape of the kernel (distance decay of edge effect)
#################################### maxDist = maximum dispersal distance of the species being modeled
#####################################dispersalRate = component setting rate of dispersal success
#####################################edgeSensitivity = sets how sensitive the species is to "specialism". Setting this to increasingly small values for either edge or interior moves the species further towards generalist
rhiConnect<-function(patches,specialism,edge,edgeIntensity,maxDist,costRaster,dispersalRate,edgeSensitivity){ 
  
  
  r<-st_rasterize(patches) # need to rasterize the patches to get at edge kernel and area
  r<-rast(r)#convert to terra library, it's better.
  #res(r)
  r<-disagg(r,fact=2)# get data to some resolution to reliably establish relative area later
  
  slClump<-patches(r,directions=4) # get patches (this will account for adjacent cells forming one patch)
  
  
  #####################this section is buffers the whole study area slightly so patches on the edge still get some edge effect
  
  extentNew<-ext(patches) #get extent
  
  
  ext<-as.polygons(extentNew,crs=crs(patches)) #polygonize the extent
  
  #use small buffer distance to simulate edge around patches at limit of study area
  extBuff<-terra::buffer(ext,width=100) 
  
  # extend the data area
  slClump<-extend(slClump,extBuff) 
  
  #set NA value to use when making the the distance raster later on
  
  slClump[is.na(slClump)] <-0 
  #get distance from edge for each patch
  distClump<-gridDistance(slClump$id_lyr.1,target=0)
  #plot(distClump)
  
  # isolate areas within cost threshold 
  rEdge <- distClump <= edge
  
  # get distances within edge threshold
  #rEdge <- rEdge * distClump
  #plot(rEdge)
  #rEdge[rEdge==0] <- NA
  # this sets the edge distance and the shape of neg. exp. curve
  distEdge= log(edgeIntensity)/edge
  # this creates an edge raster where 1 = full edge effect i.e. remove that cell(this happens below) and 0 = no edge effect
  edgeRast<-exp(distEdge*distClump) 
  plot(edgeRast)
  # add zeros back in so we can sum rasters later
  edgeRast[is.na(edgeRast)]<-0
  
  #create inverse for modelling edge species
  invEdge<-1-edgeRast
  #plot(invEdge)
  #Need to select species type. 
  if(specialism=="interior"){
    
  edgeClump<-1-(edgeRast*edgeSensitivity)} #Interior removes all edge cell values from patch cells (i.e. removes a proportion of the area)
  
  else if(specialism=="edge"){edgeClump<- 1-(invEdge*edgeSensitivity)} #For edge just leave as the edge habitat within the patch i.e. inverse of above
  
  else if(specialism=="generalist"){edgeClump=distClump>0}# for generalist all patch cells are one i.e. area not changed
  
  #plot(allRast)
  
  extClump<-extract(edgeClump, patches, fun=sum, method="simple", bind=T)
  #plot(edgeClump)
  head(extClump)
  # get cell area
  cellArea<-res(edgeClump)[1]^2 
  
  extClump$areaMod<-extClump$id_lyr.1*cellArea #create a vector for the new area
  
  
  
  
  distMat<-st_distance(st_as_sf(extClump))# this gets euclidean distances between patches (this is okay for hypothetical landscape as matrix is all assumed to be the same)
  
  #distMat<-st_distance(st_centroid(st_as_sf(extClump)))
  
  distMat<-as.matrix(distMat,nrow=nrow(patches),nrow=nrow(patches)) #create basis for probability matrix
  
  
  
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
  
  AL<-expanse(extBuff)# get study area in m2
  
  
  area<-extClump$areaMod # get area vector for PC metric
  
  #print(area)
  
  areaSum<-sum(area)
  
  
  
  PCmatNew <- outer(area,area)*pstar.mat #get product of all patch areas ij (to the power of something less that 1) and multiply by probabilities above.
  
  pcMatSum<-sum(PCmatNew)
    
    RHI <- pcMatSum/as.numeric(AL^2) #divide by total area of the study squared to get the PC metric  
    #RHI<-sqrt(pcMatSum)/as.numeric(AL)
 
  listMetrics<-list(areaSum,RHI)
  return(listMetrics)
}









###################################################### Intensity evaluation function
#set sequence of edgeIntensity values

intVec<-c(0.01,seq(0.05,0.95,by=0.05),0.99)
intVec
#create function to take values for edgeIntensity and return RHI result
intFun<-function(x,config,speciesGroup){
  
  intRHI<-rhiConnect(patches = config,specialism = speciesGroup,edge = 100,edgeIntensity = x,maxDist = 2600,dispersalRate = 0.05,edgeSensitivity = 1)
  print(x)
  
  return(intRHI)
  
}

#####Run intFun function iterating over intVec values to model changes in RHI
##############################SL

#Edge SL
intEdgeSL<-lapply(intVec,intFun,config=SL,speciesGroup="edge")


intEdgeSL<-do.call("rbind",intEdgeSL)
intEdgeSLDF<-data.frame(intEdgeSL)



#Int SL
intInteriorSL<-lapply(intVec,intFun,config=SL,speciesGroup="interior")


intInteriorSL<-do.call("rbind",intInteriorSL)
intInteriorSLDF<-data.frame(intInteriorSL)

###############SS

#Edge SS
intEdgeSS<-lapply(intVec,intFun,config=SS,speciesGroup="edge")


intEdgeSS<-do.call("rbind",intEdgeSS)
intEdgeSSDF<-data.frame(intEdgeSS)



#Int SS
intInteriorSS<-lapply(intVec,intFun,config=SS,speciesGroup="interior")


intInteriorSS<-do.call("rbind",intInteriorSS)
intInteriorSSDF<-data.frame(intInteriorSS)


############Mixed

#Edge Mixed
intEdgeMixed<-lapply(intVec,intFun,config=mixed,speciesGroup="edge")


intEdgeMixed<-do.call("rbind",intEdgeMixed)
intEdgeMixedDF<-data.frame(intEdgeMixed)



#Int Mixed
intInteriorMixed<-lapply(intVec,intFun,config=mixed,speciesGroup="interior")


intInteriorMixed<-do.call("rbind",intInteriorMixed)
intInteriorMixedDF<-data.frame(intInteriorMixed)








######################################################################







#########
intData<-data.frame(cbind(intVec,intInteriorSLDF[,2],intInteriorSSDF[,2],intInteriorMixed[,2],
                          intEdgeSLDF[,2],intEdgeSSDF[,2],intEdgeMixedDF[,2]))

colnames(intData)<-c("edgeIntensity","interiorSL","interiorSS","interiorMixed",
                         "edgeSL","EdgeSS","edgeMixed")

head(intData)

intData<-apply(intData,MARGIN = 2,FUN=unlist)

intData<-data.frame(intData)


#################################Plotting


library(ggplot2)
library(ggpubr)

slPlot<- ggplot(intData)+geom_line(aes(edgeIntensity,interiorSL,colour="Interior"),linewidth=2)+
  
  geom_line(aes(edgeIntensity,edgeSL,colour="Edge"),linewidth=2)+
  
  
  xlab("Edge Intensity")+ylab("RHI")+scale_color_manual(name = "Species group", 
                                                            
                                                            values = c("Interior" = "darkblue", "Edge" = "darkred"))#+theme(text = element_text(size = 20))   



slPlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)

##################################

head(intData)
ssPlot<- ggplot(intData)+geom_line(aes(edgeIntensity,interiorSS,colour="Interior"),linewidth=2)+
  
  geom_line(aes(edgeIntensity,EdgeSS,colour="Edge"),linewidth=2)+
  
  
  xlab("Edge Intensity")+ylab("RHI")+scale_color_manual(name = "Species group", 
                                                        
                                                        values = c("Interior" = "darkblue", "Edge" = "darkred"))#+theme(text = element_text(size = 20))   



ssPlot+theme_pubr(base_size = 22)+font("legend.text",size=29)#+labs_pubr(base_size = 20)



########################################

