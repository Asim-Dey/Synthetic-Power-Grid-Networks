

library(igraph) 
library(TDA)
#rm(list = ls())



##############################IEEE 118 Bus system ##############################


Data_IEEE300<-read.csv("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Data/IEEE_118_Bus.csv",header=TRUE)
#Data_IEEE300<-read.csv("C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Data/IEEE_118_Bus.csv",header=TRUE)



#-----Remove duplicity -----------------------------------------------------------
Edge_no_dup=unique(Data_IEEE300[,c(1,2)])

IEEE_300_edge=data.matrix(Edge_no_dup)


n33<-dim(IEEE_300_edge)[1];n33
W2<-rnorm(as.numeric(n33),5,7)
DD22<-data.frame(IEEE_300_edge,W2)


#write.table(DD22,"IEEE_edgetxt", sep = "\t",row.names = FALSE,col.names = TRUE)






IEEE_300_network=graph_from_edgelist(IEEE_300_edge,directed = F)

ed11<-c(as.vector(IEEE_300_edge[,1]),as.vector(IEEE_300_edge[,2]))
unique_ed11<-unique(ed11)
nodes_index=sort(unique_ed11)
m1<-max(nodes_index)
delete_nodes_index=setdiff(c(1:m1),nodes_index)
V(IEEE_300_network)$name=V(IEEE_300_network)
#-------------------------------------------------------------------------------------

Final_IEEE_300_network=delete_vertices(IEEE_300_network,delete_nodes_index)

node<-V(Final_IEEE_300_network);node # 300
edge<-E(Final_IEEE_300_network);edge # 409
str(Final_IEEE_300_network)

#igraph.options(vertex.size=2,  edge.arrow.size=0.9,vertex.label=NA)#vertex.label=NA
#plot(Final_IEEE_300_network)


#####################################################################################
#-----------------------Statistical Properties --------------------------------------

G0<-Final_IEEE_300_network


igraph.options(vertex.size=2, vertex.color='black', edge.width=1, edge.color='gray',edge.arrow.size=0.9,
               vertex.label=NA)#vertex.label=NA
plot(G0)


n_node<-length(V(G0));n_node
n_edge<-length(E(G0));n_edge



###################################################################################################################



########################################
# Examples of Sublevel and Power filtration on Networks
# Original source codes from:
# Cuneyt Gurcan Akcora and Ignacio Segovia-Dominguez
########################################



setwd('C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Data/') 

#setwd('C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Data/') 
options(java.parameters = "-Xmx200g") 



#write.table(IEEE_edge,"IEEE_edge.txt", sep = "\t",row.names = FALSE,col.names = TRUE)



###### Parameters
# Three columns: [From   To   Weigth], we assume these are indirected graphs
nameFileNet = "IEEE_edge.txt" 
# Type of filtration 
typeFiltration = 'subnodes' # Sublevel on nodes (subnodes), Power filtration (power)
# max dimension of the homological features to be computed. (e.g. 0 for connected components, 1 for connected components and loops, 2 for connected components, loops, voids, etc.)
maxDimension <- 1

### Next Only for Sublevel Filtration ###
# upper limit on the simplex dimension, size of the cliques to find
maxClique <-2 
# Function to be computed on each node, sublevel filtration is applied on such function
nodeFeature = 'eccentricity' 
# "degree"  "authority"  "closeness"  "betweenness"  "eccentricity"  "hub"

### Notes for Power Filtration ###
# 1) Maximum Scale for filtration is taken as the maximum value on geodesic distances
#    Such value can be modified in variable 'maxScale' in the code below...
# 2) Power filtration works on edges rather than nodes, it means computation of 
#    geodesic distances between pair of nodes. Distance between two adjacent nodes 
#    is based on third column of graph file (i.e. weigths)




###### Loading data  ###########################################


P = read.table(nameFileNet) # Edge List file
if(typeFiltration == 'subnodes'){
  
  
  # Sublevel on Nodes
  P = P[1:2] # Only takes first two columns 
  colnames(P) <- c("source", "target") 
  P$source<-paste0("v", P$source) 
  P$target<-paste0("v", P$target)
}else if(typeFiltration == 'power'){ # Power filtration
  # Takes three columns 
  colnames(P) <- c("source", "target", "weight") 
  P$source<-paste0("v", P$source) 
  P$target<-paste0("v", P$target)
}
edgeData = P 
# Graph object using igraph 



graphNet <- graph.data.frame(edgeData, directed = FALSE) # Build the graph 

###### Computing filtration
if(typeFiltration == 'subnodes'){ ###### Sublevel on Nodes
  # To compute node values...
  # choose the feature values of nodes  
  if (nodeFeature == "degree") {
    nodeValues <- apply(as_adjacency_matrix(graphNet), 1, sum)
  }else if (nodeFeature == "authority") {
    # for undirected graphs, hub scores are equal to authorithy scores
    nodeValues = authority_score(graphNet)$vector
  }else if (nodeFeature == "closeness") {
    nodeValues = closeness(graphNet, normalized=TRUE)
  }else if (nodeFeature == "betweenness") {
    nodeValues = betweenness(graphNet, normalized=TRUE)
  }else if (nodeFeature == "eccentricity") {
    nodeValues = eccentricity(graphNet)
  }else if (nodeFeature == "hub") {
    # for undirected graphs, hub scores are equal to authorithy scores
    nodeValues = hub_score(graphNet)$vector
  }else {
    message(nodeFeature," function not found... ")
  }
  # To clean values in case there are NaN of Inf
  FValsClean = nodeValues 
  FValsClean[is.na(FValsClean)] <- min(nodeValues[is.finite(nodeValues)]) # NaN = minimum value
  FValsClean[!is.finite(FValsClean)] <- max(nodeValues[is.finite(nodeValues)]) # Inf = maximum value
  # for maxClique=3 below means we are adding 0,1,2 simplices (nodes,edges,triangles) to our complex
  
  
  cmplx <- cliques(as.undirected(graphNet), min = 0, max = maxClique)
  # use sublevel=T for sublevel, sublevel=F for superlevel filtration
  # F.values are node values. At these values, node appear in the complex,
  # and their edges to other active nodes also activate ...
  FltRips <- funFiltration(FUNvalues = FValsClean, cmplx = cmplx,  sublevel = T) # Construct filtration using F.values
  #extract the persistence diagram
  # if there is a single activated vertex, the code below will give
  # a warning message.
  
  PD <- filtrationDiag(filtration = FltRips, maxdimension = maxDimension)$diagram
} else if(typeFiltration == 'power'){ ####### Power filtration 
  # Compute matrix of the distances between each node (geodesic distance / short path)
  distMatrixPower_Graph <- shortest.paths(graphNet, v=V(graphNet), to=V(graphNet))
  #diag(distMatrixPower_Graph) <- diagMatDis # Replace using the minimum value
  # To find scale
  valsDMP = unique(c(distMatrixPower_Graph))
  maxScale = max(valsDMP[is.finite(valsDMP)]) # the maximum value on geodesic distances
  # Rips Filtration (it computes complexes)
  FltRips <- ripsFiltration(X = distMatrixPower_Graph, maxdimension = maxDimension,  maxscale = maxScale, dist = "arbitrary", printProgress = FALSE)
  # Persistence Diagram of Power Filtration
  PD <- filtrationDiag(filtration = FltRips, maxdimension = maxDimension)$diagram
}

###### To show persistent Diagrams
print(PD)

#persistenceDiagram[persistenceDiagram[,3]==Inf,3] = max(FValsClean) + 0.01

PD[PD[,3]==Inf,3] = max(FValsClean)

PD


################################################################################################
############################## Simulated #################################################



#######################################################################################
#################################### PoisNN ################################################

#edge_list <- read.table("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Python/Poisson simulated data/IEEE118_PoissonNN_1.txt",header=TRUE)

ED<-as.matrix(edge_list)+1 #  vertex of name start with 0. So we add 1

n22<-dim(ED)[1];n22
W<-rnorm(as.numeric(n22),5,7)
DD11<-data.frame(ED,W)


#write.table(DD11,"IEEE_edge_N_118_Pois_10.txt", sep = "\t",row.names = FALSE,col.names = TRUE)


###########################################################

###### Parameters
# Three columns: [From   To   Weigth], we assume these are indirected graphs
nameFileNet = "IEEE_edge_N_118_Pois_10.txt" 
# Type of filtration 
typeFiltration = 'subnodes' # Sublevel on nodes (subnodes), Power filtration (power)
# max dimension of the homological features to be computed. (e.g. 0 for connected components, 1 for connected components and loops, 2 for connected components, loops, voids, etc.)
maxDimension <- 1

### Next Only for Sublevel Filtration ###
# upper limit on the simplex dimension, size of the cliques to find
maxClique <-2 
# Function to be computed on each node, sublevel filtration is applied on such function
nodeFeature = 'eccentricity' 
# "degree"  "authority"  "closeness"  "betweenness"  "eccentricity"  "hub"

### Notes for Power Filtration ###
# 1) Maximum Scale for filtration is taken as the maximum value on geodesic distances
#    Such value can be modified in variable 'maxScale' in the code below...
# 2) Power filtration works on edges rather than nodes, it means computation of 
#    geodesic distances between pair of nodes. Distance between two adjacent nodes 
#    is based on third column of graph file (i.e. weigths)




###### Loading data  ###########################################


P = read.table(nameFileNet) # Edge List file
if(typeFiltration == 'subnodes'){
  
  
  # Sublevel on Nodes
  P = P[1:2] # Only takes first two columns 
  colnames(P) <- c("source", "target") 
  P$source<-paste0("v", P$source) 
  P$target<-paste0("v", P$target)
}else if(typeFiltration == 'power'){ # Power filtration
  # Takes three columns 
  colnames(P) <- c("source", "target", "weight") 
  P$source<-paste0("v", P$source) 
  P$target<-paste0("v", P$target)
}
edgeData = P 
# Graph object using igraph 



graphNet <- graph.data.frame(edgeData, directed = FALSE) # Build the graph 

###### Computing filtration
if(typeFiltration == 'subnodes'){ ###### Sublevel on Nodes
  # To compute node values...
  # choose the feature values of nodes  
  if (nodeFeature == "degree") {
    nodeValues <- apply(as_adjacency_matrix(graphNet), 1, sum)
  }else if (nodeFeature == "authority") {
    # for undirected graphs, hub scores are equal to authorithy scores
    nodeValues = authority_score(graphNet)$vector
  }else if (nodeFeature == "closeness") {
    nodeValues = closeness(graphNet, normalized=TRUE)
  }else if (nodeFeature == "betweenness") {
    nodeValues = betweenness(graphNet, normalized=TRUE)
  }else if (nodeFeature == "eccentricity") {
    nodeValues = eccentricity(graphNet)
  }else if (nodeFeature == "hub") {
    # for undirected graphs, hub scores are equal to authorithy scores
    nodeValues = hub_score(graphNet)$vector
  }else {
    message(nodeFeature," function not found... ")
  }
  # To clean values in case there are NaN of Inf
  FValsClean = nodeValues 
  FValsClean[is.na(FValsClean)] <- min(nodeValues[is.finite(nodeValues)]) # NaN = minimum value
  FValsClean[!is.finite(FValsClean)] <- max(nodeValues[is.finite(nodeValues)]) # Inf = maximum value
  # for maxClique=3 below means we are adding 0,1,2 simplices (nodes,edges,triangles) to our complex
  cmplx <- cliques(as.undirected(graphNet), min = 0, max = maxClique)
  # use sublevel=T for sublevel, sublevel=F for superlevel filtration
  # F.values are node values. At these values, node appear in the complex,
  # and their edges to other active nodes also activate ...
  FltRips <- funFiltration(FUNvalues = FValsClean, cmplx = cmplx, sublevel = T) # Construct filtration using F.values
  #extract the persistence diagram
  # if there is a single activated vertex, the code below will give
  # a warning message.
  
  PD_Sim <- filtrationDiag(filtration = FltRips, maxdimension = maxDimension)$diagram
  
} else if(typeFiltration == 'power'){ ####### Power filtration 
  # Compute matrix of the distances between each node (geodesic distance / short path)
  distMatrixPower_Graph <- shortest.paths(graphNet, v=V(graphNet), to=V(graphNet))
  #diag(distMatrixPower_Graph) <- diagMatDis # Replace using the minimum value
  # To find scale
  valsDMP = unique(c(distMatrixPower_Graph))
  maxScale = max(valsDMP[is.finite(valsDMP)]) # the maximum value on geodesic distances
  # Rips Filtration (it computes complexes)
  FltRips <- ripsFiltration(X = distMatrixPower_Graph, maxdimension = maxDimension,  maxscale = maxScale, dist = "arbitrary", printProgress = FALSE)
  # Persistence Diagram of Power Filtration
  PD_Sim <- filtrationDiag(filtration = FltRips, maxdimension = maxDimension)$diagram
}

###### To show persistent Diagrams
print(PD_Sim)
 

PD_Sim[PD_Sim[,3]==Inf,3] = max(FValsClean)
PD_Sim


WD <- wasserstein(PD, PD_Sim, dimension = c(0,1))
WD

#######################################################################################

# PoisNN hub: 14.80009, 23.66406,13.9735,12.44698, 14.31908,11.36906,6.007974,19.928, 16.92537,14.11669
# PoisNN eccentricity:146.5,186.5,190.5, 174.5,123.5, 152,71,67.5,206,87.5


























#######################################################################################
#################################### GeoDe ################################################


edge_list <- read.table("C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/Data/edgelistIEEE18.txt")
edge_list <- read.table("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Data/edgelistIEEE18.txt")










ED<-as.matrix(edge_list)+1 #  vertex of name start with 0. So we add 1

n22<-dim(ED)[1];n22
W<-rnorm(as.numeric(n22),5,7)
DD11<-data.frame(ED,W)


#write.table(DD11,"IEEE_edge_118_2.txt", sep = "\t",row.names = FALSE,col.names = TRUE)


###########################################################

###### Parameters
# Three columns: [From   To   Weigth], we assume these are indirected graphs
nameFileNet = "IEEE_edge_118_5.txt" 
# Type of filtration 
typeFiltration = 'subnodes' # Sublevel on nodes (subnodes), Power filtration (power)
# max dimension of the homological features to be computed. (e.g. 0 for connected components, 1 for connected components and loops, 2 for connected components, loops, voids, etc.)
maxDimension <- 1

### Next Only for Sublevel Filtration ###
# upper limit on the simplex dimension, size of the cliques to find
maxClique <-2 
# Function to be computed on each node, sublevel filtration is applied on such function
nodeFeature = 'hub' 
# "degree"  "authority"  "closeness"  "betweenness"  "eccentricity"  "hub"

### Notes for Power Filtration ###
# 1) Maximum Scale for filtration is taken as the maximum value on geodesic distances
#    Such value can be modified in variable 'maxScale' in the code below...
# 2) Power filtration works on edges rather than nodes, it means computation of 
#    geodesic distances between pair of nodes. Distance between two adjacent nodes 
#    is based on third column of graph file (i.e. weigths)




###### Loading data  ###########################################


P = read.table(nameFileNet) # Edge List file
if(typeFiltration == 'subnodes'){
  
  
  # Sublevel on Nodes
  P = P[1:2] # Only takes first two columns 
  colnames(P) <- c("source", "target") 
  P$source<-paste0("v", P$source) 
  P$target<-paste0("v", P$target)
}else if(typeFiltration == 'power'){ # Power filtration
  # Takes three columns 
  colnames(P) <- c("source", "target", "weight") 
  P$source<-paste0("v", P$source) 
  P$target<-paste0("v", P$target)
}
edgeData = P 
# Graph object using igraph 



graphNet <- graph.data.frame(edgeData, directed = FALSE) # Build the graph 

###### Computing filtration
if(typeFiltration == 'subnodes'){ ###### Sublevel on Nodes
  # To compute node values...
  # choose the feature values of nodes  
  if (nodeFeature == "degree") {
    nodeValues <- apply(as_adjacency_matrix(graphNet), 1, sum)
  }else if (nodeFeature == "authority") {
    # for undirected graphs, hub scores are equal to authorithy scores
    nodeValues = authority_score(graphNet)$vector
  }else if (nodeFeature == "closeness") {
    nodeValues = closeness(graphNet, normalized=TRUE)
  }else if (nodeFeature == "betweenness") {
    nodeValues = betweenness(graphNet, normalized=TRUE)
  }else if (nodeFeature == "eccentricity") {
    nodeValues = eccentricity(graphNet)
  }else if (nodeFeature == "hub") {
    # for undirected graphs, hub scores are equal to authorithy scores
    nodeValues = hub_score(graphNet)$vector
  }else {
    message(nodeFeature," function not found... ")
  }
  # To clean values in case there are NaN of Inf
  FValsClean = nodeValues 
  FValsClean[is.na(FValsClean)] <- min(nodeValues[is.finite(nodeValues)]) # NaN = minimum value
  FValsClean[!is.finite(FValsClean)] <- max(nodeValues[is.finite(nodeValues)]) # Inf = maximum value
  # for maxClique=3 below means we are adding 0,1,2 simplices (nodes,edges,triangles) to our complex
  cmplx <- cliques(as.undirected(graphNet), min = 0, max = maxClique)
  # use sublevel=T for sublevel, sublevel=F for superlevel filtration
  # F.values are node values. At these values, node appear in the complex,
  # and their edges to other active nodes also activate ...
  FltRips <- funFiltration(FUNvalues = FValsClean, cmplx = cmplx, sublevel = T) # Construct filtration using F.values
  #extract the persistence diagram
  # if there is a single activated vertex, the code below will give
  # a warning message.
 
  PD_Sim <- filtrationDiag(filtration = FltRips, maxdimension = maxDimension)$diagram
  
} else if(typeFiltration == 'power'){ ####### Power filtration 
  # Compute matrix of the distances between each node (geodesic distance / short path)
  distMatrixPower_Graph <- shortest.paths(graphNet, v=V(graphNet), to=V(graphNet))
  #diag(distMatrixPower_Graph) <- diagMatDis # Replace using the minimum value
  # To find scale
  valsDMP = unique(c(distMatrixPower_Graph))
  maxScale = max(valsDMP[is.finite(valsDMP)]) # the maximum value on geodesic distances
  # Rips Filtration (it computes complexes)
  FltRips <- ripsFiltration(X = distMatrixPower_Graph, maxdimension = maxDimension,  maxscale = maxScale, dist = "arbitrary", printProgress = FALSE)
  # Persistence Diagram of Power Filtration
  PD_Sim <- filtrationDiag(filtration = FltRips, maxdimension = maxDimension)$diagram
}

###### To show persistent Diagrams
print(PD_Sim)


PD_Sim[PD_Sim[,3]==Inf,3] = max(FValsClean)
PD_Sim

 
WD <- wasserstein(PD, PD_Sim, dimension = c(0,1))
WD
   
# 
# GeoDe Degree:        65.94046, 115.6964, 60.94046, 49.94046, 73.19639, 63.19639
# GeoDe betweenness: 5.424816, 2.170455,  5.645916, 6.261743, 4.277916,5.212476
# GeoDe closeness:   1.671583, 1.659508, 1.567951,  1.741564, 1.73356, 1.922228

# GeoDe eccentricity: 72, 33,  60, 68, 71.5, 78.5
# GeoDe hub:  5.084858, 5.923618, 1.99576, 4.638203, 6.879884















################################# ER and PA ###############################################
################################ Simulation ################################################################

source('C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Pref_at_G.R')
source('C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Pref_at_G.R')


n=n_node
m=n_edge


# G1<- sample_gnm(n,m,directed = FALSE, loops = FALSE) #  Erdos-Renyi graphs G(n,m)

G1<-PrefAt(n,m) # Preferential Attachment Model with equal number of node and edge 



igraph.options(vertex.size=2, vertex.color='black', edge.width=1, edge.color='gray',edge.arrow.size=0.9,
               vertex.label=NA)#vertex.label=NA
plot(G1)




n_node1<-length(V(G1));n_node1
n_edge1<-length(E(G1));n_edge1



##################### Edge list ####################################################################

ED3<-as_edgelist(G1)

n44<-dim(ED3)[1];n44
W4<-rnorm(as.numeric(n44),5,7)
DD44<-data.frame(ED3,W4)


#write.table(DD44,"IEEE_edgetxt_PA_5.txt", sep = "\t",row.names = FALSE,col.names = TRUE)






######################################################################################
#######################################################################################
#######################################  TDA ################################################



###########################################################

###### Parameters
# Three columns: [From   To   Weigth], we assume these are indirected graphs
nameFileNet = "IEEE_edgetxt_PA_5.txt" 
# Type of filtration 
typeFiltration = 'subnodes' # Sublevel on nodes (subnodes), Power filtration (power)
# max dimension of the homological features to be computed. (e.g. 0 for connected components, 1 for connected components and loops, 2 for connected components, loops, voids, etc.)
maxDimension <- 1

### Next Only for Sublevel Filtration ###
# upper limit on the simplex dimension, size of the cliques to find
maxClique <-2 
# Function to be computed on each node, sublevel filtration is applied on such function
nodeFeature = 'hub' 
# "degree"  "authority"  "closeness"  "betweenness"  "eccentricity"  "hub"

### Notes for Power Filtration ###
# 1) Maximum Scale for filtration is taken as the maximum value on geodesic distances
#    Such value can be modified in variable 'maxScale' in the code below...
# 2) Power filtration works on edges rather than nodes, it means computation of 
#    geodesic distances between pair of nodes. Distance between two adjacent nodes 
#    is based on third column of graph file (i.e. weigths)




###### Loading data  ###########################################


P = read.table(nameFileNet) # Edge List file
if(typeFiltration == 'subnodes'){
  
  
  # Sublevel on Nodes
  P = P[1:2] # Only takes first two columns 
  colnames(P) <- c("source", "target") 
  P$source<-paste0("v", P$source) 
  P$target<-paste0("v", P$target)
}else if(typeFiltration == 'power'){ # Power filtration
  # Takes three columns 
  colnames(P) <- c("source", "target", "weight") 
  P$source<-paste0("v", P$source) 
  P$target<-paste0("v", P$target)
}
edgeData = P 
# Graph object using igraph 



graphNet <- graph.data.frame(edgeData, directed = FALSE) # Build the graph 

###### Computing filtration
if(typeFiltration == 'subnodes'){ ###### Sublevel on Nodes
  # To compute node values...
  # choose the feature values of nodes  
  if (nodeFeature == "degree") {
    nodeValues <- apply(as_adjacency_matrix(graphNet), 1, sum)
  }else if (nodeFeature == "authority") {
    # for undirected graphs, hub scores are equal to authorithy scores
    nodeValues = authority_score(graphNet)$vector
  }else if (nodeFeature == "closeness") {
    nodeValues = closeness(graphNet, normalized=TRUE)
  }else if (nodeFeature == "betweenness") {
    nodeValues = betweenness(graphNet, normalized=TRUE)
  }else if (nodeFeature == "eccentricity") {
    nodeValues = eccentricity(graphNet)
  }else if (nodeFeature == "hub") {
    # for undirected graphs, hub scores are equal to authorithy scores
    nodeValues = hub_score(graphNet)$vector
  }else {
    message(nodeFeature," function not found... ")
  }
  # To clean values in case there are NaN of Inf
  FValsClean = nodeValues 
  FValsClean[is.na(FValsClean)] <- min(nodeValues[is.finite(nodeValues)]) # NaN = minimum value
  FValsClean[!is.finite(FValsClean)] <- max(nodeValues[is.finite(nodeValues)]) # Inf = maximum value
  # for maxClique=3 below means we are adding 0,1,2 simplices (nodes,edges,triangles) to our complex
  cmplx <- cliques(as.undirected(graphNet), min = 0, max = maxClique)
  # use sublevel=T for sublevel, sublevel=F for superlevel filtration
  # F.values are node values. At these values, node appear in the complex,
  # and their edges to other active nodes also activate ...
  FltRips <- funFiltration(FUNvalues = FValsClean, cmplx = cmplx, sublevel = T) # Construct filtration using F.values
  #extract the persistence diagram
  # if there is a single activated vertex, the code below will give
  # a warning message.
  
  PD_ER <- filtrationDiag(filtration = FltRips, maxdimension = maxDimension)$diagram
  
} else if(typeFiltration == 'power'){ ####### Power filtration 
  # Compute matrix of the distances between each node (geodesic distance / short path)
  distMatrixPower_Graph <- shortest.paths(graphNet, v=V(graphNet), to=V(graphNet))
  #diag(distMatrixPower_Graph) <- diagMatDis # Replace using the minimum value
  # To find scale
  valsDMP = unique(c(distMatrixPower_Graph))
  maxScale = max(valsDMP[is.finite(valsDMP)]) # the maximum value on geodesic distances
  # Rips Filtration (it computes complexes)
  FltRips <- ripsFiltration(X = distMatrixPower_Graph, maxdimension = maxDimension,  maxscale = maxScale, dist = "arbitrary", printProgress = FALSE)
  # Persistence Diagram of Power Filtration
  PD_ER <- filtrationDiag(filtration = FltRips, maxdimension = maxDimension)$diagram
}








###### To show persistent Diagrams
print(PD_ER)

#persistenceDiagram[persistenceDiagram[,3]==Inf,3] = max(FValsClean) + 0.01

PD_ER[PD_ER[,3]==Inf,3] = max(FValsClean)
PD_ER


WD <- wasserstein(PD, PD_ER, dimension = c(0,1))
WD



#ER Degree:         66.94046, 66.19639, 77.94046, 75.19639, 66.94046
#ER betweenness:    6.261504, 6.25239, 6.328218, 6.294665,  6.09119
# ER Closeness:     1.544987, 1.584296, 1.572874, 1.534674, 1.548945
# ER eccentricity: 166, 156.5,  135, 162.5, 163
# ER hub:          14.78512, 16.01811, 20.31846, 16.32957, 20.89478

# PA Degree:       138.4405,  176.9405, 186.4405,  192.9405, 194.4405
# PA betweenness:  5.246792, 3.444533, 1.734693, 3.333176, 3.142678
# PA Closeness:     1.320419, 1.323743, 1.305995, 1.191078,  1.27876
# PA eccentricity: 166.5,  148.5, 139.5, 135.5, 135.5
# PA hub:         14.43425, 16.9797, 14.27741, 15.42017, 16.17929

# GeoDe Degree:        65.94046, 115.6964, 60.94046, 49.94046, 73.19639, 63.19639
# GeoDe betweenness:   5.424816, 2.170455,  5.645916, 6.261743, 4.277916,5.212476
# GeoDe closeness:     1.671583, 1.659508, 1.567951,  1.741564, 1.73356, 1.922228
# GeoDe eccentricity:  72, 33,  60, 68, 71.5, 78.5
# GeoDe hub:  5.084858, 5.923618, 1.99576, 4.638203, 6.879884













