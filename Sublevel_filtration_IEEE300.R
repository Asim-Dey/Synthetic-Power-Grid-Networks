
#install.packages("NetSwan")
library(igraph)
library(NetSwan)

library(EnvStats)
library(TDA)


##############################IEEE 300 Bus system ##############################


Data_IEEE300<-read.csv("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Data/IEEE 300 Bus.csv",header=TRUE)
#Data_IEEE300<-read.csv("C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/A/Data/IEEE 300 Bus.csv",header=TRUE)

source('C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Pref_at_G.R')
#source('C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/A/Pref_at_G.R')


#-----Remove duplicity -----------------------------------------------------------
Edge_no_dup=unique(Data_IEEE300[,c(1,2)])

#Edge_no_dup=unique(Data_IEEE300[,c(1,2,3)])
#write.table(Edge_no_dup,"IEEE_300_edge.txt", sep = "\t",row.names = FALSE,col.names = TRUE)





IEEE_300_edge=data.matrix(Edge_no_dup)

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




#

####################################################################################################################



########################################
# Examples of Sublevel and Power filtration on Networks
# Original source codes from:
# Cuneyt Gurcan Akcora and Ignacio Segovia-Dominguez
########################################



setwd('C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Data/') 

#setwd('C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Data/') 
options(java.parameters = "-Xmx200g") 


###### Parameters
# Three columns: [From   To   Weigth], we assume these are indirected graphs
nameFileNet = "IEEE_300_edge.txt" 
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
#################################### poisNN ################################################


#edge_list<-read.table("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Python/Poisson simulated data/IEEE300_PoissonNN_10.txt",header=TRUE)

#ED<-as.matrix(edge_list)+1 #  vertex of name sart with 0. So we add 1

#n22<-dim(ED)[1];n22
#W<-rnorm(as.numeric(n22),5,7)
#DD11<-data.frame(ED,W)


#write.table(DD11,"IEEE300_PoissonNN_10.txt", sep = "\t",row.names = FALSE,col.names = TRUE)


###########################################################

# Three columns: [From   To   Weigth], we assume these are indirected graphs
nameFileNet = "IEEE300_PoissonNN_8.txt" 
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

####################################################

# poisNN hub: 37.99335,47.18549 ,44.53237,  55.27763, 48.12679,46.95708, 49.9739, 58.45607 
# poisNN eccentricity:  712, 721, 721, 637, 752, 663,  813.5, 745.5, 669, 827.5




















#######################################################################################
#################################### CLC ################################################


#edge_list <- read.table("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Python/CLC simulated data/N_IEEE300_10.txt",header=TRUE)

ED<-as.matrix(edge_list)+1 #  vertex of name sart with 0. So we add 1

n22<-dim(ED)[1];n22
W<-rnorm(as.numeric(n22),5,7)
DD11<-data.frame(ED,W)


#write.table(DD11,"IEEE_CLC_N_edge_300_10.txt", sep = "\t",row.names = FALSE,col.names = TRUE)


###########################################################

###### Parameters
# Three columns: [From   To   Weigth], we assume these are indirected graphs
nameFileNet = "IEEE_CLC_N_edge_300_10.txt" 
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

####################################################

### Two parameter 
# CLC hub:  29.59726,23.94265,29.31128,33.81334,33.75989,27.10568,38.43247, 28.31472
# CLC eccentricity: 213, 404.5, 223.5, 234.5, 377.5, 345, 297.5, 229,211.5, 376.5



### One parameter 
# CLC hub:  19.30798,9.852413,23.593,17.029,23.260
# CLC eccentricity: 220.532, 299.533, 247.532, 409.032,174.532




# CLC Degree:       
# CLC betweenness: 
# CLC closeness:   

 









#######################################################################################
#################################### GeoDe ################################################



#edge_list <- read.table("C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/Data/edgelistIEEE300_1.txt")
edge_list <- read.table("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Data/edgelistIEEE300_6.txt")



ED<-as.matrix(edge_list)+1 #  vertex of name sart with 0. So we add 1

n22<-dim(ED)[1];n22
W<-rnorm(as.numeric(n22),5,7)
DD11<-data.frame(ED,W)


#write.table(DD11,"IEEE_edge_300_6.txt", sep = "\t",row.names = FALSE,col.names = TRUE)


###########################################################

###### Parameters
# Three columns: [From   To   Weigth], we assume these are indirected graphs
nameFileNet = "IEEE_edge_300_5.txt" 
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
# GeoDe Degree:       301.5, 374,  368, 353.5, 363
# GeoDe betweenness:  12.16582,  11.65484, 12.71054, 11.21288, 12.94173
# GeoDe closeness:   0.945483, 0.5599269,  0.7984239, 0.8669957, 0.9174831
# GeoDe eccentricity: 408.5,  465.5, 68.5, 312, 453.5
# GeoDe hub:          5.230362, 11.06427,  6.43406, 3.967542, 4.240161















################################# ER and PA ###############################################
################################ Simulation ################################################################

source('C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Pref_at_G.R')
source('C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Pref_at_G.R')


n=n_node
m=n_edge


#G1<- sample_gnm(n,m,directed = FALSE, loops = FALSE) #  Erdos-Renyi graphs G(n,m)

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


#write.table(DD44,"IEEE_300_edgetxt_PA_5.txt", sep = "\t",row.names = FALSE,col.names = TRUE)






######################################################################################
#######################################################################################
#######################################  TDA ################################################



###########################################################

###### Parameters
# Three columns: [From   To   Weigth], we assume these are indirected graphs
nameFileNet = "IEEE_300_edgetxt_PA_5.txt" 
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



#ER Degree:         133, 321, 316, 306.5, 228
#ER betweenness:    14.14935, 12.58464, 12.70963, 12.77958, 14.07298
# ER Closeness:     1.661531, 1.216975, 1.86038,  1.249909, 1.330589
# ER eccentricity:  495.5, 470.5, 481.5, 555, 509.5
# ER hub:            33.53718, 34.35302, 36.62647, 34.54533, 29.8837



# PA Degree:       544, 321.5, 517.5, 508, 837.5
# PA betweenness:  7.926248, 5.150593, 10.40006, 8.849849, 5.192771
# PA Closeness:    2.81699, 2.753383, 2.69641, 2.656841, 2.72996
# PA eccentricity:  523.5, 485, 525.5, 507.5, 510.5
# PA hub:           40.59975, 38.77245, 28.31533, 30.65981, 26.75761



# GeoDe Degree:       301.5, 374,  368, 353.5, 363
# GeoDe betweenness:  12.16582,  11.65484, 12.71054, 11.21288, 12.94173
# GeoDe closeness:   0.945483, 0.5599269,  0.7984239, 0.8669957, 0.9174831
# GeoDe eccentricity: 408.5,  465.5, 68.5, 312, 453.5
# GeoDe hub:          5.230362, 11.06427,  6.43406, 3.967542, 4.240161















