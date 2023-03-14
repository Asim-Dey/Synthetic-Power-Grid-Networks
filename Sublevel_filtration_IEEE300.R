
#install.packages("NetSwan")
library(igraph)
library(NetSwan)

library(EnvStats)
library(TDA)


##############################IEEE 300 Bus system ##############################


Data_IEEE300<-read.csv("IEEE 300 Bus.csv",header=TRUE)
source(Pref_at_G.R')


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


edge_list<-read.table("IEEE300_PoissonNN_10.txt",header=TRUE)


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






#######################################################################################
#################################### CLC ################################################


#edge_list <- read.table("N_IEEE300_10.txt",header=TRUE)

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










#######################################################################################
#################################### GeoDe ################################################



#edge_list <- read.table("edgelistIEEE300_1.txt")



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









################################# ER and PA ###############################################
################################ Simulation ################################################################

source('Pref_at_G.R')

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















