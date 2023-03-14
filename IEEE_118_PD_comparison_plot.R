
#install.packages("NetSwan")
library(igraph)
library(NetSwan)

library(EnvStats)
library(TDA)




##############################IEEE 300 Bus system ##############################


Data_IEEE300<-read.csv("IEEE_118_Bus.csv",header=TRUE)


source('Pref_at_G.R')


#-----Remove duplicity -----------------------------------------------------------
Edge_no_dup=unique(Data_IEEE300[,c(1,2)])

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


#####  Shortest path length ##########

AA<-shortest.paths(G0)
summary(as.numeric(AA))




######################################################################################
#######################################################################################
#######################################  TDA ################################################





AA<-as.matrix(AA)
AA[!is.finite(AA)] <- 999

summary(as.numeric(AA[AA<999]))
Mx<-max(as.numeric(AA[AA<999]))

#Mx
#dim(M11/Mx)
#summary(as.numeric(M11/Mx))

A2<- (AA/Mx)
summary(as.numeric(A2))





#----------------------------------------------------------------
#The function ripsDiag computes the persistence diagram of the 
#Rips filtration built on top of a point cloud.

maxdimension <- 1 # max dimension of the homological features to be computed
# 1 for connected components and loops

maxscale <- 0.3    # maximum value of the rips filtration.



DiagRips22 <- ripsDiag(X = A2, maxdimension, maxscale, dist = "arbitrary",
                       library = c("GUDHI", "Dionysus"), location = TRUE, printProgress = FALSE)

plot(DiagRips22[["diagram"]],first.pannel=grid())
#grid()
box()





#------------
D10<-DiagRips22$diagram
grid()
plot(D10, barcode = TRUE, main = "Barcode")
#grid()
#abline(h=118,lwd=1) # Germany 446; Spanish 473


axis(2,at=c(75,  135),tick = FALSE,
     labels=c(expression(paste(H[0])), expression(paste(H[1]))),las=1)

#-------------------------------------------

plot(DiagRips22[["diagram"]],  main = "") #, band = 2 * band[["width"]]




#plot(DiagRips22[["diagram"]], band = 2 * band[["width"]],     main = "KDE Diagram",rotated = TRUE )
#plot(DiagRips22[["diagram"]], rotated = TRUE, band = band[["width"]],   main = "Rotated Diagram")


##########################Landscapes and Silhouettes Persistence landscapes##############################

tseq <- seq(0, maxscale, length = 1000) #domain
Land <- landscape(DiagRips22[["diagram"]], dimension = 1, KK = 1, tseq)
Sil <- silhouette(DiagRips22[["diagram"]], p = 1, dimension = 1, tseq)


plot(tseq, Land, type = "l")
plot(tseq, Sil, type = "l")


############ Distance between Persistence Landscapes###############
# p-landscape distance is proportional to the p-norm of the 
# difference of the two corresponding vectors. 

Land1= landscape(Diag1, dimension=1, KK=1, tseq)
Land2= landscape(Diag2, dimension=1, KK=1, tseq) 


#If you want to compute (an approximation of) the L_p norm: 
L_p=(sum(abs(Land1-Land2)^p)/length(Land1))^(1/p)



################################################################################
######################## GeoDe #######################################################

#################################### GeoDe ################################################

edge_list <- read.table("edgelistIEEE18.txt")
ED<-as.matrix(edge_list)+1 #  vertex of name start with 0. So we add 1

G1<-graph_from_edgelist(ED,directed = FALSE)

igraph.options(vertex.size=2, vertex.color='black', edge.width=1, edge.color='gray',edge.arrow.size=0.9,
               vertex.label=NA)#vertex.label=NA
plot(G1)




n_node1<-length(V(G1));n_node1
n_edge1<-length(E(G1));n_edge1


#####  Shortest path length ##########

AA1<-shortest.paths(G1)
dim(AA1)
summary(as.numeric(AA1))




AA1<-as.matrix(AA1)
AA1[!is.finite(AA1)] <- 999

summary(as.numeric(AA1[AA1<999]))
Mx<-max(as.numeric(AA1[AA1<999]))

#Mx
#dim(M11/Mx)
#summary(as.numeric(M11/Mx))

A2<- (AA1/Mx)
summary(as.numeric(A2))

  

#----------------------------------------------------------------
#The function ripsDiag computes the persistence diagram of the 
#Rips filtration built on top of a point cloud.

maxdimension <- 1 # max dimension of the homological features to be computed
# 1 for connected components and loops

maxscale <- 0.3  # maximum value of the rips filtration.



DiagRips22 <- ripsDiag(X = A2, maxdimension, maxscale, dist = "arbitrary",
                       library = c("GUDHI", "Dionysus"), location = TRUE, printProgress = FALSE)

plot(DiagRips22[["diagram"]],first.pannel=grid())
#grid()
box()





#------------
D10<-DiagRips22$diagram
grid()
plot(D10, barcode = TRUE, main = "Barcode")
#grid()
#abline(h=118,lwd=1) # Germany 446; Spanish 473


axis(2,at=c(75,  135),tick = FALSE,
     labels=c(expression(paste(H[0])), expression(paste(H[1]))),las=1)


plot(DiagRips22[["diagram"]],  main = "") #, band = 2 * band[["width"]]



