
library(igraph)
library(NetSwan) 
library(EnvStats)


##############################IEEE 300 Bus system ##############################
Data_IEEE300<-read.csv("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Data/IEEE 300 Bus.csv",header=TRUE)
Data_IEEE300<-read.csv("C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Data/IEEE 300 Bus.csv",header=TRUE)

source('C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Pref_at_G.R')
source('C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Pref_at_G.R')

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

G0<-Final_IEEE_300_network


igraph.options(vertex.size=2, vertex.color='black', edge.width=1, edge.color='gray',edge.arrow.size=0.9,
               vertex.label=NA)#vertex.label=NA
plot(G0)


n_node<-length(V(G0))
n_edge<-length(E(G0))


deg_G0<-degree(G0)
summary(deg_G0)
MaxD<-max(deg_G0)
MinD<-min(deg_G0)

meanD_G0<-mean(deg_G0)


#---clustering coefficient, average path length, Diameter--------

AVPL_G<-average.path.length(G0)
AVPL_G 

Diameter_G<-diameter(G0)
Diameter_G 


################################################################

degree1<-degree(G0)
summary(degree1)
table(degree1)

hist(degree1,breaks=10,first.pannel=grid(),xlab="Degree",main="")
box()


meanD1<-mean(degree1);meanD1 
degree_Germany_N<-degree.distribution(G0)[2:12]
degree_Germany<-degree_Germany_N*300 ;degree_Germany #Node 445


Diameter_G
degree_Germany

####################################################################################################

CLC10 <- read.table("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Python/CLC simulated data/N_IEEE300_10.txt",header=TRUE)

CLC10<-as.matrix(CLC10)+1 # few nodes have level 0; changing them t 1

G1<-graph_from_edgelist(CLC10,directed = FALSE)

igraph.options(vertex.size=5, vertex.color='blue', edge.width=1, edge.color='black',
               edge.arrow.size=0.1,
               vertex.label=NA)#vertex.label=NA
plot(G1)


node<-V(G1);n<-length(node);n  # 
edge<-E(G1);m<-length(edge);m  #

#--------------------------------------------------------
degree1<-degree(G1)
summary(degree1)
table(degree1)


degree_g1<-degree.distribution(G1)
hist(degree1,breaks=10,first.pannel=grid(),xlab="Degree",main="")
box()
meanD1<-mean(degree1);meanD1 

degree_g1<-degree.distribution(G1)



################## KL divergence #######################
# X to Y (Q to P)
sim_A0<-degree.distribution(G1)
sim_A<-sim_A0[2:length(sim_A0)];length(sim_A)
####################################################

q1<-sim_A;length(q1)

p0<-degree_Germany_N;length(p0)
n22<-length(q1)-length(p0);n22

p<-c(p0,c(rep(0,n22)));length(p)
q1[q1==0]<-0.0001
p[p==0]<-0.0001

KLD<-sum(p*log(p/q1))
KLD


################## KL divergence #######################
# X to Y (Q to P)
sim_A0<-degree.distribution(G1)
sim_A<-sim_A0[2:length(sim_A0)];length(sim_A)
####################################################
p<-degree_Germany_N;length(p)

q0<-sim_A;length(q0)
#q0<-(sim_10/sum(sim_10))[2:12];length(q0)
n22<-length(p)-length(q0);n22
q1<-c(q0,c(rep(0,n22)));length(q1)

q1[q1==0]<-0.0001
p[p==0]<-0.0001

KLD<-sum(p*log(p/q1))
KLD
   
  
####################################################################################################  
####################################################################################################  

## Two parameters

KL12<-c(0.3208351,0.2991491,0.4453223,0.3919524,0.3756171,0.4431737, 0.3556443,.3256247,0.407112,0.3676975)
mean(KL12) # 0.3732128
sd(KL12)  # 0.04992547


## one parameter############
#KL11<-c(0.337,0.347,0.395,0.336,0.485,0.381,0.364,0.466,0.334,0.359)
#mean(KL11) # 0.3804
#sd(KL11)  # 0.05409087


