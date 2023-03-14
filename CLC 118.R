
#install.packages("NetSwan")
library(igraph)
library(NetSwan)
#install.packages("EnvStats")
library(EnvStats)


##############################IEEE 300 Bus system ##############################


Data_IEEE300<-read.csv("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Data/IEEE_118_Bus.csv",header=TRUE)
#Data_IEEE300<-read.csv("C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Data/IEEE_118_Bus.csv",header=TRUE)


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

#-----------------------Statistical Properties --------------------------------------

G0<-Final_IEEE_300_network


igraph.options(vertex.size=5, vertex.color='blue', edge.width=1, edge.color='black',
               edge.arrow.size=0.9,
               vertex.label=NA)#vertex.label=NA
plot(G0)


node<-V(G0);n<-length(node);n  # 118
edge<-E(G0);m<-length(edge);m  #179

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




#------------------
#drawing 3-nodes motifs method ---
par(mfrow=c(1,4))

for(i in 0:3){
  motifgraph <- graph.isocreate(size=3, number=i, directed=F)
  plot(motifgraph)
}

m1<-motifs(G0, 3)
m1[is.na(m1)] <- 0;m1
n02<-count_motifs(G0, 3);n02

T1<-m1[3];T1
T2<-m1[4];T2 # 0 



#--------------

par(mfrow=c(2,6))

for(i in 0:10){
  motifgraph <- graph.isocreate(size=4, number=i, directed=F)
  plot(motifgraph)
}

m2<-motifs(G0, 4)
m2[is.na(m2)] <- 0;m2
n01<-count_motifs(G0, 4);n01




#T2, V3,V4,V4
C0<-c(AVPL_G,Diameter_G)
M0<-c(m1[4],m2[8],m2[9],m2[10])
m2[11]

V1<-m2[5];V1
V2<-m2[7];V2
V3<-m2[8];V3
V4<-m2[9];V4
V5<-m2[10];V5
V6<-m2[11];V6

c(n,m,Diameter_G,AVPL_G,m1[4],m2[8],m2[9],m2[10],m2[11])
# 118.000000 179.000000  14.000000   6.308706  23.000000 132.000000  20.000000   5.000000
#   1.000000

c(T2,n,m)

########################################## Node Degree Distribution  ###########################################

degree1<-degree(G0)
summary(degree1)
table(degree1)

degree_Germany_N<-degree.distribution(G0)[2:10]
degree_Germany<-degree_Germany_N*118 ;degree_Germany #Node 445







####################################################################################################

CLC10 <- read.table("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Python/CLC simulated data/N_IEE118_10.txt",header=TRUE)
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





#####################################################################################################################
######################################################################################################
## Two parameters
KL12<-c(0.7876525,0.8363505,0.7800372, 0.640554,0.6392772,0.8425821, 0.6976287,0.6403641,0.6632167,0.5028068)

mean(KL12) # 0.703047
sd(KL12)  # 0.107465



## one parameter############
#KL11<-c(0.6466418,0.4843424,0.7259464,0.7257489,0.8641921,0.6054892, 0.6845387,0.7117918,0.7548784,0.476377)
#mean(KL11) # 0.6679947
#sd(KL11)  # 0.1200631


  



