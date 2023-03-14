#install.packages("NetSwan")
library(igraph)
library(NetSwan)

##############################Italian Power Grid system ##############################




data11 <- read.csv("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Data/Export_Output22.csv")
data11 <- read.csv("C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Data/Export_Output22.csv")


source('C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Pref_at_G.R')
source('C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Pref_at_G.R')


nodes_Italy<-data11[,c(2,3,4)][data11[,5]=="Italy",]

#write.csv(nodes_Italy,"C:/Users/akd130230/Dropbox/Power Grid Network/Data from author/processed/Italy_nodes.csv")

ID_Italy<-data11[,2][data11[,5]=="Italy"]

#data11[,17] # Fnode
#data11[,18] # Tnode
#Italy<-data11[(data11[,17] %in% ID_Italy) & (data11[,18] %in% ID_Italy)  , ]

Italy<-data11[(data11[,17] %in% ID_Italy) | (data11[,18] %in% ID_Italy)  , ] # or
#write.csv(Italy,"C:/Users/akd130230/Dropbox/Power Grid Network/Data from author/processed/Italy_edge.csv")

Italy_data_all<-Italy[,c(16,17,18)]



Italy_data_all_no_dup=unique(Italy_data_all)

Italy_data<-Italy_data_all_no_dup[,c(2,3)]   # Edge
summary(Italy_data) 

Italy_edge=data.matrix(Italy_data)
Italy_network=graph_from_edgelist(Italy_edge,directed = F)

ed11<-c(as.vector(Italy_edge[,1]),as.vector(Italy_edge[,2]))
unique_ed11<-unique(ed11)
nodes_index=sort(unique_ed11)
m1<-max(nodes_index)
delete_nodes_index=setdiff(c(1:m1),nodes_index)
V(Italy_network)$name=V(Italy_network)

new_Italy_network=delete_vertices(Italy_network,delete_nodes_index)

node<-V(new_Italy_network);node #  273
edge<-E(new_Italy_network);edge # 375
str(new_Italy_network)


G0<-new_Italy_network

igraph.options(vertex.size=2, vertex.color='black', edge.width=1, edge.color='gray',edge.arrow.size=0.9,
               vertex.label=NA)#vertex.label=NA
plot(G0)



n_node<-length(V(G0));n_node
n_edge<-length(E(G0));n_edge



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



degree1<-degree(G0)
summary(degree1)
table(degree1)
#1  2  3  4  5  6  7  8 
#54 84 68 33 21  7  3  3 

hist(degree1,breaks=10,first.pannel=grid(),xlab="Degree",main="")
box()

length(degree1)


meanD1<-mean(degree1);meanD1 



degree_Germany_N<-degree.distribution(G0)[2:9]
degree_Germany<-degree_Germany_N*273 ;degree_Germany #


Diameter_G
degree_Germany

####################################################################################################

CLC10 <- read.table("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Python/Poisson simulated data/Italy_NN_10.txt",header=TRUE)
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



##########################################################################
########################################################################################

KL11<-c(0.08992349,0.1213153,0.1094807,0.11275,0.1127617,0.08318451, 0.1140126,0.07906644,0.1029337,0.1021449)
mean(KL11) # 0.1027573
sd(KL11)   # 0.0142427

























