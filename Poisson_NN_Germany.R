
library(igraph)
library(NetSwan)


############### Data #################################################

data11 <- read.csv("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Data/Export_Output22.csv")
data11 <- read.csv("C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Data/Export_Output22.csv")


source('C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Pref_at_G.R')
source('C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Pref_at_G.R')



##############################German Grid##############################


summary(data11) 

#Romania #Poland

nodes_Germany<-data11[,c(2,3,4)][data11[,5]=="Germany",] 
length(nodes_Germany$X)
#write.csv(nodes_Germany,"C:/Users/akd130230/Dropbox/Power Grid Network/Data from author/processed/Germany_nodes.csv")

ID_Germany<-data11[,2][data11[,5]=="Germany"]

#data11[,17] # Fnode
#data11[,18] # Tnode

#Germany<-data11[(data11[,17] %in% ID_Germany) & (data11[,18] %in% ID_Germany)  , ] # Germany edge with distance

Germany<-data11[(data11[,17] %in% ID_Germany) | (data11[,18] %in% ID_Germany)  , ] # or
#write.csv(Germany,"C:/Users/akd130230/Dropbox/Power Grid Network/Data from author/processed/Germany_edge.csv")

Germany_data_all<-Germany[,c(16,17,18)]



Germany_data_all_no_dup=unique(Germany_data_all)

Germany_data<-Germany_data_all_no_dup[,c(2,3)]   # Edge
summary(Germany_data)

Germany_edge=data.matrix(Germany_data)
Germany_network=graph_from_edgelist(Germany_edge,directed = F)

ed11<-c(as.vector(Germany_edge[,1]),as.vector(Germany_edge[,2]))
unique_ed11<-unique(ed11)
nodes_index=sort(unique_ed11)
m1<-max(nodes_index)
delete_nodes_index=setdiff(c(1:m1),nodes_index)
V(Germany_network)$name=V(Germany_network)

new_Germany_network=delete_vertices(Germany_network,delete_nodes_index)

node<-V(new_Germany_network);n=length(node);n #  445
edge<-E(new_Germany_network);m=length(edge);m # 567
str(new_Germany_network)

#plot(new_Germany_network,  vertex.size=5) 

#######################################################################################################

G0<-new_Germany_network

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


degree1<-degree(G0)
summary(degree1)
table(degree1)

hist(degree1,breaks=10,first.pannel=grid(),xlab="Degree",main="")
box()

meanD1<-mean(degree1);meanD1 
degree_Germany_N<-degree.distribution(G0)[2:11]
degree_Germany<-degree_Germany_N*445 ;degree_Germany #Node 445

Diameter_G<-diameter(G0)
Diameter_G 


Diameter_G
degree_Germany


####################################################################################################

CLC10 <- read.table("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Python/Poisson simulated data/Germany_NN_10.txt",header=TRUE)

CLC10<-as.matrix(CLC10)+1 # few nodes have level 0; changing them to 1

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


#############################################################################################

KL12<-c(0.1407386,0.1560455,0.126869,0.1465704,0.1435727,0.1106051,0.140738,0.1030469,0.1557415,0.1630586)
mean(KL12) # 0.1386986
sd(KL12)   # 0.01966117














