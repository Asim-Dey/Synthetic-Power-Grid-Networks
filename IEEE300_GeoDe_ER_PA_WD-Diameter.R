
#install.packages("NetSwan")
library(igraph)
library(NetSwan)

library(EnvStats)
library(TDA)


##############################IEEE 300 Bus system ##############################


Data_IEEE300<-read.csv("IEEE 300 Bus.csv",header=TRUE)
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
 
degree1<-degree(G0)
meanD1<-mean(degree1);meanD1 # 2.726667
sdD1<-sd(degree1)


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

#=============================================================================

### Check matrix/weight values i.e., summary

cap=2# dimension cap
delta=0.05 # step size for filtration
filt_len=20 # filtration length.

d<-n_node

# writing data into file M.txt
cat(d,file='M.txt',append=F,sep = '\n')
cat(paste(0,delta,filt_len,cap,sep = ' '),file = 'M.txt',append = T,sep = '\n') # threshold: 0, 0.1, 0.2,...,1
cat(A2,file='M.txt',append = T) 


#-------------------------------------------------------------------
system('perseusWin.exe distmat M.txt Moutput')

betti_data=as.matrix(read.table('Moutput_betti.txt'))




#write.csv(betti_1, "Betti_1.csv")
#write.csv(betti_0, "Betti_0.csv")

################################################################################
#----------------------WD  ---------------------------

# dim=0
persist_data = as.matrix(read.table('Moutput_0.txt'))
persist_data[persist_data[,2] == -1, 2] = filt_len + 1
persist_data = persist_data/(filt_len + 1)
P = cbind(rep(0, nrow(persist_data)), persist_data)

# dim=1
if (file.info('Moutput_1.txt')$size>0)
{ 
  persist_data = as.matrix(read.table('Moutput_1.txt', blank.lines.skip = T))
  persist_data[persist_data[,2] == -1, 2] = filt_len + 1
  persist_data = persist_data/(filt_len + 1)
  P = rbind(P, cbind(rep(1, nrow(persist_data)), persist_data))
  
}

if (file.info('Moutput_2.txt')$size>0)
{ 
  persist_data = as.matrix(read.table('Moutput_2.txt', blank.lines.skip = T))
  persist_data[persist_data[,2] == -1, 2] = filt_len + 1
  persist_data = persist_data/(filt_len + 1)
  P= rbind(P, cbind(rep(2, nrow(persist_data)), persist_data))
  
}



PD=P




#######################################################################################
#################################### poisNN ################################################

edge_list<-read.table("IEEE300_PoissonNN_4.txt",header=TRUE)



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




######################################################################################
#######################################################################################
#######################################  TDA ################################################



AA1<-as.matrix(AA1)
AA1[!is.finite(AA1)] <- 999

summary(as.numeric(AA1[AA1<999]))
Mx<-max(as.numeric(AA1[AA1<999]))

#Mx
#dim(M11/Mx)
#summary(as.numeric(M11/Mx))

A2<- (AA1/Mx)
summary(as.numeric(A2))



cap=2# dimension cap
delta=0.05 # step size for filtration
filt_len=20 # filtration length.

d<-n_node1

# writing data into file M.txt
cat(d,file='M.txt',append=F,sep = '\n')
cat(paste(0,delta,filt_len,cap,sep = ' '),file = 'M.txt',append = T,sep = '\n') # threshold: 0, 0.1, 0.2,...,1
cat(A2,file='M.txt',append = T) 


#-------------------------------------------------------------------

system('perseusWin.exe distmat M.txt Moutput')


betti_data=as.matrix(read.table('Moutput_betti.txt'))



#write.csv(betti_1, "Betti_1.csv")
#write.csv(betti_0, "Betti_0.csv")

################################################################################
#----------------------WD  ---------------------------

# dim=0
persist_data = as.matrix(read.table('Moutput_0.txt'))
persist_data[persist_data[,2] == -1, 2] = filt_len + 1
persist_data = persist_data/(filt_len + 1)
P_GD = cbind(rep(0, nrow(persist_data)), persist_data)

# dim=1
if (file.info('Moutput_1.txt')$size>0)
{ 
  persist_data = as.matrix(read.table('Moutput_1.txt', blank.lines.skip = T))
  persist_data[persist_data[,2] == -1, 2] = filt_len + 1
  persist_data = persist_data/(filt_len + 1)
  P_GD = rbind(P_GD, cbind(rep(1, nrow(persist_data)), persist_data))
  
}

if (file.info('Moutput_2.txt')$size>0)
{ 
  persist_data = as.matrix(read.table('Moutput_2.txt', blank.lines.skip = T))
  persist_data[persist_data[,2] == -1, 2] = filt_len + 1
  persist_data = persist_data/(filt_len + 1)
  P_GD= rbind(P_GD, cbind(rep(2, nrow(persist_data)), persist_data))
  
}


P_GD


######### Compute WD between PDs ##################################

PD


WD <- wasserstein(PD, P_GD, dimension = c(0,1))
WD



#######################################################################################
#################################### CLC  ################################################

edge_list <- read.table("N_IEEE300_7.txt",header=TRUE)


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




######################################################################################
#######################################################################################
#######################################  TDA ################################################



AA1<-as.matrix(AA1)
AA1[!is.finite(AA1)] <- 999

summary(as.numeric(AA1[AA1<999]))
Mx<-max(as.numeric(AA1[AA1<999]))

#Mx
#dim(M11/Mx)
#summary(as.numeric(M11/Mx))

A2<- (AA1/Mx)
summary(as.numeric(A2))



cap=2# dimension cap
delta=0.05 # step size for filtration
filt_len=20 # filtration length.

d<-n_node1

# writing data into file M.txt
cat(d,file='M.txt',append=F,sep = '\n')
cat(paste(0,delta,filt_len,cap,sep = ' '),file = 'M.txt',append = T,sep = '\n') # threshold: 0, 0.1, 0.2,...,1
cat(A2,file='M.txt',append = T) 


#-------------------------------------------------------------------

system('perseusWin.exe distmat M.txt Moutput')


betti_data=as.matrix(read.table('Moutput_betti.txt'))



#write.csv(betti_1, "Betti_1.csv")
#write.csv(betti_0, "Betti_0.csv")

################################################################################
#----------------------WD  ---------------------------

# dim=0
persist_data = as.matrix(read.table('Moutput_0.txt'))
persist_data[persist_data[,2] == -1, 2] = filt_len + 1
persist_data = persist_data/(filt_len + 1)
P_GD = cbind(rep(0, nrow(persist_data)), persist_data)

# dim=1
if (file.info('Moutput_1.txt')$size>0)
{ 
  persist_data = as.matrix(read.table('Moutput_1.txt', blank.lines.skip = T))
  persist_data[persist_data[,2] == -1, 2] = filt_len + 1
  persist_data = persist_data/(filt_len + 1)
  P_GD = rbind(P_GD, cbind(rep(1, nrow(persist_data)), persist_data))
  
}

if (file.info('Moutput_2.txt')$size>0)
{ 
  persist_data = as.matrix(read.table('Moutput_2.txt', blank.lines.skip = T))
  persist_data[persist_data[,2] == -1, 2] = filt_len + 1
  persist_data = persist_data/(filt_len + 1)
  P_GD= rbind(P_GD, cbind(rep(2, nrow(persist_data)), persist_data))
  
}


P_GD


######### Compute WD between PDs ##################################

PD


WD <- wasserstein(PD, P_GD, dimension = c(0,1))
WD



#######################################################################################
#################################### GeoDe ################################################


edge_list <- read.table("edgelistIEEE300_1.txt")

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




######################################################################################
#######################################################################################
#######################################  TDA ################################################


AA1<-as.matrix(AA1)
AA1[!is.finite(AA1)] <- 999

summary(as.numeric(AA1[AA1<999]))
Mx<-max(as.numeric(AA1[AA1<999]))

#Mx
#dim(M11/Mx)
#summary(as.numeric(M11/Mx))

A2<- (AA1/Mx)
summary(as.numeric(A2))



cap=2# dimension cap
delta=0.05 # step size for filtration
filt_len=20 # filtration length.

d<-n_node1

# writing data into file M.txt
cat(d,file='M.txt',append=F,sep = '\n')
cat(paste(0,delta,filt_len,cap,sep = ' '),file = 'M.txt',append = T,sep = '\n') # threshold: 0, 0.1, 0.2,...,1
cat(A2,file='M.txt',append = T) 


#-------------------------------------------------------------------

system('perseusWin.exe distmat M.txt Moutput')
betti_data=as.matrix(read.table('Moutput_betti.txt'))



#write.csv(betti_1, "Betti_1.csv")
#write.csv(betti_0, "Betti_0.csv")

################################################################################
#----------------------WD  ---------------------------

# dim=0
persist_data = as.matrix(read.table('Moutput_0.txt'))
persist_data[persist_data[,2] == -1, 2] = filt_len + 1
persist_data = persist_data/(filt_len + 1)
P_GD = cbind(rep(0, nrow(persist_data)), persist_data)

# dim=1
if (file.info('Moutput_1.txt')$size>0)
{ 
  persist_data = as.matrix(read.table('Moutput_1.txt', blank.lines.skip = T))
  persist_data[persist_data[,2] == -1, 2] = filt_len + 1
  persist_data = persist_data/(filt_len + 1)
  P_GD = rbind(P_GD, cbind(rep(1, nrow(persist_data)), persist_data))
  
}

if (file.info('Moutput_2.txt')$size>0)
{ 
  persist_data = as.matrix(read.table('Moutput_2.txt', blank.lines.skip = T))
  persist_data[persist_data[,2] == -1, 2] = filt_len + 1
  persist_data = persist_data/(filt_len + 1)
  P_GD= rbind(P_GD, cbind(rep(2, nrow(persist_data)), persist_data))
  
}


P_GD


######### Compute WD between PDs ##################################

PD


WD <- wasserstein(PD, P_GD, dimension = c(0,1))
WD




################################# ER and PA ###############################################
################################ Simulation ################################################################



n=n_node
m=n_edge


G1<- sample_gnm(n,m,directed = FALSE, loops = FALSE) #  Erdos-Renyi graphs G(n,m)

#G1<-PrefAt(n,m) # Preferential Attachment Model with equal number of node and edge 



igraph.options(vertex.size=2, vertex.color='black', edge.width=1, edge.color='gray',edge.arrow.size=0.9,
               vertex.label=NA)#vertex.label=NA
plot(G1)




n_node1<-length(V(G1));n_node1
n_edge1<-length(E(G1));n_edge1


#####  Shortest path length ##########

AA1<-shortest.paths(G1)
dim(AA1)


summary(as.numeric(AA1))



######################################################################################
#######################################################################################
#######################################  TDA ################################################

AA1<-as.matrix(AA1)
AA1[!is.finite(AA1)] <- 999

summary(as.numeric(AA1[AA1<999]))
Mx<-max(as.numeric(AA1[AA1<999]))

#Mx
#dim(M11/Mx)
#summary(as.numeric(M11/Mx))

A2<- (AA1/Mx)
summary(as.numeric(A2))



cap=2# dimension cap
delta=0.05 # step size for filtration
filt_len=20 # filtration length.

d<-n_node1

# writing data into file M.txt
cat(d,file='M.txt',append=F,sep = '\n')
cat(paste(0,delta,filt_len,cap,sep = ' '),file = 'M.txt',append = T,sep = '\n') # threshold: 0, 0.1, 0.2,...,1
cat(A2,file='M.txt',append = T) 


#-------------------------------------------------------------------

system('perseusWin.exe distmat M.txt Moutput')


betti_data=as.matrix(read.table('Moutput_betti.txt'))



#write.csv(betti_1, "Betti_1.csv")
#write.csv(betti_0, "Betti_0.csv")

################################################################################
#----------------------WD  ---------------------------

# dim=0
persist_data = as.matrix(read.table('Moutput_0.txt'))
persist_data[persist_data[,2] == -1, 2] = filt_len + 1
persist_data = persist_data/(filt_len + 1)
P_ER = cbind(rep(0, nrow(persist_data)), persist_data)

# dim=1
if (file.info('Moutput_1.txt')$size>0)
{ 
  persist_data = as.matrix(read.table('Moutput_1.txt', blank.lines.skip = T))
  persist_data[persist_data[,2] == -1, 2] = filt_len + 1
  persist_data = persist_data/(filt_len + 1)
  P_ER = rbind(P_ER, cbind(rep(1, nrow(persist_data)), persist_data))
  
}

if (file.info('Moutput_2.txt')$size>0)
{ 
  persist_data = as.matrix(read.table('Moutput_2.txt', blank.lines.skip = T))
  persist_data[persist_data[,2] == -1, 2] = filt_len + 1
  persist_data = persist_data/(filt_len + 1)
  P_ER = rbind(P_ER, cbind(rep(2, nrow(persist_data)), persist_data))
  
} 

P_ER


######### Compute WD between PDs ##################################



WD <- wasserstein(PD, P_ER, dimension = c(0,1))
WD





