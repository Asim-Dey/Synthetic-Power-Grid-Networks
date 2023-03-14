
#install.packages("NetSwan")
library(igraph)
library(NetSwan)

library(EnvStats)
library(TDA)




##############################IEEE 300 Bus system ##############################
Data_IEEE300<-read.csv("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Data/IEEE_118_Bus.csv",header=TRUE)



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

#=============================================================================

#setwd("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/")
setwd("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Current/")



### Check matrix/weight values i.e., summary

cap=2# dimension cap
delta=0.1 # step size for filtration
filt_len=10 # filtration length.
  
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
#################################### PoisNN ################################################


edge_list <- read.table("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Python/Poisson simulated data/IEEE118_PoissonNN_10.txt",header=TRUE)


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



#setwd("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/B")
#setwd("C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/B")


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
delta=0.1 # step size for filtration
filt_len=10 # filtration length.

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


P_GD_1=P


######### Compute WD between PDs ##################################

PD

P_GD_1

WD <- wasserstein(PD, P_GD_1, dimension = c(0,1))
WD


# Diameter: 0.7272727,0.7727273, 0.7727273,0.8636364,1.272727,0.7727273, 0.9090909,0.7727273,0.7727273, 1.318182


































#######################################################################################
#################################### GeoDe ################################################


#edge_list <- read.table("C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/Data/edgelistIEEE18.txt")
edge_list <- read.table("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Data/edgelistIEEE18_9.txt")


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



setwd("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/B")
#setwd("C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/B")


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
delta=0.1 # step size for filtration
filt_len=10 # filtration length.

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


P_GD_1=P


######### Compute WD between PDs ##################################

PD

P_GD_1

WD <- wasserstein(PD, P_GD_1, dimension = c(0,1))
WD




#  1.000,1.545455,2.090909, 1.545455,1.545455,1.545455,2.090909, 1.545455,  1.545455, ,1.545455










################################# ER and PA ###############################################
################################ Simulation ################################################################



n=n_node
m=n_edge
    
    
   # G1<- sample_gnm(n,m,directed = FALSE, loops = FALSE) #  Erdos-Renyi graphs G(n,m)
    
    G1<-PrefAt(n,m) # Preferential Attachment Model with equal number of node and edge 


    
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
    
    
    
    #setwd("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/B")
    setwd("C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/B")
    
    

    
    
    
    
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
    delta=0.1 # step size for filtration
    filt_len=10 # filtration length.
    
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
      P = rbind(P, cbind(rep(2, nrow(persist_data)), persist_data))
      
    } 
    
    
    P_PA3=P
    
  ######### Compute WD between PDs ##################################
    
    PD
    P_PA3
    
    WD <- wasserstein(PD, P_PA3, dimension = c(0,1))
    WD
    

# ER: 13.81818,12.45455, 12.90909,15.63636,4.363636,14.27273,12.90909,3.272727, 15.18182,12.90909

# PA: 10.63636, 10.63636, 1.000,10.63636,10.63636,10.63636,10.63636,1.045455,10.63636,10.63636,





