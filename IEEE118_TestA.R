
#install.packages("NetSwan")
library(igraph)
library(NetSwan)


library(EnvStats)

##############################IEEE 300 Bus system ##############################


Data_IEEE300<-read.csv("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Data/IEEE_118_Bus.csv",header=TRUE)
Data_IEEE300<-read.csv("C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Data/IEEE_118_Bus.csv",header=TRUE)



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
#-----------------------Statistical Properties --------------------------------------

G0<-Final_IEEE_300_network


igraph.options(vertex.size=2, vertex.color='black', edge.width=1, edge.color='gray',edge.arrow.size=0.9,
               vertex.label=NA)#vertex.label=NA
plot(G0)


#options = {'node_color': 'black',  'node_size': 10, 'width': 0.5,'edge_color': 'gray'}
#nx.draw(G)
#nx.draw_random(G, **options)


#---Degree Distribution--------------



node<-V(G0);n<-length(node);n  # 118
edge<-E(G0);m<-length(edge);m  #179
#str(G1)





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

#--------------------------------------------------------------------

par(mfrow=c(2,6))

for(i in 0:10){
  motifgraph <- graph.isocreate(size=4, number=i, directed=F)
  plot(motifgraph)
}

m2<-motifs(G0, 4)
m2[is.na(m2)] <- 0;m2
n01<-count_motifs(G0, 4);n01


#T2, V3,V4,V4
C0<-c(AVPL_G,Diameter_G);C0
M0<-c(m1[4],m2[8],m2[9],m2[10]);M0




#######################################################################################
#################################### eNN ################################################

G2 <- read.table("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Data/Test_IEEE118_eNN_Motifs.txt", header=FALSE)

names(G2)<-c('nodes','edges','Diameter','AVPL','triangles','square','k13','tent','kite','k4')

head(G2)
dim(G2)


## T2-Triangle; V3/M3-Tent; V4/M4-Square; V5/M5-kite 

#############################################################################################
####################################################################################################


nrep=200

k1<-0;k2<-0;k3<-0;k4<-0
k5<-0;k6<-0;k7<-0


#for(j in 1:nrep){  #j=1
#set.seed(123)

G1<-G2[sample(nrow(G2), 2000), ] # 200
dim(G1)


#----------- Dimaeter------------

Diameter<-G1$Diameter

MO_Dim<-mean(Diameter)               # Mean degree
sO_Dim<-sd(Diameter)                 # SD of Occurance


z_Dim11<-(Diameter_G-MO_Dim)/sO_Dim;    # Z score   

Zrand_Dim<-(Diameter-MO_Dim)/sO_Dim
UQ_Dim11<-quantile(Zrand_Dim, 0.975)        # 1.256084
LQ_Dim11<-quantile(Zrand_Dim, 0.025)        # -2.264109 


#if(abs(z_Dim11)>UQ_Dim11)k1<-k1+1



#----------- AVPL------------

AVPL<-G1$AVPL

MO_ASPL<-mean(AVPL)               # Mean degree
sO_ASPL<-sd(AVPL)                 # SD of Occurance





z_AVPL11<-(AVPL_G-MO_ASPL)/sO_ASPL;    # Z score   


Zrand_AVPL<-(AVPL-MO_ASPL)/sO_ASPL
UQ_AVPL11<-quantile(Zrand_AVPL, 0.975)        # 1.256084
LQ_AVPL11<-quantile(Zrand_AVPL, 0.025)        # -2.264109 

#if(abs(z_AVPL11)>UQ_AVPL11)k2<-k2+1



#--------------------------T2---------------------------------------------------------------------------------

#T2, V3,V4,V5
#c(m1[4],m2[8],m2[9],m2[10])
#T2=m1[4]

T2<-G1$triangles

MO_T2<-mean(T2)               # Mean Occurance
sO_T2<-sd(T2)                 # SD of Occurance



z_T2<-(m1[4]-MO_T2)/sO_T2;    # Z score   


Zrand_T2<-(T2-MO_T2)/sO_T2
UQ_T2<-quantile(Zrand_T2, 0.975)        # 1.256084
LQ_T2<-quantile(Zrand_T2, 0.025)        # -2.264109 

#if(abs(z_T2)>UQ_T2)k3<-k3+1


#--------------------------V3---------------------------------------------------------------------------------

#V3=m2[8]

V3<-G1$tent

MO_V3<-mean(V3)               # Mean Occurance
sO_V3<-sd(V3)                 # SD of Occurance



z_V3<-(m2[8]-MO_V3)/sO_V3;    # Z score   


Zrand_V3<-(V3-MO_V3)/sO_V3
UQ_V3<-quantile(Zrand_V3, 0.975)        # 1.256084
LQ_V3<-quantile(Zrand_V3, 0.025)        # -2.264109 

#if(abs(z_V3)>UQ_V3)k4<-k4+1


#--------------------------V4---------------------------------------------------------------------------------
#V4=m2[9]

V4<-G1$square
MO_V4<-mean(V4)               # Mean Occurance
sO_V4<-sd(V4)                 # SD of Occurance


z_V4<-(m2[9]-MO_V4)/sO_V4;    # Z score   


Zrand_V4<-(V4-MO_V4)/sO_V4
UQ_V4<-quantile(Zrand_V4, 0.975)        # 1.256084
LQ_V4<-quantile(Zrand_V4, 0.025)        # -2.264109 

#if(abs(z_V4)>UQ_V4)k5<-k5+1



#--------------------------V5---------------------------------------------------------------------------------
#V5=m2[10]

V5<-G1$kite

MO_V5<-mean(V5)               # Mean Occurance
sO_V5<-sd(V5)                 # SD of Occurance


z_V5<-(m2[10]-MO_V5)/sO_V5;    # Z score   


Zrand_V5<-(V5-MO_V5)/sO_V5
UQ_V5<-quantile(Zrand_V5, 0.975)        # 1.256084
LQ_V5<-quantile(Zrand_V5, 0.025)        # -2.264109 

#if(abs(z_V5)>UQ_V4)k6<-k6+1

#print(j)





#}






#----------------------------------Table---------------------------------
#--- Occurance
Motif<-factor(c("AVPL", "Diameter","T2","V3","V4","V5"))

Ob_Occur<-round(c(AVPL_G,Diameter_G,m1[4],m2[8],m2[9],m2[10]),3)

Mean<-c(MO_ASPL,MO_Dim,MO_T2,MO_V3,MO_V4,MO_V5)
sd1<-c(sO_ASPL,sO_Dim,sO_T2,sO_V3,sO_V4,sO_V5)

Zm<-round(c(z_AVPL11, z_Dim11,z_T2,z_V3,z_V4,z_V5),3)
LQA<-round(c(LQ_AVPL11,LQ_Dim11,LQ_T2,LQ_V3,LQ_V4,LQ_V5),3)
UQA<-round(c(UQ_AVPL11,UQ_Dim11,UQ_T2,UQ_V3,UQ_V4,UQ_V5),3)


data.frame ('Stat'=Motif, 'Observed'=Ob_Occur, 'Mean'=Mean,'s'=sd1,'Zm'=Zm,'LL'=LQA,'UL'=UQA) 

#     Stat   Observed   Mean        s       Zm     LL    UL
#     AVPL    6.309   5.10878  0.2752473  4.359 -1.703 2.257
# Diameter   14.000  11.08900  1.1509741  2.529 -1.815 2.529
#       T2   23.000  88.23450 13.4695666 -4.843 -1.799 2.136
#       V3  132.000 445.52850 88.8356070 -3.529 -1.661 2.223
#       V4   20.000  20.86500  6.8486661 -0.126 -1.586 2.210
#       V5    5.000  78.92500 20.6111325 -3.587 -1.646 2.187





#######################################################################################
#################################### CLC 1 param ################################################

G2 <- read.table("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Data/Test_IEEE118_1parar_Motifs.txt", header=FALSE)[1:2000,]

names(G2)<-c('nodes','edges','Diameter','AVPL','triangles','square','k13','tent','kite','k4')

head(G2)
dim(G2)


## T2-Triangle; V3/M3-Tent; V4/M4-Square; V5/M5-kite 

#############################################################################################
####################################################################################################


nrep=200

k1<-0;k2<-0;k3<-0;k4<-0
k5<-0;k6<-0;k7<-0


#for(j in 1:nrep){  #j=1
#set.seed(123)

G1<-G2[sample(nrow(G2), 2000), ] # 200
dim(G1)


#----------- Dimaeter------------

Diameter<-G1$Diameter

MO_Dim<-mean(Diameter)               # Mean degree
sO_Dim<-sd(Diameter)                 # SD of Occurance


z_Dim11<-(Diameter_G-MO_Dim)/sO_Dim;    # Z score   

Zrand_Dim<-(Diameter-MO_Dim)/sO_Dim
UQ_Dim11<-quantile(Zrand_Dim, 0.975)        # 1.256084
LQ_Dim11<-quantile(Zrand_Dim, 0.025)        # -2.264109 


#if(abs(z_Dim11)>UQ_Dim11)k1<-k1+1



#----------- AVPL------------

AVPL<-G1$AVPL

MO_ASPL<-mean(AVPL)               # Mean degree
sO_ASPL<-sd(AVPL)                 # SD of Occurance





z_AVPL11<-(AVPL_G-MO_ASPL)/sO_ASPL;    # Z score   


Zrand_AVPL<-(AVPL-MO_ASPL)/sO_ASPL
UQ_AVPL11<-quantile(Zrand_AVPL, 0.975)        # 1.256084
LQ_AVPL11<-quantile(Zrand_AVPL, 0.025)        # -2.264109 

#if(abs(z_AVPL11)>UQ_AVPL11)k2<-k2+1



#--------------------------T2---------------------------------------------------------------------------------

#T2, V3,V4,V5
#c(m1[4],m2[8],m2[9],m2[10])
#T2=m1[4]

T2<-G1$triangles

MO_T2<-mean(T2)               # Mean Occurance
sO_T2<-sd(T2)                 # SD of Occurance



z_T2<-(m1[4]-MO_T2)/sO_T2;    # Z score   


Zrand_T2<-(T2-MO_T2)/sO_T2
UQ_T2<-quantile(Zrand_T2, 0.975)        # 1.256084
LQ_T2<-quantile(Zrand_T2, 0.025)        # -2.264109 

#if(abs(z_T2)>UQ_T2)k3<-k3+1


#--------------------------V3---------------------------------------------------------------------------------

#V3=m2[8]

V3<-G1$tent

MO_V3<-mean(V3)               # Mean Occurance
sO_V3<-sd(V3)                 # SD of Occurance



z_V3<-(m2[8]-MO_V3)/sO_V3;    # Z score   


Zrand_V3<-(V3-MO_V3)/sO_V3
UQ_V3<-quantile(Zrand_V3, 0.975)        # 1.256084
LQ_V3<-quantile(Zrand_V3, 0.025)        # -2.264109 

#if(abs(z_V3)>UQ_V3)k4<-k4+1


#--------------------------V4---------------------------------------------------------------------------------
#V4=m2[9]

V4<-G1$square
MO_V4<-mean(V4)               # Mean Occurance
sO_V4<-sd(V4)                 # SD of Occurance


z_V4<-(m2[9]-MO_V4)/sO_V4;    # Z score   


Zrand_V4<-(V4-MO_V4)/sO_V4
UQ_V4<-quantile(Zrand_V4, 0.975)        # 1.256084
LQ_V4<-quantile(Zrand_V4, 0.025)        # -2.264109 

#if(abs(z_V4)>UQ_V4)k5<-k5+1



#--------------------------V5---------------------------------------------------------------------------------
#V5=m2[10]

V5<-G1$kite

MO_V5<-mean(V5)               # Mean Occurance
sO_V5<-sd(V5)                 # SD of Occurance


z_V5<-(m2[10]-MO_V5)/sO_V5;    # Z score   


Zrand_V5<-(V5-MO_V5)/sO_V5
UQ_V5<-quantile(Zrand_V5, 0.975)        # 1.256084
LQ_V5<-quantile(Zrand_V5, 0.025)        # -2.264109 

#if(abs(z_V5)>UQ_V4)k6<-k6+1

#print(j)





#}






#----------------------------------Table---------------------------------
#--- Occurance
Motif<-factor(c("AVPL", "Diameter","T2","V3","V4","V5"))

Ob_Occur<-round(c(AVPL_G,Diameter_G,m1[4],m2[8],m2[9],m2[10]),3)

Mean<-c(MO_ASPL,MO_Dim,MO_T2,MO_V3,MO_V4,MO_V5)
sd1<-c(sO_ASPL,sO_Dim,sO_T2,sO_V3,sO_V4,sO_V5)

Zm<-round(c(z_AVPL11, z_Dim11,z_T2,z_V3,z_V4,z_V5),3)
LQA<-round(c(LQ_AVPL11,LQ_Dim11,LQ_T2,LQ_V3,LQ_V4,LQ_V5),3)
UQA<-round(c(UQ_AVPL11,UQ_Dim11,UQ_T2,UQ_V3,UQ_V4,UQ_V5),3)


data.frame ('Stat'=Motif, 'Observed'=Ob_Occur, 'Mean'=Mean,'s'=sd1,'Zm'=Zm,'LL'=LQA,'UL'=UQA) 


#     Stat   Observed Mean          s        Zm     LL    UL
#     AVPL    6.309  7.347615  0.5517621 -1.883 -1.409 1.472
# Diameter   14.000 18.383500  2.5234127 -1.737 -1.737 1.829
#       T2   23.000 11.017500  3.8839073  3.085 -1.807 1.798
#       V3  132.000 47.909500 17.6484289  4.765 -1.695 2.215
#       V4   20.000  9.969000  3.2025251  3.132 -1.552 2.203
#       V5    5.000  2.563000  1.9819170  1.230 -1.293 2.239














#######################################################################################
#################################### CLC 2 param ################################################

G2 <- read.table("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Data/Test_IEEE118_2parar_Motifs.txt", header=FALSE)
names(G2)<-c('nodes','edges','Diameter','AVPL','triangles','square','k13','tent','kite','k4')

head(G2)
dim(G2)


## T2-Triangle; V3/M3-Tent; V4/M4-Square; V5/M5-kite 

#############################################################################################
####################################################################################################


nrep=200

k1<-0;k2<-0;k3<-0;k4<-0
k5<-0;k6<-0;k7<-0


#for(j in 1:nrep){  #j=1
#set.seed(123)

G1<-G2[sample(nrow(G2), 2000), ] # 200
dim(G1)


#----------- Dimaeter------------

Diameter<-G1$Diameter

MO_Dim<-mean(Diameter)               # Mean degree
sO_Dim<-sd(Diameter)                 # SD of Occurance


z_Dim11<-(Diameter_G-MO_Dim)/sO_Dim;    # Z score   

Zrand_Dim<-(Diameter-MO_Dim)/sO_Dim
UQ_Dim11<-quantile(Zrand_Dim, 0.975)        # 1.256084
LQ_Dim11<-quantile(Zrand_Dim, 0.025)        # -2.264109 


#if(abs(z_Dim11)>UQ_Dim11)k1<-k1+1



#----------- AVPL------------

AVPL<-G1$AVPL

MO_ASPL<-mean(AVPL)               # Mean degree
sO_ASPL<-sd(AVPL)                 # SD of Occurance





z_AVPL11<-(AVPL_G-MO_ASPL)/sO_ASPL;    # Z score   


Zrand_AVPL<-(AVPL-MO_ASPL)/sO_ASPL
UQ_AVPL11<-quantile(Zrand_AVPL, 0.975)        # 1.256084
LQ_AVPL11<-quantile(Zrand_AVPL, 0.025)        # -2.264109 

#if(abs(z_AVPL11)>UQ_AVPL11)k2<-k2+1



#--------------------------T2---------------------------------------------------------------------------------

#T2, V3,V4,V5
#c(m1[4],m2[8],m2[9],m2[10])
#T2=m1[4]

T2<-G1$triangles

MO_T2<-mean(T2)               # Mean Occurance
sO_T2<-sd(T2)                 # SD of Occurance



z_T2<-(m1[4]-MO_T2)/sO_T2;    # Z score   


Zrand_T2<-(T2-MO_T2)/sO_T2
UQ_T2<-quantile(Zrand_T2, 0.975)        # 1.256084
LQ_T2<-quantile(Zrand_T2, 0.025)        # -2.264109 

#if(abs(z_T2)>UQ_T2)k3<-k3+1


#--------------------------V3---------------------------------------------------------------------------------

#V3=m2[8]

V3<-G1$tent

MO_V3<-mean(V3)               # Mean Occurance
sO_V3<-sd(V3)                 # SD of Occurance



z_V3<-(m2[8]-MO_V3)/sO_V3;    # Z score   


Zrand_V3<-(V3-MO_V3)/sO_V3
UQ_V3<-quantile(Zrand_V3, 0.975)        # 1.256084
LQ_V3<-quantile(Zrand_V3, 0.025)        # -2.264109 

#if(abs(z_V3)>UQ_V3)k4<-k4+1


#--------------------------V4---------------------------------------------------------------------------------
#V4=m2[9]

V4<-G1$square
MO_V4<-mean(V4)               # Mean Occurance
sO_V4<-sd(V4)                 # SD of Occurance


z_V4<-(m2[9]-MO_V4)/sO_V4;    # Z score   


Zrand_V4<-(V4-MO_V4)/sO_V4
UQ_V4<-quantile(Zrand_V4, 0.975)        # 1.256084
LQ_V4<-quantile(Zrand_V4, 0.025)        # -2.264109 

#if(abs(z_V4)>UQ_V4)k5<-k5+1



#--------------------------V5---------------------------------------------------------------------------------
#V5=m2[10]

V5<-G1$kite

MO_V5<-mean(V5)               # Mean Occurance
sO_V5<-sd(V5)                 # SD of Occurance


z_V5<-(m2[10]-MO_V5)/sO_V5;    # Z score   


Zrand_V5<-(V5-MO_V5)/sO_V5
UQ_V5<-quantile(Zrand_V5, 0.975)        # 1.256084
LQ_V5<-quantile(Zrand_V5, 0.025)        # -2.264109 

#if(abs(z_V5)>UQ_V4)k6<-k6+1

#print(j)





#}






#----------------------------------Table---------------------------------
#--- Occurance
Motif<-factor(c("AVPL", "Diameter","T2","V3","V4","V5"))

Ob_Occur<-round(c(AVPL_G,Diameter_G,m1[4],m2[8],m2[9],m2[10]),3)

Mean<-c(MO_ASPL,MO_Dim,MO_T2,MO_V3,MO_V4,MO_V5)
sd1<-c(sO_ASPL,sO_Dim,sO_T2,sO_V3,sO_V4,sO_V5)

Zm<-round(c(z_AVPL11, z_Dim11,z_T2,z_V3,z_V4,z_V5),3)
LQA<-round(c(LQ_AVPL11,LQ_Dim11,LQ_T2,LQ_V3,LQ_V4,LQ_V5),3)
UQA<-round(c(UQ_AVPL11,UQ_Dim11,UQ_T2,UQ_V3,UQ_V4,UQ_V5),3)


data.frame ('Stat'=Motif, 'Observed'=Ob_Occur, 'Mean'=Mean,'s'=sd1,'Zm'=Zm,'LL'=LQA,'UL'=UQA) 

#      Stat Observed    Mean      s         Zm     LL    UL
#     AVPL    6.309  5.961035  0.3419486  1.017 -1.816 2.044
# Diameter   14.000 14.446000  1.6505897 -0.270 -1.482 2.153
#       T2   23.000  9.107000  2.9971724  4.635 -1.704 1.966
#       V3  132.000 42.897000 14.4159738  6.181 -1.796 2.019
#       V4   20.000  6.078500  2.7533849  5.056 -1.481 2.151
#       V5    5.000  2.646000  2.0767379  1.134 -1.274 2.578












#######################################################################################
#################################### GeoDe ################################################

GG1 <- read.csv("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Data/IEEE118_6000.txt", header=FALSE)
GG2 <- read.csv("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Data/IEEE118_6000B.txt", header=FALSE)
#GG3 <- read.csv("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Data/IEEE118_6000C.txt", header=FALSE)

G2<-rbind(GG1,GG2)

dim(G2) 

names(G2)<-c("m", "t","nodes", "edges", "vee", "triangles", "path3", "square", 
             "k13",  "tent", "kite", "k4","Diameter", "AVPL")

head(G2)
dim(G2)

## T2-Triangle; V3/M3-Tent; V4/M4-Square; V5/M5-kite 

#############################################################################################
####################################################################################################


nrep=200

k1<-0;k2<-0;k3<-0;k4<-0
k5<-0;k6<-0;k7<-0


#for(j in 1:nrep){  #j=1
  #set.seed(123)
  
  G1<-G2[sample(nrow(G2), 2000), ] # 200
  dim(G1)
  
  
  #----------- Dimaeter------------
  
  Diameter<-G1$Diameter
  
  MO_Dim<-mean(Diameter)               # Mean degree
  sO_Dim<-sd(Diameter)                 # SD of Occurance
  

  z_Dim11<-(Diameter_G-MO_Dim)/sO_Dim;    # Z score   
  
  Zrand_Dim<-(Diameter-MO_Dim)/sO_Dim
  UQ_Dim11<-quantile(Zrand_Dim, 0.975)        # 1.256084
  LQ_Dim11<-quantile(Zrand_Dim, 0.025)        # -2.264109 
  
  
  #if(abs(z_Dim11)>UQ_Dim11)k1<-k1+1
  
  
  
  #----------- AVPL------------
  
  AVPL<-G1$AVPL
  
  MO_ASPL<-mean(AVPL)               # Mean degree
  sO_ASPL<-sd(AVPL)                 # SD of Occurance
  
  

  
  
  z_AVPL11<-(AVPL_G-MO_ASPL)/sO_ASPL;    # Z score   
  
  
  Zrand_AVPL<-(AVPL-MO_ASPL)/sO_ASPL
  UQ_AVPL11<-quantile(Zrand_AVPL, 0.975)        # 1.256084
  LQ_AVPL11<-quantile(Zrand_AVPL, 0.025)        # -2.264109 
  
  #if(abs(z_AVPL11)>UQ_AVPL11)k2<-k2+1
  
 
  
  #--------------------------T2---------------------------------------------------------------------------------
  
  #T2, V3,V4,V5
  #c(m1[4],m2[8],m2[9],m2[10])
  #T2=m1[4]
  
  T2<-G1$triangles
  
  MO_T2<-mean(T2)               # Mean Occurance
  sO_T2<-sd(T2)                 # SD of Occurance
  
  
  
  z_T2<-(m1[4]-MO_T2)/sO_T2;    # Z score   
  
  
  Zrand_T2<-(T2-MO_T2)/sO_T2
  UQ_T2<-quantile(Zrand_T2, 0.975)        # 1.256084
  LQ_T2<-quantile(Zrand_T2, 0.025)        # -2.264109 
  
  #if(abs(z_T2)>UQ_T2)k3<-k3+1
  
 
  #--------------------------V3---------------------------------------------------------------------------------
  
  #V3=m2[8]
  
  V3<-G1$tent
  
  MO_V3<-mean(V3)               # Mean Occurance
  sO_V3<-sd(V3)                 # SD of Occurance
  
  
  
  z_V3<-(m2[8]-MO_V3)/sO_V3;    # Z score   
  
  
  Zrand_V3<-(V3-MO_V3)/sO_V3
  UQ_V3<-quantile(Zrand_V3, 0.975)        # 1.256084
  LQ_V3<-quantile(Zrand_V3, 0.025)        # -2.264109 
  
  #if(abs(z_V3)>UQ_V3)k4<-k4+1
  
  
  #--------------------------V4---------------------------------------------------------------------------------
  #V4=m2[9]
  
  V4<-G1$square
  MO_V4<-mean(V4)               # Mean Occurance
  sO_V4<-sd(V4)                 # SD of Occurance
  
  
   z_V4<-(m2[9]-MO_V4)/sO_V4;    # Z score   
  
  
  Zrand_V4<-(V4-MO_V4)/sO_V4
  UQ_V4<-quantile(Zrand_V4, 0.975)        # 1.256084
  LQ_V4<-quantile(Zrand_V4, 0.025)        # -2.264109 
  
  #if(abs(z_V4)>UQ_V4)k5<-k5+1
  
  
  
  #--------------------------V5---------------------------------------------------------------------------------
  #V5=m2[10]
  
  V5<-G1$kite
  
  MO_V5<-mean(V5)               # Mean Occurance
  sO_V5<-sd(V5)                 # SD of Occurance
  
  
  z_V5<-(m2[10]-MO_V5)/sO_V5;    # Z score   
  
  
  Zrand_V5<-(V5-MO_V5)/sO_V5
  UQ_V5<-quantile(Zrand_V5, 0.975)        # 1.256084
  LQ_V5<-quantile(Zrand_V5, 0.025)        # -2.264109 
  
  #if(abs(z_V5)>UQ_V4)k6<-k6+1
  
  #print(j)
  
 
  
  
  
#}






#----------------------------------Table---------------------------------
#--- Occurance
Motif<-factor(c("AVPL", "Diameter","T2","V3","V4","V5"))

Ob_Occur<-round(c(AVPL_G,Diameter_G,m1[4],m2[8],m2[9],m2[10]),3)

Mean<-c(MO_ASPL,MO_Dim,MO_T2,MO_V3,MO_V4,MO_V5)
sd1<-c(sO_ASPL,sO_Dim,sO_T2,sO_V3,sO_V4,sO_V5)

Zm<-round(c(z_AVPL11, z_Dim11,z_T2,z_V3,z_V4,z_V5),3)
LQA<-round(c(LQ_AVPL11,LQ_Dim11,LQ_T2,LQ_V3,LQ_V4,LQ_V5),3)
UQA<-round(c(UQ_AVPL11,UQ_Dim11,UQ_T2,UQ_V3,UQ_V4,UQ_V5),3)


data.frame ('Stat'=Motif, 'Observed'=Ob_Occur, 'Mean'=Mean,'s'=sd1,'Zm'=Zm,'LL'=LQA,'UL'=UQA) 


#      Stat Observed    Mean         s     Zm     LL    UL
#     AVPL    6.309   5.943686  0.525153  0.695 -1.476 2.248
# Diameter   14.000  13.723000  2.044593  0.135 -1.332 2.092
#       T2   23.000  23.445000  5.939917 -0.075 -1.758 1.950
#       V3  132.000 121.298000 35.019171  0.306 -1.693 2.219
#       V4   20.000  15.821500  5.386478  0.776 -1.638 2.075
#       V5    5.000   9.119500  4.915415 -0.838 -1.448 2.417




################################# ER and PA ###############################################
################################ Simulation ################################################################


#nrep=200

#k1<-0;k2<-0;k3<-0;k4<-0
#k5<-0;k6<-0;k7<-0


#for(j in 1:nrep){  #j=1
  #set.seed(123)
  
  n=n # node
  m=m # edge
  
  B=2000
  
  
  M<-numeric(B)
  AVPL<-numeric(B)
  Diameter<-numeric(B)
  
  V1<-numeric(B);V2<-numeric(B);V3<-numeric(B);V4<-numeric(B)
  V5<-numeric(B);V6<-numeric(B);T1<-numeric(B);T2<-numeric(B)
  
  
  C_V1<-numeric(B);C_V2<-numeric(B);C_V3<-numeric(B);C_V4<-numeric(B)
  C_V5<-numeric(B);C_V6<-numeric(B);C_T1<-numeric(B);C_T2<-numeric(B)
  
  for(i in 1:B) {  #i=1
    
 
    #g1<- sample_gnm(n,m,directed = FALSE, loops = FALSE) #  Erdos-Renyi graphs G(n,m)
    
    g1<-PrefAt(n,m) # Preferential Attachment Model with equal number of node and edge 
    
    #igraph.options(vertex.size=2,  edge.arrow.size=0.9,vertex.label=NA)#vertex.label=NA
    #plot(g1,  vertex.size=5,vertex.color="red")
    
    
    ################################### Statistics ################################# 
    
    D1<-degree(g1)
    M[i]<-mean(D1) 
    
    
    
    AVPL[i]<-average.path.length(g1)
    
    Diameter[i]<-diameter(g1)
    
    #-----------------------------------------------------------------------------------------------  
    
    mg3<-motifs(g1, 4)
    mg3[is.na(mg3)] <- 0
    n1<-count_motifs(g1, 4)
    
    #------------------------------------------------------------------
    
    V1[i]<-mg3[5];C_V1[i]<-mg3[5]/n1
    
    V2[i]<-mg3[7];C_V2[i]<-mg3[7]/n1
    
    V3[i]<-mg3[8];C_V3[i]<-mg3[8]/n1
    V4[i]<-mg3[9];C_V4[i]<-mg3[9]/n1
    V5[i]<-mg3[10];C_V5[i]<-mg3[10]/n1
    V6[i]<-mg3[11];C_V6[i]<-mg3[11]/n1
    
    #--------------------------------------
    
    mg1<-motifs(g1, 3)
    mg1[is.na(mg1)] <- 0
    n2<-count_motifs(g1, 3)
    
    T1[i]<-mg1[3];C_T1[i]<-mg1[3]/n2
    T2[i]<-mg1[4];C_T2[i]<-mg1[4]/n2
    
    
    #----------------------------------------------------------------------
    # print(i)
    
  }
  
  
  
  #----------- Dimaeter------------------------------------------------
  
  
  MO_Dim<-mean(Diameter)               # Mean degree
  sO_Dim<-sd(Diameter)                 # SD of Occurance
  
  
  z_Dim11<-(Diameter_G-MO_Dim)/sO_Dim;    # Z score   
  
  Zrand_Dim<-(Diameter-MO_Dim)/sO_Dim
  UQ_Dim11<-quantile(Zrand_Dim, 0.975)        # 1.256084
  LQ_Dim11<-quantile(Zrand_Dim, 0.025)        # -2.264109 
  
  
  #if(abs(z_Dim11)>UQ_Dim11)k1<-k1+1
  
  
  
  
  
  
  
  #----------- AVPL----------------------------------------------------
  MO_ASPL<-mean(AVPL)               # Mean degree
  sO_ASPL<-sd(AVPL)                 # SD of Occurance
  
  z_AVPL11<-(AVPL_G-MO_ASPL)/sO_ASPL;    # Z score   
  
  
  Zrand_AVPL<-(AVPL-MO_ASPL)/sO_ASPL
  UQ_AVPL11<-quantile(Zrand_AVPL, 0.975)        # 1.256084
  LQ_AVPL11<-quantile(Zrand_AVPL, 0.025)        # -2.264109 
  
  #if(abs(z_AVPL11)>UQ_AVPL11)k2<-k2+1
  
  
  
  #--------------------------T2--------------------------------------------
  #T2, V3,V4,V5
  #c(m1[4],m2[8],m2[9],m2[10])
  #T2=m1[4]
  
  
  MO_T2<-mean(T2)               # Mean Occurance
  sO_T2<-sd(T2)                 # SD of Occurance
  
  z_T2<-(m1[4]-MO_T2)/sO_T2;    # Z score   
  
  
  Zrand_T2<-(T2-MO_T2)/sO_T2
  UQ_T2<-quantile(Zrand_T2, 0.975)        # 1.256084
  LQ_T2<-quantile(Zrand_T2, 0.025)        # -2.264109 
  
  #if(abs(z_T2)>UQ_T2)k3<-k3+1
  
  
  #--------------------------V3----------------------------------------------
  #V3=m2[8]
  
  MO_V3<-mean(V3)               # Mean Occurance
  sO_V3<-sd(V3)                 # SD of Occurance
  
  z_V3<-(m2[8]-MO_V3)/sO_V3;    # Z score   
  
  
  Zrand_V3<-(V3-MO_V3)/sO_V3
  UQ_V3<-quantile(Zrand_V3, 0.975)        # 1.256084
  LQ_V3<-quantile(Zrand_V3, 0.025)        # -2.264109 
  
  #if(abs(z_V3)>UQ_V3)k4<-k4+1
  
  
  #--------------------------V4----------------------------------------------
  #V4=m2[9]
  
  
  MO_V4<-mean(V4)               # Mean Occurance
  sO_V4<-sd(V4)                 # SD of Occurance
  
  z_V4<-(m2[9]-MO_V4)/sO_V4;    # Z score   
  
  
  Zrand_V4<-(V4-MO_V4)/sO_V4
  UQ_V4<-quantile(Zrand_V4, 0.975)        # 1.256084
  LQ_V4<-quantile(Zrand_V4, 0.025)        # -2.264109 
  
  #if(abs(z_V4)>UQ_V4)k5<-k5+1
  
  
  
  #--------------------------V5---------------------------------------------
  #V5=m2[10]
  
  MO_V5<-mean(V5)               # Mean Occurance
  sO_V5<-sd(V5)                 # SD of Occurance
  
  z_V5<-(m2[10]-MO_V5)/sO_V5;    # Z score   
  
  
  Zrand_V5<-(V5-MO_V5)/sO_V5
  UQ_V5<-quantile(Zrand_V5, 0.975)        # 1.256084
  LQ_V5<-quantile(Zrand_V5, 0.025)        # -2.264109 
  
  #if(abs(z_V5)>UQ_V4)k6<-k6+1
  
  
  
  #print(j)
  
#}








#----------------------------------Table---------------------------------


  Motif<-factor(c("AVPL", "Diameter","T2","V3","V4","V5"))
  
  Ob_Occur<-round(c(AVPL_G,Diameter_G,m1[4],m2[8],m2[9],m2[10]),3)
  
  Mean<-c(MO_ASPL,MO_Dim,MO_T2,MO_V3,MO_V4,MO_V5)
  sd1<-c(sO_ASPL,sO_Dim,sO_T2,sO_V3,sO_V4,sO_V5)
  
  Zm<-round(c(z_AVPL11, z_Dim11,z_T2,z_V3,z_V4,z_V5),3)
  LQA<-round(c(LQ_AVPL11,LQ_Dim11,LQ_T2,LQ_V3,LQ_V4,LQ_V5),3)
  UQA<-round(c(UQ_AVPL11,UQ_Dim11,UQ_T2,UQ_V3,UQ_V4,UQ_V5),3)
  
  
  data.frame ('Stat'=Motif, 'Observed'=Ob_Occur, 'Mean'=Mean,'s'=sd1,'Zm'=Zm,'LL'=LQA,'UL'=UQA) 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
############  ER ###########################

  #      Stat Observed      Mean          s     Zm     LL    UL
  #     AVPL    6.309   4.211095  0.1244334 16.857 -2.019 2.055
  #  Diameter   14.000  9.418500  0.9585495  4.780 -1.480 2.693
  #       T2   23.000  4.573000  2.2176857  8.309 -1.611 1.996
  #       V3  132.000 38.205500 19.4780637  4.815 -1.602 2.248
  #       V4   20.000  9.581000  3.4392373  3.029 -1.623 2.157
  #       V5    5.000  0.520000  0.8584687  5.219 -0.606 2.889
  
  
############  PA ##############


  #     Stat Observed       Mean          s     Zm     LL    UL
  #     AVPL    6.309   4.097325  0.1093167 20.229 -1.920 2.066
  # Diameter   14.000   9.234000  0.8767145  5.436 -1.408 2.014
  #       T2   23.000  11.464500  3.1663469  3.643 -1.726 2.064
  #       V3  132.000 184.656500 57.4637150 -0.916 -1.752 2.095
  #       V4   20.000  25.142000  6.3374220 -0.811 -1.758 2.187
  #       V5    5.000   7.013500  4.5068252 -0.447 -1.334 2.438













