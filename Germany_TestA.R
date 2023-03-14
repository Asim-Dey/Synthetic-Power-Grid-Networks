
#install.packages("NetSwan")
library(igraph)
library(NetSwan)


data11 <- read.csv("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Data/Export_Output22.csv")
data11 <- read.csv("C:/Users/adey/OneDrive//Synthentic Network/Random Graph Model/R codes/2020/Old/Data/Export_Output22.csv")


source('C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Pref_at_G.R')
source('C:/Users/adey/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Pref_at_G.R')


##############################German Grid##############################

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



#-------------------------------------Remove duplicity -----------------------------------------------------
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

#####################################################################################
#-----------------------Statistical Properties --------------------------------------

G0<-new_Germany_network

igraph.options(vertex.size=2,  edge.arrow.size=0.9,vertex.label=NA)#vertex.label=NA
plot(G0,  vertex.size=2)

#n0=1000 # node
#m0=1300 # edge

#G0<-PrefAt(n0,m0) # Preferential Attachment Model with equal number of node and edge 

#G0<-sample_gnm(n0,m0,directed = FALSE, loops = FALSE) #  Erdos-Renyi graphs G(n,m)





#---Degree Distribution--------------

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
C0<-c(AVPL_G,Diameter_G);C0
M0<-c(m1[4],m2[8],m2[9],m2[10]);M0


#######################################################################################
#################################### eNN ################################################


G2 <- read.table("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Data/Test_Germany_eNN_Motifs.txt", header=FALSE)[1:2000,]
G2 <- read.table("C:/Users/adey/OneDrive/Synthentic Network/Random Graph Model/Data/Test_Germany_eNN_Motifs.txt", header=FALSE)[1:2000,]

names(G2)<-c('nodes','edges','Diameter','AVPL','triangles','square','k13','tent','kite','k4')

head(G2)
dim(G2)




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


#    Stat   Observed    Mean           s      Zm     LL    UL
#     AVPL   11.989    7.62529   0.2677859  16.294 -1.924 1.735
# Diameter   31.000   16.39050   1.1487687  12.718 -2.951 1.401
#       T2   29.000  231.05300  18.6603938 -10.828 -2.468 1.444
#       V3  165.000 1123.57200 126.6197777  -7.570 -2.042 1.567
#       V4   22.000   53.02200  11.0334212  -2.812 -2.540 1.720
#       V5    9.000  169.38400  25.8283584  -6.210 -2.260 1.727





#######################################################################################
#################################### CLC 1 param ################################################


G2 <- read.table("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Data/Test_Germany_1parar_Motifs.txt", header=FALSE)
G2 <- read.table("C:/Users/adey/OneDrive/Synthentic Network/Random Graph Model/Data/Test_Germany_1parar_Motifs.txt", header=FALSE)

names(G2)<-c('nodes','edges','Diameter','AVPL','triangles','square','k13','tent','kite','k4')

head(G2)
dim(G2)




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

#     Stat   Observed   Mean      s        Zm     LL    UL
#     AVPL   11.989  13.87087  1.252714 -1.502 -1.342 1.660
# Diameter   31.000  36.18050  5.381301 -0.963 -1.706 1.639
#       T2   29.000  57.08700  9.727371 -2.887 -1.962 2.253
#       V3  165.000 286.86300 46.829560 -2.602 -1.855 1.690
#       V4   22.000  37.27100  7.661521 -1.993 -1.732 2.314
#       V5    9.000  25.28150  8.059264 -2.020 -1.648 2.447












#######################################################################################
#################################### CLC 2 param ################################################


G2 <- read.table("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Data/Test_Germany_2parar_Motifs.txt", header=FALSE)
G2 <- read.table("C:/Users/adey/OneDrive/Synthentic Network/Random Graph Model/Data/Test_Germany_2parar_Motifs.txt", header=FALSE)

names(G2)<-c('nodes','edges','Diameter','AVPL','triangles','square','k13','tent','kite','k4')

head(G2)
dim(G2)




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

#   Stat     Observed   Mean      s       Zm     LL    UL
#     AVPL   11.989  10.90248  0.8666369 1.253 -1.503 1.647
# Diameter   31.000  28.06950  3.9704599 0.738 -1.781 1.746
#       T2   29.000  25.46250  5.1655528 0.685 -1.832 2.040
#       V3  165.000 129.10450 26.7020613 1.344 -1.802 2.131
#       V4   22.000  18.75550  4.9244131 0.659 -1.778 2.085
#       V5    9.000   7.29050  3.4791036 0.491 -1.808 2.216











#######################################################################################
#################################### GeoDe ################################################

GG1 <- read.csv("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Data/German_3000.txt", header=FALSE)
GG2 <- read.csv("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Data/German_3000B.txt", header=FALSE)


GG1 <- read.csv("C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/Data/German_3000.txt", header=FALSE)
GG2 <- read.csv("C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/Data/German_3000B.txt", header=FALSE)



G2<-rbind(GG1,GG2)
dim(G2)

names(G2)<-c("m", "t","nodes", "edges", "vee", "triangles", "path3", "square", 
             "k13",  "tent", "kite", "k4","Diameter", "AVPL")

head(G2)
dim(G2)

#############################################################################################
####################################################################################################

nrep=200

k1<-0;k2<-0;k3<-0;k4<-0
k5<-0;k6<-0;k7<-0


#for(j in 1:nrep){  #j=1
#set.seed(123)

G1<-G2[sample(nrow(G2), 1500), ] # 200
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


#     Stat Observed       Mean         s     Zm     LL    UL
#     AVPL   11.989  11.783317  1.052763  0.195 -1.538 2.322
# Diameter   31.000  28.901333  4.073655  0.515 -1.449 2.479
#       T2   29.000  30.004000  6.253483 -0.161 -1.920 2.078
#       V3  165.000 148.540667 34.905481  0.472 -1.749 2.148
#       V4   22.000  22.360000  6.348433 -0.057 -1.632 2.306
#       V5    9.000   5.716667  3.122641  1.051 -1.510 2.332

 
################################# ER and PA ###############################################
################################ Simulation ################################################################


#nrep=200

#k1<-0;k2<-0;k3<-0;k4<-0
#k5<-0;k6<-0;k7<-0


#for(j in 1:nrep){  #j=1
#set.seed(123)



n=n_node # node
m=n_edge # edge

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




############ ER #######################


#    Stat Observed      Mean          s     Zm     LL    UL
#     AVPL   11.989  6.248849  0.1487396 38.590 -1.919 1.977
# Diameter   31.000 14.689000  1.3206631 12.351 -1.279 2.507
#       T2   29.000  2.787000  1.6712353 15.685 -1.668 1.937
#       V3  165.000 20.824000 13.2140202 10.911 -1.576 2.359
#       V4   22.000  5.189500  2.3998895  7.005 -1.746 2.004
#       V5    9.000  0.055000  0.2449592 36.516 -0.225 3.858





############ PA ######################################

#     Stat Observed       Mean          s     Zm     LL    UL
#     AVPL   11.989   5.628892  0.1135115 56.028 -1.777 2.084
# Diameter   31.000  13.465500  1.0589478 16.558 -1.384 2.393
#       T2   29.000  14.750000  3.7072054  3.844 -1.821 2.225
#       V3  165.000 338.680000 96.9859822 -1.791 -1.729 2.189
#       V4   22.000  37.893500  7.7753705 -2.044 -1.658 2.071
#       V5    9.000   8.190000  5.0841739  0.159 -1.414 2.323







