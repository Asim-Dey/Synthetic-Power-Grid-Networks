#install.packages("NetSwan")
library(igraph)
library(NetSwan)

##############################Italian Power Grid system ##############################




data11 <- read.csv("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Data/Export_Output22.csv")
data11 <- read.csv("C:/Users/adey/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Data/Export_Output22.csv")


source('C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Pref_at_G.R')
source('C:/Users/adey/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Old/Pref_at_G.R')


nodes_Italy<-data11[,c(2,3,4)][data11[,5]=="Italy",]

#write.csv(nodes_Italy,"C:/Users/akd130230/Dropbox/Power Grid Network/Data from author/processed/Italy_nodes.csv")

ID_Italy<-data11[,2][data11[,5]=="Italy"]

#data11[,17] # Fnode
#data11[,18] # Tnode

#Italy<-data11[(data11[,17] %in% ID_Italy) & (data11[,18] %in% ID_Italy)  , ]

Italy<-data11[(data11[,17] %in% ID_Italy) | (data11[,18] %in% ID_Italy)  , ] # or
#write.csv(Italy,"C:/Users/akd130230/Dropbox/Power Grid Network/Data from author/processed/Italy_edge.csv")

Italy_data_all<-Italy[,c(16,17,18)]



#-------------------------------------Remove duplicity -----------------------------------------------------

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

#plot(new_Italy_network, layout=layout_with_fr, vertex.size=5)

#####################################################################################
#-----------------------Statistical Properties --------------------------------------

G0<-new_Italy_network

igraph.options(vertex.size=2, vertex.color='black', edge.width=1, edge.color='gray',edge.arrow.size=0.9,
               vertex.label=NA)#vertex.label=NA
plot(G0)



#---Degree Distribution--------------

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


######################################################################################
#################################### eNN ################################################


G2 <- read.table("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Data/Test_Italy_eNN_Motifs.txt", header=FALSE)[1:2000,]
G2 <- read.table("C:/Users/adey/OneDrive/Synthentic Network/Random Graph Model/Data/Test_Italy_eNN_Motifs.txt", header=FALSE)[1:2000,]

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


#     Stat   Observed   Mean         s      Zm     LL    UL
#     AVPL    9.971   6.67714  0.2879699 11.439 -1.796 2.545
# Diameter   28.000  14.52400  1.3846959  9.732 -1.101 1.788
#       T2   23.000 161.29100 14.3869319 -9.612 -2.244 2.065
#       V3  133.000 765.88400 92.2398267 -6.861 -1.809 2.365
#       V4   29.000  35.39850  7.5172437 -0.851 -1.649 2.075
#       V5    6.000 121.89800 19.6512514 -5.898 -1.878 1.990





######################################################################################
#################################### CLC 1 param ################################################


G2 <- read.table("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Data/Test_Italy_1parar_Motifs.txt", header=FALSE)[1:2000,]
G2 <- read.table("C:/Users/adey/OneDrive/Synthentic Network/Random Graph Model/Data/Test_Italy_1parar_Motifs.txt", header=FALSE)[1:2000,]

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



#     stat    Observed  Mean        s       Zm     LL    UL
#     AVPL    9.971  11.26477  0.8997255 -1.438 -1.761 1.384
# Diameter   28.000  29.58550  4.1532294 -0.382 -2.067 1.544
#       T2   23.000  32.33600  6.4879234 -1.439 -1.747 1.798
#       V3  133.000 151.72050 27.5371300 -0.680 -1.842 1.971
#       V4   29.000  20.13150  4.4712642  1.983 -1.819 1.760
#       V5    6.000  12.59650  5.1657555 -1.277 -1.471 1.820







#######################################################################################
#################################### CLC 2 param ################################################


G2 <- read.table("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Data/Test_Italy_2parar_Motifs.txt", header=FALSE)
G2 <- read.table("C:/Users/adey/OneDrive/Synthentic Network/Random Graph Model/Data/Test_Italy_2parar_Motifs.txt", header=FALSE)

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

#      Stat Observed  Mean      s        Zm     LL    UL
#     AVPL    9.971  9.82726  0.7987042 0.180 -1.524 1.706
# Diameter   28.000 25.13450  3.7044934 0.774 -1.656 1.853
#       T2   23.000 16.30250  3.9308675 1.704 -1.858 1.958
#       V3  133.000 74.48700 18.3034789 3.197 -1.884 2.159
#       V4   29.000 12.95350  4.0354033 3.976 -1.723 2.242
#       V5    6.000  3.71800  2.3061082 0.990 -1.612 2.290




#######################################################################################
#################################### GeoDe ################################################

GG1 <- read.csv("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Data/Italy_6000C.txt", header=FALSE)
GG2 <- read.csv("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Data/Italy_6000D.txt", header=FALSE)
GG3 <- read.csv("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Data/Italy_6000D1.txt", header=FALSE)




GG1 <- read.csv("C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/Data/Italy_6000C.txt", header=FALSE)
GG2 <- read.csv("C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/Data/Italy_6000D.txt", header=FALSE)
GG3 <- read.csv("C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/Data/Italy_6000D1.txt", header=FALSE)





G2<-rbind(GG1,GG2,GG3)
dim(G2)




names(G2)<-c("m", "t","nodes", "edges", "vee", "triangles", "path3", "square", 
             "k13",  "tent", "kite", "k4","Diameter", "AVPL")

head(G2)
dim(G2)

#write.csv(dd,"C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/A/Data/ER_Italy.csv")



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



#    Stat Observed      Mean          s     Zm     LL    UL
#     AVPL    9.971   9.19660  0.8107862  0.956 -1.535 2.290
# Diameter   28.000  22.17533  3.1831871  1.830 -1.312 2.458
#       T2   23.000  24.43200  5.7299736 -0.250 -1.821 2.193
#       V3  133.000 122.94000 32.1840998  0.313 -1.786 2.115
#       V4   29.000  17.98800  5.6183835  1.960 -1.778 2.138
#       V5    6.000   5.68400  3.4337547  0.092 -1.364 2.131











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


#   Stat Observed      Mean          s     Zm     LL    UL
#     AVPL    9.971  5.385842  0.1339493 34.233 -1.871 2.004
# Diameter   28.000 12.418000  1.1467055 13.588 -1.237 2.252
#       T2   23.000  3.387000  1.8544407 10.576 -1.826 1.948
#       V3  133.000 26.949000 15.8079203  6.709 -1.705 2.281
#       V4   29.000  6.910000  2.7354417  8.075 -1.795 2.226
#       V5    6.000  0.135500  0.4089295 14.341 -0.331 2.114




############# PA ######################################

#     Stat Observed      Mean          s     Zm     LL    UL
#     AVPL    9.971   4.99613  0.1072302 46.397 -1.903 2.070
# Diameter   28.000  11.72100  1.0003296 16.274 -1.720 2.278
#       T2   23.000  14.01000  3.5643025  2.522 -1.967 2.242
#       V3  133.000 285.69850 81.5475849 -1.873 -1.725 2.199
#       V4   29.000  34.03450  7.2663409 -0.693 -1.794 2.060
#       V5    6.000   8.01250  5.0062836 -0.402 -1.401 2.195  











