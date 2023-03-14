


#install.packages("NetSwan")
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
#-----------------------Statistical Properties --------------------------------------

G0<-Final_IEEE_300_network


igraph.options(vertex.size=2, vertex.color='black', edge.width=1, edge.color='gray',edge.arrow.size=0.9,
               vertex.label=NA)#vertex.label=NA
plot(G0)


#options = {'node_color': 'black',  'node_size': 10, 'width': 0.5,'edge_color': 'gray'}
#nx.draw(G)
#nx.draw_random(G, **options)


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
C0<-c(AVPL_G,Diameter_G)
M0<-c(m1[4],m2[8],m2[9],m2[10])
m2[11]

T2<-m1[4]
c(n_node,n_edge,T2)



##########################################################################################################
########################################## Node Degree Distribution  ###########################################



degree1<-degree(G0)
summary(degree1)
table(degree1)

hist(degree1,breaks=10,first.pannel=grid(),xlab="Degree",main="")
box()

# 1   2  3  4  5  6  7  8  9 11 
# 69 76 84 42 14  6  5  2  1  1 

meanD1<-mean(degree1);meanD1 



degree_Germany_N<-degree.distribution(G0)[2:12]
degree_Germany<-degree_Germany_N*300 ;degree_Germany #Node 445


#--------Poisson degree distribution----------------

degree_of_Germany=as.numeric(degree1) 

library(MASS)
m1=fitdistr(degree_of_Germany,densfun = "Poisson");m1 #lambda: 2.54831461 : Mean=lambda

m2=fitdistr(degree_of_Germany,densfun = "exponential");m2 #gamma:  0.39241623  : Mean=(1/gamma)=2.548315



plot(dpois(x=1:11, lambda=2.54831461), type="b", ylab = "p(k)",xlab = "Number of edge (Degree), k",
     main = "",ylim = c(0,0.4),col="red",first.pannel=grid())

lines(degree_Germany_N,type="o", col="black")
lines(dexp(x=1:11, rate=0.54447115),type="o", col="blue")

#plot(degree_Germany_N,type="o", col="black")


legend("topright",  #cex=1.1, y.intersp=.9,
       c('Observed degree distribution ', 'Poisson degree distribution','Exponential degree distribution'),
       lwd=c(1,1,1), 
       lty=c(1,1,1),
       col=c('black',  'red','blue')) 



########################################################################################

degree_Germany_N
# 0.230000000 0.253333333 0.280000000 0.140000000 0.046666667 0.020000000 0.016666667
# 0.006666667 0.003333333 0.000000000 0.003333333


sim_1<-c(0, 44, 95, 82, 46, 27, 4, 2, 0, 1,0,0) # 0:10--0:12
sim_2<-c(0, 45, 96, 79, 46, 20, 8, 4, 1,0,0,0);0:9
sim_3<-c(0, 57, 85, 65, 53, 23, 7, 1, 1,0,0,0);0:9
sim_4<-c(0, 52, 96, 73, 57, 19, 4, 4, 2,0,0,0);0:9
sim_5<-c(0, 56, 71, 85, 55, 24, 12, 3, 2,0,0,0); 0:9
sim_6<-c(0, 58, 85, 80, 48, 18, 5, 2, 1,0,0,0);0:9
sim_7<-c(0, 54, 112, 69, 38, 13, 7, 4, 1,0,0,0);0:9
sim_8<-c(0, 44, 95, 68, 56, 20, 10, 2, 2,0,0,0);0:9
sim_9<-c(0, 53, 75, 79, 60, 28, 9, 1, 1, 1,0,0);0:10
sim_10<-c(0, 63, 85, 78, 41, 22, 13, 3, 1,0,0,0);0:9


#density1<- sim_1/sum(sim_1)
#density1
#density2<-c(density1,0)


plot((sim_10/sum(sim_10))[2:12], type="o", ylab = "p(k)",xlab = "Degree, k",
     main = "",ylim = c(0,0.4),col="black",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
     first.pannel=grid())



lines((sim_1/sum(sim_1))[2:12],type="o", col="black")
lines((sim_2/sum(sim_2))[2:12],type="o", col="black")

lines((sim_3/sum(sim_3))[2:12],type="o", col="black")
lines((sim_4/sum(sim_4))[2:12],type="o", col="black")
lines((sim_5/sum(sim_5))[2:12],type="o", col="black")
lines((sim_6/sum(sim_6))[2:12],type="o", col="black")
lines((sim_7/sum(sim_7))[2:12],type="o", col="black")
lines((sim_8/sum(sim_8))[2:12],type="o", col="black")

lines((sim_9/sum(sim_9))[2:12],type="o", col="black")

lines(degree_Germany_N,type="o", col="red",lwd=1)

legend("topright",  #cex=1.1, y.intersp=.9,
       c('Observed', 'GeoDe'),
       lwd=c(1,1), 
       lty=c(1,1),
       col=c('red','black'),cex=1.5) 




#-----------------------------------





m1=fitdistr(degree_of_Germany,densfun = "Poisson");m1 #lambda:   2.72666667  : Mean=lambda
m2=fitdistr(degree_of_Germany,densfun = "exponential");m2 #gamma:    0.36674817   : Mean=(1/gamma)=2.548315



plot(dpois(x=1:11, lambda=2.72666667), type="b", ylab = "p(k)",xlab = "Degree, k",
     main = "",ylim = c(0,0.4),col="red",first.pannel=grid())

lines(degree_Germany_N,type="o", col="black")
lines((sim_1/sum(sim_1))[2:12],type="o", col="green")

lines(dexp(x=1:11, rate=0.36674817),type="o", col="blue")



legend("topright",  #cex=1.1, y.intersp=.9,
       c('Observed', 'Poisson','Exponential', 'GeoDe'),
       lwd=c(1,1,1,1), 
       lty=c(1,1,1,1),
       col=c('black',  'red','blue','green')) 


################## KL divergence #######################
# X to Y (Q to P)
sim_A0<-degree.distribution(g1)
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
sim_A0<-degree.distribution(g1)
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






########################################################################################


################################ Simulation ################################################################

#source('C:/Users/asimi/Downloads/R codes/2020/Pref_at_G.R')
set.seed(123)

n=n_node
m=n_edge

B=500


M<-numeric(B)
AVPL<-numeric(B)
Diameter<-numeric(B)

V1<-numeric(B);V2<-numeric(B);V3<-numeric(B);V4<-numeric(B)
V5<-numeric(B);V6<-numeric(B);T1<-numeric(B);T2<-numeric(B)


C_V1<-numeric(B);C_V2<-numeric(B);C_V3<-numeric(B);C_V4<-numeric(B)
C_V5<-numeric(B);C_V6<-numeric(B);C_T1<-numeric(B);C_T2<-numeric(B)

for(i in 1:B) {  #i=1
  
  #g1<-degree.sequence.game(deg_G0, method="simple.no.multiple")
  g1<- sample_gnm(n,m,directed = FALSE, loops = FALSE) #  Erdos-Renyi graphs G(n,m)
  
  #g1<-PrefAt(n,m) # Preferential Attachment Model with equal number of node and edge 

  igraph.options(vertex.size=2, vertex.color='black', edge.width=1, edge.color='gray',edge.arrow.size=0.9,
                 vertex.label=NA) #vertex.label=NA
  plot(g1)
  length(V(g1))
  length(E(g1))

  
 
  
  
  
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
  
  
}





   

dd<-data.frame(AVPL,Diameter,T2,V3,V4,V5)

#write.csv(dd,"C:/Users/asimi/Downloads/R codes/2020/PA_IEEE300.csv")
#write.csv(dd,"C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/Data/D0ER.csv")

D0<-read.csv("C:/Users/asimi/Downloads/R codes/2020/D0ER.csv")
attach(D0)




 

######################################################################################################################

D0PfA<-read.csv("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/A/Data/PA_IEEE300.csv")
D0PfA<-read.csv("C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/A/Data/PA_IEEE300.csv")
dim(D0PfA)

D0_ER<-read.csv("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/A/Data/ER_IEEE300.csv")
D0_ER<-read.csv("C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/A/Data/ER_IEEE300.csv")
dim(D0_ER)

D0GeoDe<- read.csv("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/B/Data/IEEE300_200B.txt", header=FALSE)
D0GeoDe<- read.csv("C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/B/Data/IEEE300_200B.txt", header=FALSE)
dim(D0GeoDe)

names(D0GeoDe)<-c("m", "t","nodes", "edges", "vee", "triangles", "path3", "square", 
                  "k13",  "tent", "kite", "k4","Diameter", "AVPL")









######################################################################


#C0<-c(AVPL_G,Diameter_G)
#T2, V3,V4,V4
#M0<-c(m1[4],m2[8],m2[9],m2[10])

X1<-D0PfA[1:200,]$V5
X2<-D0_ER[1:200,]$V5
X3<-D0GeoDe[1:200,]$kite


plot.ts(X1, ylim=c(0,25),ylab="M5/Kite",first.pannel=grid(),xlab="Iterations")
lines(X2,col="red")
lines(X3,col="green")
abline(h=m2[10], col="blue",lty=2,lwd=2) 

legend('topright', #inset = c(0,0.1),
       #title="Diameter",
       c("IEEE 300", "PA","ER","DeoGe"),
       lty=c(2,1,1,1),lwd=c(2,1,1,1),
       col=c("blue", "black","red","green"),cex=0.75)  


 
################################################################################



# two sample mean test
t<-t.test(X3,mu=m2[10])
t














#par(mfrow=c(1,2))
plot(ecdf(X1) ,cex=0.5,col="black",ylab="Cumulative density",xlim=c(0,60),
     xlab="V5/Kite",main="",first.pannel=grid())

lines(ecdf(X2),cex=0.5, col="red")
lines(ecdf(X3),cex=0.5, col="green")

#grid()

legend('topright', #inset = c(0,0.1),
       #title="Diameter",
       c("FA","ER","DeoGe"),
       lty=c(1,1,1),
       col=c("black","red","green"),cex=0.75)   





# two sample mean test
t<-t.test(X3,X1)
t




#Two-sample Kolmogorov-Smirnov test
ks<-ks.test(X1,X3)
ks


# two sample mean test
t<-t.test(X1,X3)
t



################## Multivariate motif analysis ######################################


#C0<-c(AVPL_G,Diameter_G)
#T2, V3,V4,V4
#M0<-c(m1[4],m2[8],m2[9],m2[10])


X1_PfA<-D0PfA[1:200,]$T2
X2_PfA<-D0PfA[1:200,]$V3
X3_PfA<-D0PfA[1:200,]$V4
X4_PfA<-D0PfA[1:200,]$V5


X1_ER<-D0_ER[1:200,]$T2
X2_ER<-D0_ER[1:200,]$V3
X3_ER<-D0_ER[1:200,]$V4
X4_ER<-D0_ER[1:200,]$V5


X1_GeoDe<-D0GeoDe[1:200,]$triangles
X2_GeoDe<-D0GeoDe[1:200,]$tent
X3_GeoDe<-D0GeoDe[1:200,]$square
X4_GeoDe<-D0GeoDe[1:200,]$kite


#------------------------ boxplot---------------------

DF2 <- data.frame(
  x = c(c(X1_ER,X1_PfA, X1_GeoDe), c(X2_ER, X2_PfA, X2_GeoDe), c(X3_ER, X3_PfA, X3_GeoDe),
        c(X4_ER, X4_PfA, X4_GeoDe)),
  y = rep(c("T2", "M3","M4","M5"), each = 3*length(X2_ER)),
  z = rep(rep(c("ER", "PA","GeoDe"), each=length(X2_ER)), 4),
  stringsAsFactors = FALSE
)


library(ggplot2)
ggplot(DF2, aes(y, x, fill=factor(z))) +
  geom_boxplot()


ggplot(DF2, aes(y, x)) +
  geom_boxplot(aes(fill = z)) + 
  ggtitle("Google searches")+ 
  theme(plot.title = element_text(hjust = 0.5)) #+ 
#theme_classic()

#########################

DF_T2 <- data.frame(
  x = c(c(X1_ER,X1_PfA, X1_GeoDe)),
  y = rep(c("T2"), each = 3*length(X2_ER)),
  z = rep(rep(c("ER", "PA","GeoDe"), each=length(X1_ER)), 4),
  stringsAsFactors = FALSE
)


DF_M3 <- data.frame(
  x = c(c(X2_ER,X2_PfA, X2_GeoDe)),
  y = rep(c("M3"), each = 3*length(X2_ER)),
  z = rep(rep(c("ER", "PA","GeoDe"), each=length(X2_ER)), 4),
  stringsAsFactors = FALSE
)

DF_M4 <- data.frame(
  x = c(c(X3_ER,X3_PfA, X3_GeoDe)),
  y = rep(c("M4"), each = 3*length(X2_ER)),
  z = rep(rep(c("ER", "PA","GeoDe"), each=length(X3_ER)), 4),
  stringsAsFactors = FALSE
)


DF_M5 <- data.frame(
  x = c(c(X4_ER,X4_PfA, X4_GeoDe)),
  y = rep(c("M5"), each = 3*length(X2_ER)),
  z = rep(rep(c("ER", "PA","GeoDe"), each=length(X4_ER)), 4),
  stringsAsFactors = FALSE
)


 
ggplot(DF_M5, aes(y, x)) +
  geom_boxplot(aes(fill = z)) + 
 # ggtitle(expression(paste("M"[5])))+ 
  theme(plot.title = element_text(hjust = 0.5)) #+ 
#theme_classic()


#------------------------ boxplot---------------------



M0<-c(m1[4],m2[8],m2[9],m2[10]);M0
X_GeoDe<-c(mean(X1_GeoDe),mean(X2_GeoDe),mean(X3_GeoDe),mean(X4_GeoDe))
X_PfA<-c(mean(X1_PfA),mean(X2_PfA),mean(X3_PfA),mean(X4_PfA))
X_ER<-c(mean(X1_ER),mean(X2_ER),mean(X3_ER),mean(X4_ER))


 


plot.ts(M0, col="blue",ylim=c(0,300),lty=2,lwd=1,ylab="",first.pannel=grid(),
        xaxt="no",xlab="")
lines(X_GeoDe,col="green",lwd=1)
lines(X_ER,col="red",lwd=1)
lines(X_PfA,col="black",lwd=1)


legend('topright', #inset = c(0,0.1),
       #title="Diameter",
       c("Observed", "PA","ER","DeoGe"),
       lty=c(2,1,1,1),lwd=c(1,1,1,1),
       col=c("blue", "black","red","green"),cex=0.75)  


axis(1, at=c(1,  2, 3,4), 
     labels= c("T2",  "M3", "M4","M5"),cex=1)




################## KL divergence #######################
# X to Y (Q to P)

p<-(M0/sum(M0));length(p)

#q0<-sim_A;length(q0)
q1<-(X_GeoDe/sum(X_GeoDe));length(q1)


KLD<-sum(p*log(p/q1))
KLD





############ two mean test ############

w1<-wilcox.test(M0, X_PfA)
w1


# two sample mean test
t<-t.test(X3,mu=m2[10])
t


############### Distribution Test ######################


#par(mfrow=c(1,2))
plot(ecdf(M0) ,cex=0.5,col="black",ylab="Cumulative density",xlim=c(0,400),
     xlab="Motifs",main="",first.pannel=grid())

lines(ecdf(X_GeoDe),cex=0.5, col="red")
#grid()

legend('bottomright', #inset = c(0,0.1),
       #title="Diameter",
       c("Observed","GeoDe"),
       lty=c(1,1),   
       col=c("black","red"),cex=0.75)   

# IEEE300_Dist_motifs_comb_PfA




#Two-sample Kolmogorov-Smirnov test
#compute p-values by Monte Carlo simulation, for discrete goodness-of-fit tests only

ks1<-ks.test(M0,X_GeoDe,simulate.p.value=TRUE, B=2000);ks1  # D = 0.25, p-value = 1

ks2<-ks.test(M0,X_PfA,simulate.p.value=TRUE, B=2000);ks2   # D = 0.25, p-value = 1.00

ks3<-ks.test(M0,X_ER,simulate.p.value=TRUE, B=2000);ks3    # D = 0.75, p-value = 0.2286





















