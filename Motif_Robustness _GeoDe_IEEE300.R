

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


G0<-Final_IEEE_300_network

igraph.options(vertex.size=2, vertex.color='black', edge.width=1, edge.color='gray',edge.arrow.size=0.9,
               vertex.label=NA)# vertex.label=NA
plot(G0)

n_node<-length(V(G0));n_node
n_edge<-length(E(G0));n_edge

#################################### GeoDe Instance ################################################


edge_list <- read.table("edgelistIEEE300_1.txt")


ED<-as.matrix(edge_list)+1 #  vertex of name start with 0. So we add 1

G1<-graph_from_edgelist(ED,directed = FALSE)

igraph.options(vertex.size=2, vertex.color='black', edge.width=1, edge.color='gray',edge.arrow.size=0.9,
               vertex.label=NA)#vertex.label=NA
plot(G1)


n_node1<-length(V(G1));n_node1
n_edge1<-length(E(G1));n_edge1


################################# ER and PA ###############################################
################################ Simulation ################################################################



n=n_node
m=n_edge


G_ER<- sample_gnm(n,m,directed = FALSE, loops = FALSE) #  Erdos-Renyi graphs G(n,m)

G_PA<-PrefAt(n,m) # Preferential Attachment Model with equal number of node and edge 



igraph.options(vertex.size=2, vertex.color='black', edge.width=1, edge.color='gray',edge.arrow.size=0.9,
               vertex.label=NA)#vertex.label=NA
plot(G_ER)
plot(G_PA)



##################################################################################
##################################################################################
######################### Robustness #################################

g<-G1 # G0, G1, G_ER, G_PA 



node<-length(V(g));node
edge<-length(E(g));edge

#-------- Motif, size=4 ---------------
m02<-motifs(g, 4)
m02[is.na(m02)] <- 0
n22<-count_motifs(g, 4);n22

V1_0<-m02[5];C_V1_0<-m02[5]/n22; C_V1_0[is.na(C_V1_0)] <-0
V2_0<-m02[7];C_V2_0<-m02[7]/n22; C_V2_0[is.na(C_V2_0)] <-0
V3_0<-m02[8];C_V3_0<-m02[8]/n22;C_V3_0[is.na(C_V3_0)] <-0
V4_0<-m02[9];C_V4_0<-m02[9]/n22;C_V4_0[is.na(C_V4_0)] <-0
V5_0<-m02[10];C_V5_0<-m02[10]/n22;C_V5_0[is.na(C_V5_0)] <-0
V6_0<-m02[11];C_V6_0<-m02[11]/n22;C_V6_0[is.na(C_V6_0)] <-0

Tot_V0<-sum(V1_0,V2_0,V3_0,V5_0,V6_0)
#-----------------------------------------------
n<-length(V(g))
noeud<-V(g);  arc<-E(g)


fr<-numeric(n)
V1<-numeric(n);V2<-numeric(n);V3<-numeric(n);V4<-numeric(n)
V5<-numeric(n);V6<-numeric(n);T1<-numeric(n);T2<-numeric(n)
Tot_V<-numeric(n)

C_V1<-numeric(n);C_V2<-numeric(n);C_V3<-numeric(n);C_V4<-numeric(n)
C_V5<-numeric(n);C_V6<-numeric(n);C_T1<-numeric(n);C_T2<-numeric(n)

n1<-numeric(n);n2<-numeric(n)
#----------------------------------degree based attack-------------------------------
mat<- matrix(ncol=2,nrow=n, 0) 
mat[,1]<-1:n

deg<-degree(g)
mat[,2]<-deg
matri<-mat[order(mat[,2]),] 




g2<-g

#-------------------------------------------------------------------------------

for(i in 1:n){       #i=5
  
  v=n+1-i
  g2<-delete_vertices(g2, matri[v,1])
  # plot(g2, layout=layout_with_fr, vertex.size=5)
  
  
  #------Motif size=4 ---------------------------------------------------
  m2<-motifs(g2, 4)
  m2[is.na(m2)] <- 0
  
  n02<-count_motifs(g2, 4)
  n2[i]<- n02
   
  V1[i]<-m2[5];C_V1[i]<-m2[5]/n02
  
  V2[i]<-m2[7];C_V2[i]<-m2[7]/n02
  
  V3[i]<-m2[8];C_V3[i]<-m2[8]/n02
  V4[i]<-m2[9];C_V4[i]<-m2[9]/n02
  V5[i]<-m2[10];C_V5[i]<-m2[10]/n02
  V6[i]<-m2[11];C_V6[i]<-m2[11]/n02
  
  Tot_V[i]<-sum(m2[5],m2[7],m2[8],m2[9],m2[10],m2[11])
  #-------------------------------------------------
  fr[i]<-i/n
  
  matri[matri[,1]>matri[v,1],1]<-matri[matri[,1]>matri[v,1],1]-1 
}

C_V1<-V1/n22;C_V1[is.na(C_V1)] <-0
C_V2<-V2/n22;C_V2[is.na(C_V2)] <-0
C_V3<-V3/n22;C_V3[is.na(C_V3)] <-0 
C_V4<-V4/n22;C_V4[is.na(C_V4)] <-0
C_V5<-V5/n22;C_V5[is.na(C_V5)] <-0

Time<-0:n 

motif_degree<-data.frame(Time,c(0,fr),c(Tot_V0,Tot_V),c(V1_0,V1),c(V2_0,V2),c(V3_0,V3),c(V4_0,V4),c(V5_0,V5),c(C_V1_0,C_V1),c(C_V2_0,C_V2),c(C_V3_0,C_V3),c(C_V4_0,C_V4),c(C_V5_0,C_V5))

colnames(motif_degree) <- c("Time","fr","Tot_V", "V1", "V2", "V3", "V4", "V5", "C_V1", "C_V2", "C_V3", "C_V4", "C_V5") 





#################################### Plots #####################################


d0<- read.csv("MotifIEEE300_original.csv")
head(d0)

d0_GeoDE<- read.csv("MotifIEEE300_GeoDe.csv")
head(d0_GeoDE)




 


plot(d0$fr,d0$C_V3, type='l', lty=4,lwd=2, col='black',xlab="Fraction of nodes removed",ylim=c(),
     xlim=c(0,0.2), ylab="Concentration",first.pannel=grid())    

#lines(d0$fr,d0$C_V5, type='l',lty=4,lwd=2,  col='black')  

lines(d0_GeoDE$fr,d0_GeoDE$C_V3, type='l',lty=1,lwd=2,  col='red')  
#lines(d0_GeoDE$fr,d0_GeoDE$C_V5, type='l',lty=4,lwd=2,  col='red') 


lines(d0$fr,d0$C_V4, type='l',lty=4,lwd=2,  col='green')  
lines(d0_GeoDE$fr,d0_GeoDE$C_V4, type='l',lty=1,lwd=2,  col='blue') 



lLab<-c(expression("Observed M"[3]),expression("GeoDe M"[3]), expression("Observed M"[4]), expression("GeoDe M"[4]))
legend('topright', title="IEEE 300",legend=lLab,
       lty=c(4,1,4,1),  lwd=c(2,2,2,2),
       col=c("black","red","green","blue"))













