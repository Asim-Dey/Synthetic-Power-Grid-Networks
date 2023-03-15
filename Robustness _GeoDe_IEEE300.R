

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
               vertex.label=NA)#vertex.label=NA
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

g<-G0 # G0, G1, G_ER, G_PA 

c0 = components(g)
GC_0<-max(c0$csize);GC_0

AVPL_0<-average.path.length(g) 
#---------------------------------------------------------------------

n<-length(V(g))

fr<-numeric(n)
GC<-numeric(n)
AVPL<-numeric(n)


n1<-numeric(n);n2<-numeric(n)
#----------------------------------degree based attack------------------
mat<- matrix(ncol=2,nrow=n, 0) 
mat[,1]<-1:n

deg<-degree(g)                  # Degree based attacks
#deg<-betweenness(g)              # Betweenness based attacks
mat[,2]<-deg
matri<-mat[order(mat[,2]),] 



g2<-g
#-----------------------------------------------------------------------------

for(i in 1:n){       #i=5
  
  v=n+1-i
  g2<-delete_vertices(g2, matri[v,1])
  
  
  c = components(g2)    #largest Connected components/Giant component 
  GC[i]<-max(c$csize) 
  
  AVPL[i]<-average.path.length(g2) # Average pacth length
  
  
  fr[i]<-i/n   # fraction of node is removes
  matri[matri[,1]>matri[v,1],1]<-matri[matri[,1]>matri[v,1],1]-1 
  
  
  }




Giant_Comp_APL<-data.frame(c(0,fr),c(GC_0,GC),c(AVPL_0,AVPL))
colnames(Giant_Comp_APL) <- c("fr", "GC","APL") 


#################################### Plots #####################################


d0<- read.csv("Giant_fr_IEEE300Original.csv")
head(d0)

d0_GeoDE<- read.csv("Giant_fr_IEEE300GeoDe.csv")
head(d0_GeoDE)





plot(d0$fr,d0$GC/d0$GC[1], type="l", lty=4, lwd=2, col='black',xlab="Fraction of nodes removed",
     ylab ="Normalized Giant Component" ,ylim=c(),xlim=c(0.0,0.4),first.pannel=grid())

lines(d0_GeoDE$fr,d0_GeoDE$GC/d0_GeoDE$GC[1],  pch = 17,lty=1,lwd=2,col='red') 
#lines(d0_ER$fr,d0_ER$GC/d0_ER$GC[1], type='l', pch = 19, lty=1,lw=2,col='green') 
#lines(d0_PA$fr,d0_PA$GC/d0_PA$GC[1], type='l', pch = 15,lty=1,lwd=2,col='yellow')


legend('topright',c("Observed", "GeoDe"), title="IEEE 300",
       lty=c(4,1),  lwd=c(2,2),  pch = c(NA,NA),col=c("black", "red")) 


