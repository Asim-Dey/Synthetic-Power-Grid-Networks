
#install.packages("NetSwan")
library(igraph) 
library(NetSwan)


data11 <- read.csv("C:/Users/asimi/OneDrive/Power Grid Network/Data from author/processed/Export_Output22.csv")
data11 <- read.csv("D:/02/Dropbox/Power Grid Network/Data from author/processed/Export_Output22.csv")
summary(data11)


country<-c("Italy")
nc<-length(country)


for (j in 1:nc){ #j=1
  
  
  #-------------------------------------------------------------------------------------
  
  nodes_G2<-data11[,c(2,3,4)][data11[,5]==country[j],]
  length(nodes_G2$X)
  ID_G2<-data11[,2][data11[,5]==country[j]]
  
  G2<-data11[(data11[,17] %in% ID_G2) | (data11[,18] %in% ID_G2)  , ] # or
  G2_data_all<-G2[,c(16,17,18)]
  
  G2_data_all_no_dup=unique(G2_data_all)
  G2_data<-G2_data_all_no_dup[,c(2,3)]   # Edge
  G2_edge=data.matrix(G2_data)
  G2_network=graph_from_edgelist(G2_edge,directed = F)
  
  ed11<-c(as.vector(G2_edge[,1]),as.vector(G2_edge[,2]))
  unique_ed11<-unique(ed11)
  nodes_index=sort(unique_ed11)
  m1<-max(nodes_index)
  delete_nodes_index=setdiff(c(1:m1),nodes_index)
  V(G2_network)$name=V(G2_network)
  new_G2_network=delete_vertices(G2_network,delete_nodes_index)
  
  
  
}




G0<-new_G2_network

igraph.options(vertex.size=2, vertex.color='black', edge.width=1, edge.color='gray',edge.arrow.size=0.9,
               vertex.label=NA)#vertex.label=NA
plot(G0)

n_node<-length(V(G0));n_node
n_edge<-length(E(G0));n_edge

#################################### GeoDe Instance ################################################


edge_list <- read.table("C:/Users/akd130230/OneDrive/Synthentic Network/Random Graph Model/Data/edgelistItaly_1.txt")
edge_list <- read.table("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Data/edgelistItaly_1.txt")



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

c0 = components(g)
GC_0<-max(c0$csize);GC_0

AVPL_0<-average.path.length(g) 
#---------------------------------------------------------------------

n<-length(V(g))

fr<-numeric(n)
GC<-numeric(n)
AVPL<-numeric(n)


n1<-numeric(n);n2<-numeric(n)
#----------------------------------Betweeness based attack-----------------
mat<- matrix(ncol=2,nrow=n, 0) 
mat[,1]<-1:n

#deg<-degree(g)                  # Degree based attacks
deg<-betweenness(g)              # Betweenness based attacks
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


#write.csv(Giant_Comp_APL, paste0("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Robustness motifs and others/GC Data/Giant_fr_Italy_GeoDe_Betweeness.csv"))




#################################### Plots #####################################


d0<- read.csv("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Robustness motifs and others/GC Data/Giant_fr_Italy_Original_Betweeness.csv")
head(d0)

d0_GeoDE<- read.csv("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Robustness motifs and others/GC Data/Giant_fr_Italy_GeoDe_Betweeness.csv")
head(d0_GeoDE)


#d0_ER<- read.csv("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Robustness motifs and others/GC Data/Giant_fr_IEEE300ER.csv")
#head(d0_ER)


#d0_PA<- read.csv("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/R codes/2020/Robustness motifs and others/GC Data/Giant_fr_IEEE300PA.csv")
#head(d0_PA)



plot(d0$fr,d0$GC/d0$GC[1], type="l", lty=4, lwd=2, col='black',xlab="Fraction of nodes removed",
     ylab ="Normalized Giant Component" ,ylim=c(),xlim=c(0.0,0.4),first.pannel=grid())

lines(d0_GeoDE$fr,d0_GeoDE$GC/d0_GeoDE$GC[1],  pch = 17,lty=1,lwd=2,col='red') 
#lines(d0_ER$fr,d0_ER$GC/d0_ER$GC[1], type='l', pch = 19, lty=1,lw=2,col='green') 
#lines(d0_PA$fr,d0_PA$GC/d0_PA$GC[1], type='l', pch = 15,lty=1,lwd=2,col='yellow')


legend('topright',c("Observed", "GeoDe"),title="Italian power grid",
       lty=c(4,1,1,1),  lwd=c(2,2),  pch = c(NA,NA),col=c("black", "red")) 










