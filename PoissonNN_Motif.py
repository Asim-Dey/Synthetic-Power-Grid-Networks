import networkx as nx
import numpy as np
from collections import Counter
import scipy
import networkx as nx
from matplotlib import pyplot as plt
import numpy as np
from itertools import combinations
from itertools import product
import pandas as pd


rng = np.random.default_rng()


# Source: Network Models for Power Grids: A Generative Approach
# Authors: Deepjyoti Deka and Sriram Vishwanath
# DOI: ???
# Variables:
# degrees is a list containing the number of desired vertices of degree d, i.e. degrees[d] = # vertices of degree d
def PoissonNN(degrees):
    n = sum(degrees[1:])
    G =nx.Graph()

    radius = np.sqrt(rng.random(size=n))
    angle = 2*np.pi*rng.random(size=n)
    G.add_nodes_from(range(n))
    dist = [d/n for d in degrees]
    entry_deg = rng.choice(len(degrees),size =n, p = dist)

    # adujust degrees to account for future degrees added
    prev_deg = entry_deg[0]
    for i in range(1,n):
        entry_deg[i] = max(1,int(np.round(entry_deg[i] - (prev_deg/i)*np.log(n/(n-i)))))
        prev_deg += entry_deg[i]

    entry_deg = list(entry_deg)
    entry_deg.reverse()
    print(entry_deg)
    for i in range(1,n):
        print(i)
        if entry_deg[i] >= i:
            G.add_edges_from([(j,i) for j in range(i)])
            continue
        distances = [ (radius[i]**2 + radius[j]**2 - 2*radius[i]*radius[j]*np.cos(min(np.abs(angle[i]-angle[j]),2*np.pi - np.abs(angle[i]-angle[j]))),j) for j in range(i)]
        distances.sort()
        G.add_edges_from([(distances[j][1],i) for j in range(entry_deg[i])])

    return G
     
 ##############################################################################
 
 

# include a zero in the beginning
G1=PoissonNN([0,121, 123, 120,  34,  25,  15, 3, 2, 1, 1])



trials=2000

df = pd.DataFrame(columns=['nodes','edges','Diameter','AVPL','triangles','square','k13','tent','kite','k4'])

for t in range(trials):
  print(t)
  G=PoissonNN([0,121, 123, 120, 34, 25, 15, 3, 2, 1, 1]) # degree distribution --change here
  #nodes=G.number_of_nodes()
  df.loc[t, ['nodes']] = G.number_of_nodes()

  #edges = G.number_of_edges()
  df.loc[t, ['edges']] = G.number_of_edges()
  
  dist = dict(nx.shortest_path_length(G))
  Diameter1=max(dist[u][v] for (u,v) in combinations(G,2))
  df.loc[t, ['Diameter']] = max(dist[u][v] for (u,v) in combinations(G,2))
  
  #AVPL1=np.mean([dist[u][v] for (u,v) in combinations(G,2)])
  df.loc[t, ['AVPL']] = np.mean([dist[u][v] for (u,v) in combinations(G,2)])


  triangles = 0
  vee = 0
  path3 = 0
  square = 0
  k13 = 0
  tent = 0
  kite = 0
  k4 = 0

  for v in G.nodes:
      for (x,y) in combinations(G.neighbors(v),2):
          if G.has_edge(x,y):
              triangles += 1
          else:
              vee += 1
              square += sum(1 for z in G.neighbors(x) if ((not z == v) and (G.has_edge(z,y)) and (not G.has_edge(z,v))))#len([z for z in G.neighbors(x) if G.has_edge(z,y)])
      for (x,y,z) in combinations(G.neighbors(v),3):
          count = sum(1 for (a,b) in combinations([x,y,z],2) if G.has_edge(a,b))
          if count == 0:
              k13 += 1
          elif count == 1:
              tent += 1
          elif count == 2:
              kite += 1
          else:
              k4 += 1
  for (u,v) in G.edges:
      path3 += sum([1 for (z,y) in product(G.neighbors(u),G.neighbors(v)) if not any([G.has_edge(z,v), G.has_edge(z,y), G.has_edge(u,y)])])
  
  triangles = triangles/3
  square = square/4
  kite = kite/2
  k4 = k4/4
    
  df.loc[t, ['triangles']] = triangles
  df.loc[t, ['square']] = square
  df.loc[t, ['k13']] = k13
  df.loc[t, ['tent']] = tent
  df.loc[t, ['kite']] = kite
  df.loc[t, ['k4']] = k4


  print(df)
  np.shape(df) ### Dimension 

## T2-Triangle; V3/M3-Tent; V4/M4-Square; V5/M5-kite 

######### saving in a matrix ###########

 # matrix in a text file
  mat = np.matrix(df)
  with open('Test_Germany_eNN_Motifs.txt','a') as f:
   for line in mat:
    np.savetxt(f, line, fmt='%.2f')
        




 











