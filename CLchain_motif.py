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


def log_normal_stats(alpha,beta,dmax):
    prob = [ np.exp(-(np.log(d)/alpha)**beta) for d in range(1,dmax+1)]
    prob = [0] + [ p/sum(prob) for p in prob]
    mu = sum( d*prob[d] for d in range(dmax+1))
    sigma = np.sqrt(sum( (d-mu)**2*prob[d] for d in range(dmax+1)))
    return (mu,sigma)

# Given n return the inputs for CLC based on the parameters given in the paper
# Need to have n >= 16
def CLC_defaults(n, dmean = 2.425, dstd = 0.1846):
    # diameter
    c = 1.3101
    k = 0.574
    diameter = int(np.floor(c*(n**k)))
    # max degree
    gamma = 4
    c = 1.517
    dmax = int(np.ceil(c*n**(1/gamma)))
    print(diameter,dmax)
    # parameters alpha, beta
    min_val = (float('inf'),float('inf'))
    min_a = None
    for a in np.linspace(0.001,5,1000):
        #blower = 25/(np.log(np.log(dmax)) - np.log(a))
        #for b in np.linspace(blower/10,10*blower,1000):
        for b in np.linspace(0.001,5,1000):
            (mu,sigma) = log_normal_stats(a,b,dmax)
            if ((mu-dmean)**2, (sigma - dstd)**2) < min_val:
                min_val = ((mu-dmean)**2, (sigma - dstd)**2)
                min_a = (a,b)
    (a,b) = min_a
    prob = [ np.exp(-(np.log(d)/a)**b) for d in range(1,dmax+1)]
    prob = [0] + [ p/sum(prob) for p in prob]
    print(min_a)
    print(min_val)
    print(prob)
    degrees = Counter(rng.choice(dmax+1,size=n,p=prob))
    max_deg = max(degrees.keys())
    degrees = [degrees.get(d,0) for d in range(max_deg+1)]

    return (diameter,degrees)
    

# Source: "A generative graph model for electrical infrastructure networks"
# Authors: Aksoy, Purvine, Cotilla-Sanchez, and Halappanavar
# DOI: 10.1093/comnet/cny016
# Variables:
# degrees is a list containing the number of desired vertices of degree d, i.e. degrees[d] = # vertices of degree d
# diameter is an integer that is the desired diameter of the network
def CLC(degrees, diameter):
    G =nx.Graph()
    
    nu = sum(degrees[1:])
    n = nu
    max_deg = len(degrees)-1
    delta = int(np.round(diameter - 2*np.log(n/(diameter+1))))
    
    # inflate
    while sum(d*degrees[d] - degrees[d]*np.exp(-d) for d in range(1,max_deg+1)) < n:
        j = 1+rng.choice(max_deg,p=[count/n for count in degrees[1:]])
        degrees[j] += 1
        n += 1

    # Select diameter path
    if sum(degrees[3:]) >= delta+1:
        vert_degree = list(rng.permutation(sum([[d]*degrees[d] for d in range(3,max_deg+1)],[]))) + list(rng.permutation(sum([[d]*degrees[d] for d in [1,2]],[])))
        alpha = min(sum(degrees[3:]) - (delta+1),delta+1) - 1
    else:
        vert_degree = list(rng.permutation(sum([[d]*degrees[d] for d in range(2,max_deg+1)],[]))) + [1]*degrees[1]
        alpha = min(sum(degrees[2:]) - (delta+1),delta+1) - 1

    if nu/(delta+1) < max_deg +1:
        block_options = list(rng.choice(delta+1,size = nu//(delta+1), replace = False))
    else:
        block_options = list(range(delta+1))

    # Starting part of subdiagonal
    beta = (delta+1)//2 - (alpha-1)//2

    blocks = list(range(delta+1)) + list(range(beta,beta+alpha+1)) + list(rng.choice(block_options,size = n - (delta+1) - (alpha+1)))

    
    # Start with the diameter and subdiameter path edges
    G.add_edges_from([(i,i+1) for i in range(delta)])
    G.add_edges_from([(delta+1+j,delta+2+j) for j in range(alpha+1)])

    # Build rest of edges

    for i in range(delta+1):
        B = [j for j in range(n) if blocks[j] == i]
        m = sum(vert_degree[j] for j in B)
        endpoints = rng.choice(B,size = m, p = [vert_degree[j]/m for j in B])
        G.add_edges_from([(u,v) for (u,v) in zip(endpoints[::2],endpoints[1::2]) if not u == v])


    giant = max([C for C in nx.connected_components(G)], key = len)
    return G.subgraph(giant)
        
######################################################################################################            

 #import os
 #os.getcwd()
 #os.chdir("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Python/CLC simulated data")
 #os.chdir("C:/Users/asim.kumer/OneDrive/Synthentic Network/Random Graph Model/Python/CLC simulated data")
 #os.getcwd()
 
 
#DD=CLC_defaults(300, dmean = 2.425, dstd = 0.1846)
#DD


#  IEEE 118: D-14, DD=[7 56 19 15 11  6  2  1  1], Edges-178
#  IEEE 300: D-24, DD=[69 76 84 42 14  6  5  2  1  0  1], Edges-409
#  Italy: D-28, DD=[54 84 68 33 21  7  3  3], Node-273, Edges-375
#  Germany: D-31, DD=[121 123 120  34  25  15   3   2   1   1], Node-445, Edges-567


trials=2000

df = pd.DataFrame(columns=['nodes','edges','Diameter','AVPL','triangles','square','k13','tent','kite','k4'])

for t in range(trials):
  print(t)
  G = CLC([121, 123, 120,  34,  25,  15,   3,   2,   1,   1], 31) # Two parameters
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
#https://stackoverflow.com/questions/48314657/how-to-store-values-from-loop-to-a-dataframe

 #import os
 #os.getcwd()
 #os.chdir("C:/Users/asimi/OneDrive/Synthentic Network/Random Graph Model/Python/CLC simulated data")
 #os.chdir("C:/Users/asim.kumer/OneDrive/Synthentic Network/Random Graph Model/Python/CLC simulated data")
 #os.getcwd()
 
 # matrix in a text file
  mat = np.matrix(df)
  with open('Test_Germany_2parar_Motifs.txt','a') as f:
   for line in mat:
    np.savetxt(f, line, fmt='%.2f')
        
        
        
        
        
        
        
        
        
        
        