import networkx as nx
import numpy as np

rng = np.random.default_rng()


# Source: Network Models for Power Grids: A Generative Approach
# Authors: Deepjyoti Deka and Sriram Vishwanath
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
 
 
G1=PoissonNN([0,121, 123, 120,  34,  25,  15, 3, 2, 1, 1])

nodes=G1.number_of_nodes()
edges = G1.number_of_edges()
edges
nodes
#G.degree()''

options = {'node_color': 'black',  'node_size': 10, 'width': 0.5,'edge_color': 'gray'}
nx.draw(G1, **options)

####### Edgelist ###################


nx.generate_edgelist(G1, data=False)
nx.write_edgelist(G1, "Germany_NN_10.txt",data=False)


 











