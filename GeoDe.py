import scipy
import networkx as nx
from matplotlib import pyplot as plt
import numpy as np
from itertools import combinations
from itertools import product


########################################################################################################

def ER(n,p):
    idx = np.random.geometric(p) - 1
    while idx < n*(n-1)//2:
        yield idx
        idx += np.random.geometric(p)




def DelaunaySphere(n):
    pos = np.random.normal(size = (3,n))
    pos = np.transpose(pos/np.linalg.norm(pos,axis=0))

    #sv = scipy.spatial.SphericalVoronoi(pos,1,[0,0,0])
    #sv.sort_vertices_of_regions()

    # Adapted from scipy SphericalVoronoi
    hull = scipy.spatial.ConvexHull(pos)

    G = nx.Graph()
    G.add_nodes_from(range(n))
    nx.set_node_attributes(G,{ i : pos[i] for i in range(n)}, name = 'position')

    
    triangles = pos[hull.simplices]
    
    for (a,b,c) in hull.simplices:
        G.add_edges_from([(a,b), (b,c),(c,a)])
    return(G)




def random_walk_NoN(n, m = None, p = None):
    if p == None:
        p = 0.80 ### change between 0,1 by, for example, 0.01.
        # 0.50,0.66, 0.75,0.85,0.95
    if m == None:
        m = n

    # This takes an index and maps it to a pair (i,j) where i > j.
    def int2pair(x):
        row = int(np.floor(( 1+ np.sqrt(1+8*x))/2))
        col = x - row*(row-1)//2
        return (row,col)

    G = DelaunaySphere(n)
    DinvA = np.zeros((n,n))
    for (u,v) in G.edges:
        DinvA[u][v] = 1/G.degree[u]
        DinvA[v][u] = 1/G.degree[v]

    P = (m/n)*p*np.matmul(DinvA, np.linalg.inv(np.eye(n) - (1-p)*DinvA))
    P = P + np.transpose(P)
    entries = [P[i][j] for (i,j) in combinations(range(n),2)]
    rho = max(entries)

    H = nx.Graph()
    H.add_nodes_from(G.nodes)
    H.add_edges_from( [(i,j) for (i,j) in map(int2pair,ER(n,rho)) if np.random.random() < P[i][j]/rho])

    return H
        
########################################################################################################
    


"""UnionFind.py

Union-find data structure. Based on Josiah Carlson's code,
http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/215912
with significant additional changes by D. Eppstein.

Modified by Stephen J. Young 5/4/20
-- Remove keys from weights dictionary if they are not roots
-- This allows for easy access to the number of elements len(self.parents) and number of sets len(self.weights)
"""

class UnionFind:
    """Union-find data structure.

    Each unionFind instance X maintains a family of disjoint sets of
    hashable objects, supporting the following two methods:

    - X[item] returns a name for the set containing the given item.
      Each set is named by an arbitrarily-chosen one of its members; as
      long as the set remains unchanged it will keep the same name. If
      the item is not yet part of a set in X, a new singleton set is
      created for it.

    - X.union(item1, item2, ...) merges the sets containing each item
      into a single larger set.  If any item is not yet part of a set
      in X, it is added to X as one of the members of the merged set.
    """

    def __init__(self):
        """Create a new empty union-find structure."""
        self.weights = {}
        self.parents = {}

    def __getitem__(self, object):
        """Find and return the name of the set containing the object."""

        # check for previously unknown object
        if object not in self.parents:
            self.parents[object] = object
            self.weights[object] = 1
            return object

        # find path of objects leading to the root
        path = [object]
        root = self.parents[object]
        while root != path[-1]:
            path.append(root)
            root = self.parents[root]

        # compress the path and return
        for ancestor in path:
            self.parents[ancestor] = root
        return root

            
    def __iter__(self):
        """Iterate through all items ever found or unioned by this structure."""
        return iter(self.parents)

    def union(self, *objects):
        """Find the sets containing the objects and merge them all."""
        roots = [self[x] for x in objects]
        heaviest = max([(self.weights[r],r) for r in roots])[1]
        for r in roots:
            if r != heaviest:
                self.weights[heaviest] += self.weights[r]
                del self.weights[r] # included so the only roots have a weight
                self.parents[r] = heaviest

if __name__ == '__main__':
    import seaborn as sns
    trials =200 # Change nomuber of simutaion here 
    sns.set_palette('hls',16)
    for m in [600]:
#for m in [1000,1100,1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000]:
        data = []
        for t in range(trials):
            #print(m, "\t", t)
            H = random_walk_NoN(465, m = m) #  change m here -------
            #print(H.nodes, H.edges)
            for C in nx.connected_components(H):
                if len(C) < 200:   # change here for subgraph
                    continue
                #print("\t",len(C))
                G = H.subgraph(C)
                print(G)
                
                
                #options = {'node_color': 'black',  'node_size': 10, 'width': 0.5,'edge_color': 'gray'}
                #nx.draw(G)
                #nx.draw_random(G, **options)
                
                dist = dict(nx.shortest_path_length(G))
                nodes = len(C)
                print(nodes)
                
                edges = G.number_of_edges()
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
                            square += len([z for z in G.neighbors(x) if G.has_edge(z,y)])
                    for (x,y,z) in combinations(G.neighbors(v),3):
                        count = sum(1 for (a,b) in combinations([x,y,z],2) if H.has_edge(a,b))
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


                with open('net1.txt','a') as output:
                    output.write( ', '.join(map(repr, [m, t, nodes, edges, vee, triangles, path3, square, k13,
              tent, kite, k4, max(dist[u][v] for (u,v) in combinations(C,2)), np.mean([dist[u][v] for (u,v) in combinations(C,2)])])) + "\n")
                    
                #Diameter=max(dist[u][v] for (u,v) in combinations(C,2))
                #AVPL=np.mean([dist[u][v] for (u,v) in combinations(C,2)])])

###################################################################################################################
                    
                    
