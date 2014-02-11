'''
This module implements the disparity filter to compute a significance score of edge weights in networks
'''

import networkx as nx
import numpy as np
from scipy import integrate


def disparity_filter(G):
    ''' Compute significance scores (alpha) for weighted edges in G as defined in Serrano et al. 2009
        Args
            G: Weighted NetworkX graph
        Returns
            Weighted graph with a significance score (alpha) assigned to each edge
        References
            M. A. Serrano et al. (2009) Extracting the Multiscale backbone of complex weighted networks. PNAS, 106:16, pp. 6483-6488.
    '''
    B = nx.Graph()
    for u in G:
        k = len(G[u])
        if k > 1:
            sum_w = sum(G[u][v]['weight'] for v in G[u])
            for v in G[u]:
                w = G[u][v]['weight']
                p_ij = float(w)/sum_w
                alpha_ij = 1 - (k-1) * integrate.quad(lambda x: (1-x)**(k-2), 0, p_ij)[0]
                B.add_edge(u, v, weight=w, alpha=float('%.4f' % alpha_ij))
    return B
                
            
if __name__ == '__main__':
    G = nx.barabasi_albert_graph(1000, 5)
    for u, v in G.edges():
        G[u][v]['weight'] = np.random.randint(1,100)
    alpha = 0.05
    G = disparity_filter(G)
    G2 = nx.Graph([(u, v, d) for u, v, d in G.edges(data=True) if d['alpha'] < alpha])
    print 'alpha = %s' % alpha
    print 'original: nodes = %s, edges = %s' % (G.number_of_nodes(), G.number_of_edges())
    print 'backbone: nodes = %s, edges = %s' % (G2.number_of_nodes(), G2.number_of_edges())
    print G2.edges(data=True)