import networkx as nx
import numpy as np
from scipy import integrate


def backbone(G, alpha):
    ''' Extract backbone of weighted networks
        Args:
            G: NetworkX graph
            alpha: Significance level alpha, e.g., alpha in [0.01, 0.5]
        Return:
            NetworkX graph of backbone network
    '''

    backbone_G = nx.Graph()
    for u in G:
        k = len(G[u])
        if k > 1:
            sum_w = sum(G[u][v]['weight'] for v in G[u])
            for v in G[u]:
                w = G[u][v]['weight']
                p_ij = float(w)/sum_w
                alpha_ij = 1 - (k-1) * integrate.quad(lambda x: (1-x)**(k-2), 0, p_ij)[0]
                if alpha_ij < alpha:
                    backbone_G.add_edge(u, v, weight = w)
    return backbone_G
                
            
if __name__ == '__main__':
    g = nx.barabasi_albert_graph(100, 5)
    for i, j in g.edges():
        g[i][j]['weight'] = np.random.randint(1,100)
    alpha = 0.1
    bbg = backbone(g, alpha)
    print 'alpha = %s' % alpha
    print 'original: nodes = %s, edges = %s' % (g.number_of_nodes(), g.number_of_edges())
    print 'backbone: nodes = %s, edges = %s' % (bbg.number_of_nodes(), bbg.number_of_edges())