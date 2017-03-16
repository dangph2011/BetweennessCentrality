import networkx as nx
import matplotlib.pyplot as plt
import random


def main():
    G = nx.read_edgelist("../data/test2.txt", create_using = nx.DiGraph(), nodetype = int)
    print nx.info(G)
    betweenness = dict.fromkeys(G, 0.0)  # b[v]=0 for v in G
    nodes = G
    for s in nodes:
        # single source shortest paths
        #if weight is None:  # use BFS
        S, P, sigma = _single_source_shortest_path_basic(G, s)
        #else:  # use Dijkstra's algorithm
            #S, P, sigma = _single_source_dijkstra_path_basic(G, s, weight)
        # accumulation

        betweenness = _accumulate_basic(betweenness, S, P, sigma, s)
    # rescaling

    print betweenness
    #
	#Creates the snapshot of the network

	# spring_pos = nx.spring_layout(G_fb)
	# plt.axis("off")
	# nx.draw_networkx(G_fb, pos = spring_pos, with_labels = False, node_size = 35)
	# plt.savefig("results/fbnetwork.png", format = "PNG")
	# plt.clf()

	# plotBetweeness(G_fb, spring_pos)
	# plt.clf()
	# detectCommunities(G_fb, spring_pos)
	# plt.clf()
	# plt.close()

def _single_source_shortest_path_basic(G, s):
    S = []
    P = {}
    for v in G:
        P[v] = []
    sigma = dict.fromkeys(G, 0.0)    # sigma[v]=0 for v in G
    D = {}
    sigma[s] = 1.0
    D[s] = 0
    Q = [s]
    while Q:   # use BFS to find shortest paths
        v = Q.pop(0)
        S.append(v)
        Dv = D[v]
        sigmav = sigma[v]
        for w in G[v]:
            if w not in D:
                Q.append(w)
                D[w] = Dv + 1
            if D[w] == Dv + 1:   # this is a shortest path, count paths
                sigma[w] += sigmav
                P[w].append(v)  # predecessors
    return S, P, sigma

def _accumulate_basic(betweenness, S, P, sigma, s):
    delta = dict.fromkeys(S, 0)
    while S:
        w = S.pop()
        coeff = (1.0 + delta[w]) / sigma[w]
        for v in P[w]:
            delta[v] += sigma[v] * coeff
        if w != s:
            betweenness[w] += delta[w]
    return betweenness

if __name__ == '__main__':
	main()
