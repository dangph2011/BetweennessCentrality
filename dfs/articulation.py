import networkx as nx
import matplotlib.pyplot as plt
import random


def main():
    # G = nx.read_edgelist("../data/test2.txt", create_using = nx.DiGraph(), nodetype = int)
    G = nx.read_edgelist("../data/wiki-Vote.txt", nodetype = int)
    print nx.info(G)

    # articulation = nx.articulation_points(G)
    articulation = list(nx.articulation_points(G))
    articulation.sort()

    print articulation
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

if __name__ == '__main__':
	main()
