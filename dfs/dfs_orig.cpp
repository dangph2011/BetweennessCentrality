#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

//Travel color, WHITE: unvisited, GRAY: visited but not finish, BLACK: finished
enum eColor {
	WHITE,
	GRAY,
	BLACK
};

//print graph
void print_graph(std::vector<std::vector<unsigned> > graph) {
	for (auto &it : graph){
		for (auto &it1 : it) {
			std::cout << &it - &graph[0] << " " << it1 << std::endl;
		}
	}
}

uint32_t g_time = 0;
//struct Node (Vertex)
typedef struct sNode{
	eColor color;
	uint32_t discovery_time;
	uint32_t finish_time;
	uint32_t predecessor;
	sNode ()
	: color(WHITE), discovery_time(0), finish_time(0) {}
}sNode;

void dfsVisit(std::vector<std::vector<uint32_t> > p_graph, uint32_t u, std::vector<sNode> &p_node) {
	//g_time += 1;
	p_node[u].discovery_time = ++g_time;
	p_node[u].color = GRAY;
	for (auto &it : p_graph[u]) {
		if (p_node[it].color == WHITE) {
			p_node[it].predecessor = u;
			dfsVisit(p_graph, it, p_node);
		}
	}
	p_node[u].color = BLACK;
	//g_time += 1;
	p_node[u].finish_time = ++g_time;
}

int main() {
	std::ifstream m_fin;
 	std::string m_line;
	bool m_directed = true;
	std::vector<std::vector<uint32_t> > m_graph;
	uint32_t u, v;

	std::vector<std::pair<uint32_t, uint32_t> > m_edge_list;
	m_fin.open("../data/wiki-Vote.txt");

	if (m_fin.fail()) {
	    std::cout << "Error: Opening file";
		return 0;
	}

	while (m_fin >> u >> v) {
		m_edge_list.push_back(std::pair<uint32_t, uint32_t>(u,v));
	}

	for (auto &it : m_edge_list) {

		if (m_graph.size() < (it.first + 1)) {
			m_graph.resize(it.first + 1);
		}
		if (m_graph.size() < (it.second + 1)) {
			m_graph.resize(it.second + 1);
		}

		m_graph[it.first].push_back(it.second);
		if (!m_directed) {
			m_graph[it.second].push_back(it.first);
		}
	}

	// //sort adjacency list
	for (auto &it : m_graph) {
		std::sort(it.begin(), it.end());
	}

	//init vector of Node in graph by size
	std::vector<sNode> m_node(m_graph.size());

	for (auto i = 0; i < m_graph.size(); i++) {
		if (m_node[i].color == WHITE) {
			dfsVisit(m_graph, i, m_node);
		}
	}

	for (auto i = 0; i < m_graph.size(); i++){
		std::cout << "Node Id=" << i << " DT=" << m_node[i].discovery_time << " FT=" << m_node[i].finish_time << "\n";
	}

	// for (auto i = 0; i < m_graph.size(); i++){
    //     std::cout << "Node Id=" << i << " Predecessor: ";
    //     for (auto &it : m_node[i].predecessor) {
    //         std::cout << it << " ";
    //     }
	// 	std::cout << std::endl;
    // }

	return 0;
}
