#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <stack>
#include <set>
#include <unistd.h>

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
	int32_t predecessor;
	uint32_t num_child;
	uint32_t low;
	sNode ()
	: color(WHITE), discovery_time(0), finish_time(0), predecessor(-1), num_child(0), low(0) {}
}sNode;

int main() {
	std::ifstream m_fin;
 	std::string m_line;
	bool m_directed = false;
	std::vector<std::vector<uint32_t> > m_graph;
	uint32_t u, v;

	std::vector<std::pair<uint32_t, uint32_t> > m_edge_list;
	m_fin.open("../data/test2.txt");

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
	std::set<uint32_t> m_articulation;
	std::set<std::pair<uint32_t, uint32_t> > m_bridge;
	std::stack<std::pair<uint32_t, uint32_t> > m_list_edge;

    std::stack<uint32_t> S;
	std::stack<uint32_t> m_node_component;
	bool check = true;
	for (auto i = 0; i < m_graph.size(); i++) {
        if (m_node[i].color == WHITE) {
            uint32_t s = i;
            S.push(s);
			m_node_component.push(s);
            //std::cerr << "test\n";
			m_node[s].color = GRAY;
			m_node[s].discovery_time = ++g_time;
			m_node[s].low = m_node[s].discovery_time;

            while (!S.empty()) {
                uint32_t u = S.top();
                //S.pop();

				check = true;
                for (auto &it : m_graph[u]) {
                    if (m_node[it].color == WHITE) {
                        //m_node[it].discovery_time = ++g_time;
                        S.push(it);
						m_node_component.push(it);
						m_node[it].color = GRAY;
						m_node[it].discovery_time = ++g_time;
						m_node[it].low = m_node[it].discovery_time;
                        m_node[it].predecessor = u;
						m_node[u].num_child++;
						m_list_edge.push(std::pair<uint32_t, uint32_t>(u, it));
						//std::cout << "edges0: " << u << " " << it << std::endl;
						check = false;
						break;
                    } else if (it != m_node[u].predecessor && u != m_node[it].predecessor) {
						//TODO check if u it is travelled
						m_node[u].low = std::min(m_node[u].low, m_node[it].discovery_time);
						m_list_edge.push(std::pair<uint32_t, uint32_t>(u, it));
						//std::cout << "edges1: " << u << " " << it << std::endl;
					}
                }

				if (check) {
					if (m_node[u].color == GRAY) {
						m_node[u].color = BLACK;
                    	m_node[u].finish_time = ++g_time;
						if (m_node[u].predecessor > -1) {

							auto l_predecessor = m_node[u].predecessor;
							m_node[l_predecessor].low = std::min(m_node[l_predecessor].low, m_node[u].low);

							// std::cout << "edges1: " << u << " " << l_predecessor << std::endl;
							// m_list_edge.push(std::pair<uint32_t, uint32_t>(u, l_predecessor));

							if (l_predecessor > -1 && m_node[u].low > m_node[l_predecessor].discovery_time) {
								m_bridge.insert(std::pair<uint32_t, uint32_t>(l_predecessor, u));
								std::cout << "bridges: " << l_predecessor << "--" << u << " " << std::endl;

								auto l_pop_edge = m_list_edge.top();
								while (l_pop_edge.first != l_predecessor || l_pop_edge.second != u) {
									std::cout << l_pop_edge.first << "--" << l_pop_edge.second << " ";
									m_list_edge.pop();
									l_pop_edge = m_list_edge.top();
								}
								//std::cout << l_pop_edge.first << "--" << l_pop_edge.second << " ";
								m_list_edge.pop();
								std::cout << std::endl;

								//std::cout << "Node Component\n";
								auto l_pop_node = m_node_component.top();
								while (l_pop_node != l_predecessor && l_pop_node != u) {
									std::cout << l_pop_node << " ";
									m_node_component.pop();
									l_pop_node = m_node_component.top();
								}
								std::cout << l_pop_node << " ";
								m_node_component.pop();
								std::cout << std::endl;
							}
						}
					}
					S.pop();
				}
            }
        }

		//std::cout << "i=" << i << "\n";
		while (!m_list_edge.empty()) {
			auto l_pop_edge = m_list_edge.top();
			std::cout << l_pop_edge.first << "--" << l_pop_edge.second << " ";
			m_list_edge.pop();
		}
		std::cout << std::endl;

		//std::cout << "Node Component\n";
		//auto l_pop_node = m_node_component.top();
		while (!m_node_component.empty()) {
			auto l_pop_node = m_node_component.top();
			std::cout << l_pop_node << " ";
			m_node_component.pop();
		}
		std::cout << std::endl;
	}

	//articulation vertex
	// for (auto i = 0; i < m_graph.size(); i++){
	// 	if (m_node[i].predecessor == -1) {
	// 		if (m_node[i].num_child > 1) {
	// 			m_articulation.insert(i);
	// 		}
	// 	} else {
	// 		//std::cout << i << " " << m_node[i].predecessor << " " << m_node[i].low << " " << m_node[m_node[i].predecessor].discovery_time << std::endl;
	// 		if (m_node[m_node[i].predecessor].predecessor > -1 && m_node[i].low >= m_node[m_node[i].predecessor].discovery_time) {
	// 			m_articulation.insert(m_node[i].predecessor);
	// 		}
	// 	}
	// }

	for (auto i = 0; i < m_graph.size(); i++){
		std::cout << m_node[i].low << " " << m_node[i].discovery_time << " " << m_node[i].finish_time << "\n";
	}

	//Find bridge
	// for (auto i = 0; i < m_graph.size(); i++) {
	// 	if (m_node[i].predecessor == -1)
	// 		continue;
	// 	if (m_node[i].low > m_node[m_node[i].predecessor].discovery_time) {
	// 		m_bridge.insert(std::pair<uint32_t, uint32_t>(m_node[i].predecessor, i));
	// 	}
	// }

	//print bridge
	std::cout << "Bridge:\n";
	for (auto &it : m_bridge){
		std::cout << "\t" << it.first << " " << it.second << "\n";
	}

	return 0;
}
