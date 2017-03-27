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
	std::set<uint32_t> m_articulation;
	std::set<std::pair<uint32_t, uint32_t> > m_bridge;

    std::stack<uint32_t> S;
	bool check = true;
	for (auto i = 0; i < m_graph.size(); i++) {
        if (m_node[i].color == WHITE) {
            uint32_t s = i;
            S.push(s);
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
						m_node[it].color = GRAY;
						m_node[it].discovery_time = ++g_time;
						m_node[it].low = m_node[it].discovery_time;
                        m_node[it].predecessor = u;
						m_node[u].num_child++;
						check = false;
						break;
                    } else if (it != m_node[u].predecessor) {
						m_node[u].low = std::min(m_node[u].low, m_node[it].discovery_time);
					}
                }

				if (check) {
					S.pop();
					if (m_node[u].color == GRAY) {
						m_node[m_node[u].predecessor].low = std::min(m_node[m_node[u].predecessor].low, m_node[u].low);
						m_node[u].color = BLACK;
                    	m_node[u].finish_time = ++g_time;
					}
				}
            }
        }
	}

	//articulation vertex
	for (auto i = 0; i < m_graph.size(); i++){
		if (m_node[i].predecessor == -1) {
			if (m_node[i].num_child > 1) {
				m_articulation.insert(i);
			}
		} else {
			for (auto &it : m_graph[i]) {
				if (m_node[it].low > m_node[i].discovery_time) {
					m_articulation.insert(i);
				}
			}
		}
	}

	// for (auto i = 0; i < m_graph.size(); i++){
	// 	std::cout << m_node[i].low << " " << m_node[i].discovery_time << " " << m_node[i].finish_time << "\n";
	// }

	//Find bridge
	for (auto i = 0; i < m_graph.size(); i++) {
		for (auto &it : m_graph[i]) {
			if (m_node[i].low > m_node[it].discovery_time) {
				m_bridge.insert(std::pair<uint32_t, uint32_t>(it, i));
			}
		}
	}

	//print bridge
	std::cout << "Bridge:\n";
	for (auto &it : m_bridge){
		std::cout << "\t" << it.first << " " << it.second << "\n";
	}

	return 0;
}
