#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <stack>

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
	std::vector<uint32_t> predecessor;
	sNode ()
	: color(WHITE), discovery_time(0), finish_time(0) {}
}sNode;

int main() {
	std::ifstream m_fin;
 	std::string m_line;
	bool m_directed = true;
	std::vector<std::vector<uint32_t> > m_graph;
	uint32_t u, v;

	std::vector<std::pair<uint32_t, uint32_t> > m_edge_list;
	m_fin.open("../data/test1.txt");

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

    std::stack<uint32_t> S;
	bool check = true;
	for (auto i = 0; i < m_graph.size(); i++) {
        if (m_node[i].color == WHITE) {
            uint32_t s = i;
            S.push(s);
            //std::cerr << "test\n";
            while (!S.empty()) {
                uint32_t u = S.top();
                //S.pop();
                //std::cerr << "test=" << g_time << " " << u << " " << S.size() << "\n";
                if (m_node[u].color == WHITE) {
                    m_node[u].color = GRAY;
                    m_node[u].discovery_time = ++g_time;
                }
				
				check = true;
                for (auto &it : m_graph[u]) {
                    if (m_node[it].color == WHITE) {
                        //m_node[it].discovery_time = ++g_time;
                        S.push(it);
                        m_node[it].predecessor.push_back(u);
						check = false;
                    }
                }

				if (check) {
					S.pop();
					if (m_node[u].color == GRAY) {
						m_node[u].color = BLACK;
                    	m_node[u].finish_time = ++g_time;
					}
				}
            }

            // m_node[s].color = BLACK;
            // m_node[s].discovery_time = ++g_time;
        }
	}

	for (auto i = 0; i < m_graph.size(); i++){
		std::cout << "Node Id=" << i << " DT=" << m_node[i].discovery_time << " FT=" << m_node[i].finish_time << "\n";
	}

    for (auto i = 0; i < m_graph.size(); i++){
        std::cout << "Node Id=" << i << " Predecessor: ";
        for (auto &it : m_node[i].predecessor) {
            std::cout << it << " ";
        }
		std::cout << std::endl;
    }

	return 0;
}
