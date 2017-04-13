#include <iostream>
#include <vector>
#include <set>
#include <fstream>
#include <stack>
#include <queue>
#include <iomanip>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <unistd.h>


std::vector<double> g_bet_cen;
//Travel color, WHITE: unvisited, GRAY: visited but not finish, BLACK: finished
enum eColor {
	WHITE,
	GRAY,
	BLACK
};

/**
    Get current time in milisecond.
    @return the current time in milisecond
*/
double getCurrentTimeMlsec() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	//second
	return tv.tv_sec + tv.tv_usec * 1e-6;
	//milisecond
	//return (tv.tv_sec + tv.tv_usec / 1000000.0) * 1000;
}

//set bit at index in the array addr
static inline void set_bit(uint32_t index, uint32_t *addr)
{
	uint32_t p = index >> 5, x = index & 31; //& ~ %31
	addr[p] |= 1 << x;
}

//test bit at index in the array addr
static inline int test_bit(uint32_t index, const uint32_t *addr)
{
	uint32_t p = index >> 5, x = index & 31;
	return (addr[p] >> x) & 1;
}

//struct Node (Vertex)
typedef struct Vertex{
	eColor color;
	uint32_t discovery_time;
	uint32_t finish_time;
	int32_t predecessor;
	uint32_t num_child;
	uint32_t low;
	Vertex ()
	: color(WHITE), discovery_time(0), finish_time(0), predecessor(-1), num_child(0), low(0) {}
}Vertex;

typedef struct GraphConnected {
	int32_t graph_id;
	//Edge from u to v
	//From graph
	int32_t u_start;
	//to graph
	int32_t v_end;

	bool traveled;
	GraphConnected()
	: graph_id(-1), u_start(-1), v_end(-1), traveled(false) {}

	GraphConnected(int32_t p_graph_id, int32_t p_u_start, int32_t p_v_end) {
		graph_id = p_graph_id;
		u_start = p_u_start;
		v_end = p_v_end;
		traveled = false;
	}

}GraphConnected;

class GraphComponent{
	public:
	std::vector<uint32_t> graph_id_;
	std::vector<std::vector<uint32_t> > edge_list_;
	bool directed_;
	GraphComponent(bool p_directed = false) {
		directed_ = p_directed;
	}

	void addEdge(uint32_t u, uint32_t v) {
		if (edge_list_.size() < (u + 1)) {
			edge_list_.resize(u + 1);
		}
		if (edge_list_.size() < (v + 1)) {
			edge_list_.resize(v + 1);
		}

		edge_list_[u].push_back(v);
		if (!directed_) {
			edge_list_[v].push_back(u);
		}
	}

	void delEdge(uint32_t u, uint32_t v) {
		for (auto &it : edge_list_[u]) {
			if (it == v) {
				edge_list_[u].erase(edge_list_[u].begin() + (&it - &edge_list_[u][0]));
				break;
			}
		}

		if (!directed_) {
			for (auto &it : edge_list_[v]) {
				if (it == u) {
					edge_list_[v].erase(edge_list_[v].begin() + (&it - &edge_list_[v][0]));
					break;
				}
			}
		}
	}

	//print graph
    void printEdges() {
		//
		std::cout << "List Edges\n";
    	for (auto &it : edge_list_){
    		for (auto &it1 : it) {
    			std::cout << &it - &edge_list_[0] << " " << it1 << std::endl;
    		}
    	}
    }
};

class Graph{
	public:
	int32_t graph_id_;
	//int32_t graph_id_connect_;
	uint32_t v_num_;
	uint32_t e_num_;
    std::vector<std::vector<uint32_t> > edge_list_;
    std::set<uint32_t> vertex_list_;
	std::vector<uint32_t> reach_;
	//relatively predecessor/parent
	GraphConnected graph_connected_;
	//int32_t v_connect_;

    bool directed_ = false;

	Graph(){
		graph_id_ = -1;
		//v_connect_ = -1;
		//graph_id_connect_ = -1;
	}

	Graph(uint32_t p_size) {
		edge_list_.resize(p_size);
		reach_.resize(p_size, 1);
		//initReach();
		graph_id_ = -1;
		//v_connect_ = -1;
		//graph_id_connect_ = -1;
    }

    Graph(uint32_t p_size, uint32_t p_graph_id) {
		edge_list_.resize(p_size);
		reach_.resize(p_size, 1);
		graph_id_ = p_graph_id;

		//v_connect_ = -1;
		//graph_id_connect_ = -1;
    }

	void initReach() {
		for (auto &it : reach_) {
			it = 1;
		}
	}

	void setGraphConnected(uint32_t p_graph_id, uint32_t pu_start, uint32_t pv_end) {
		graph_connected_.graph_id = p_graph_id;
		graph_connected_.u_start = pu_start;
		graph_connected_.v_end = pv_end;
	}

	void setGraphConnectedId(uint32_t p_graph_id) {
		graph_connected_.graph_id = p_graph_id;
	}

	uint32_t getTotalReach() {
		uint32_t total = 0;
		for (auto &it : vertex_list_) {
			total += reach_[it];
		}
		return total;
	}

    uint32_t getNumberOfVectex() {
        return vertex_list_.size();
    }

    uint32_t getNumberOfEdge() {
        return edge_list_.size();
    }

	// int32_t getVConnect() {
	// 	return v_connect_;
	// }
	//
	// void setVConnect(int32_t p_vertex_id) {
	// 	v_connect_ = p_vertex_id;
	// }

	bool getDirected() {
		return directed_;
	}

	void setDirected(bool p_directed) {
		directed_ = p_directed;
	}

	int32_t getGraphId(){
		return graph_id_;
	}

	void setGraphId(int32_t p_id) {
		graph_id_ = p_id;
	}

	void setReach(uint32_t p_vertex, uint32_t p_value) {
		if (reach_.size() < p_vertex) {
			reach_.resize(p_vertex+1, 1);
		}
		reach_[p_vertex] = p_value;
	}

	uint32_t getReach(uint32_t p_vertex) {
		return reach_[p_vertex];
	}

	// void initAndSetReach() {
	// 	reach_.resize(edge_list_.size());
	// 	for (auto &it : reach_) {
	// 		it = 1;
	// 	}
	// }

	void addVertex(uint32_t u) {
		if (edge_list_.size() < (u + 1)) {
			edge_list_.resize(u + 1);
		}

		if (reach_.size() < (u + 1)){
			reach_.resize(u+1, 1);
		}

		vertex_list_.insert(u);
	}

	void addEdge(uint32_t u, uint32_t v) {
		if (edge_list_.size() < (u + 1)) {
			edge_list_.resize(u + 1);
			reach_.resize(u + 1, 1);
		}
		if (edge_list_.size() < (v + 1)) {
			edge_list_.resize(v + 1);
			reach_.resize(v + 1, 1);
		}

		edge_list_[u].push_back(v);
		//reach_[u] = 1;
		//reach_[v] = 1;
		if (!directed_) {
			edge_list_[v].push_back(u);
		}
	}

	void delEdge(uint32_t u, uint32_t v) {
		for (auto &it : edge_list_[u]) {
			if (it == v) {
				edge_list_[u].erase(edge_list_[u].begin() + (&it - &edge_list_[u][0]));
				break;
			}
		}

		if (!directed_) {
			for (auto &it : edge_list_[v]) {
				if (it == u) {
					edge_list_[v].erase(edge_list_[v].begin() + (&it - &edge_list_[v][0]));
					break;
				}
			}
		}
	}

	void sortEdges() {
		// //sort adjacency list
		for (auto &it : edge_list_) {
			std::sort(it.begin(), it.end());
		}
	}

	void sortAndRemoveDuplicateEdges() {
		// //sort adjacency list
		for (auto &it : edge_list_) {
			std::sort(it.begin(), it.end());
			it.erase( unique( it.begin(), it.end() ), it.end() );
		}
	}

	//print Vertex, Edge and graph Id connected
	void printBridgeInfo(){
		std::cout << "GraphID=" << graph_connected_.graph_id << " Start=" << graph_connected_.u_start << " End=" << graph_connected_.v_end << std::endl;
	}

    //print graph
    void printGraph() {
		//
		printVertices();
    	printEdges();
    }

	//print graph
    void printEdges() {
		//
		std::cout << "List Edges\n";
    	for (auto &it : edge_list_){
    		for (auto &it1 : it) {
    			std::cout << &it - &edge_list_[0] << " " << it1 << std::endl;
    		}
    	}
    }

	//print graph
    void printVertices() {
		//
		std::cout << "List Vertices\n";
    	for (auto &it : vertex_list_){
    		std::cout << it << " " << reach_[it] << std::endl;
    	}
    }
};

bool readEdgeList(Graph &g, std::string p_file_name, bool p_directed = false){
	g.directed_ = p_directed;
	std::ifstream f_in;
	f_in.open(p_file_name);
	if (f_in.fail()) {
		std::cerr << "Error: Opening file\n";
		return false;
	}
	uint32_t u, v;

	std::vector<std::pair<uint32_t, uint32_t> > edge_list;

	while (f_in >> u >> v){
		edge_list.push_back(std::pair<unsigned, unsigned>(u,v));
	}

	//Store in adjacency list
	for (auto &it : edge_list) {
		g.vertex_list_.insert(it.first);
		g.vertex_list_.insert(it.second);

		//std::cout << graph.size();
		if (g.edge_list_.size() < (it.first + 1)) {
			g.edge_list_.resize(it.first + 1);
			g.reach_.resize(it.first + 1, 1);
		}
		if (g.edge_list_.size() < (it.second + 1)) {
			g.edge_list_.resize(it.second + 1);
			g.reach_.resize(it.second + 1, 1);
		}

		g.edge_list_[it.first].push_back(it.second);
		//g.reach_[it.first] = 1;
		//g.reach_[it.second] = 1;
		if (!g.directed_) {
			g.edge_list_[it.second].push_back(it.first);
		}
	}

	//g.initAndSetReach();
	// //sort adjacency list
	g.sortAndRemoveDuplicateEdges();
	// for (auto &it : g.edge_list_) {
	// 	std::sort(it.begin(), it.end());
	// }
	return true;
}

bool readMetis(Graph &g, std::string p_file_name, bool p_directed = false) {
	g.setDirected(p_directed);
	std::ifstream f_in;
	std::string line;
	uint32_t v;
	uint32_t count = 1;
	f_in.open(p_file_name);
	if (f_in.fail()) {
		std::cerr << "Error: Opening file\n";
		return false;
	}
	getline(f_in, line);
	std::istringstream iss1(line);
	iss1 >> g.v_num_ >> g.e_num_;
	//std::cout << "Vertex=" << v_num_ << " Edge=" << e_num_ << std::endl;
	g.edge_list_.resize(g.v_num_+1);
	g.reach_.resize(g.v_num_+1, 1);
	//std::cout << "Vertex=" << v_num_ << " Edge=" << e_num_ << std::endl;
	while (getline(f_in, line)) {
		std::istringstream iss(line);
		//std::cout << "Line=" << line << std::endl;
		if (!line.empty()) {
			g.vertex_list_.insert(count);
			//g.reach_[count] = 1;
		}

		while (iss >> v) {
			//std::cout << "\tv=" << v << " " << std::endl;
			//g.vertex_list_.insert(count);
			g.vertex_list_.insert(v);
			g.edge_list_[count].push_back(v);
			//g.reach_[v] = 1;
			if (!g.directed_) {
				g.edge_list_[v].push_back(count);
			}
		}

		count++;
		//std::cout << count << "\n";
	}

	//g.initAndSetReach();
	//sort and erase duplicate
	g.sortAndRemoveDuplicateEdges();
	// for (auto &it : g.edge_list_) {
	// 	std::sort(it.begin(), it.end());
	// 	it.erase( unique( it.begin(), it.end() ), it.end() );
	// }

	//std::cout << "END\n";
	return true;
}

void betweennessCentrality(Graph &g, std::vector<double> &l_bet_cen) {
	//l_bet_cen.resize(g.edge_list_.size());

	//set all the value of betweenness centrality of all vectices equal to 0 before computing
	// for (auto &it : l_bet_cen){
	// 	it = 0;
	// }

	//travel all vectice
	for (auto &i : g.vertex_list_) {
		//for (size_t j = 0; j < vertex_list_[i].size(); j++) {
			uint32_t s = i;
			//std::cerr << vertex_list_.size() << std::endl;
			//sleep(1);

			uint32_t *check_map = new uint32_t [g.edge_list_.size()/32 + 1];
			std::memset(check_map, 0, (g.edge_list_.size()/32+1)*sizeof(uint32_t));

			std::stack<uint32_t> S; //visited vertices
			std::vector<std::vector<uint32_t> > predecessor(g.edge_list_.size());
			std::queue<uint32_t> Q; //unvisited vertices
			std::vector<int32_t> number_of_shortest_path(g.edge_list_.size());
			for (auto &it : number_of_shortest_path){
				it = 0;
			}

			number_of_shortest_path[s] = 1;

			std::vector<int64_t> distance(g.edge_list_.size());
			for (auto &it : distance){
				it = -1;
			}

			distance[s] = 0;
			Q.push(s);
			set_bit(s, check_map);
			while (!Q.empty()) {
				uint32_t v = Q.front();
				Q.pop();
				S.push(v);
				//find all neighbor w of v
				for (size_t k = 0; k < g.edge_list_[v].size(); k++) {
					uint32_t w = g.edge_list_[v][k];
					//check if w is visited vertex

					if (!test_bit(w, check_map)) {
						set_bit(w, check_map);	//set w is visited vertex
						if (distance[w] < 0) {
							Q.push(w);
							distance[w] = distance[v] + 1;
						}
					}

					if (distance[w] == distance[v] + 1) {
						number_of_shortest_path[w] += number_of_shortest_path[v];
						predecessor[w].push_back(v);
					}
				}
			}

			std::vector<double> dependency(g.edge_list_.size());
			for (auto j = 0; j < dependency.size(); j++) {
				dependency[j] = g.getReach(j) - 1;
			}

			//S returns vertices in order of non-increasing distance from s
			while (!S.empty()) {
				uint32_t w = S.top();
				S.pop();
				for (auto &it : predecessor[w]) {
					int v = it;
					dependency[v] += ((double)number_of_shortest_path[v] / (double)number_of_shortest_path[w]) * (1 + (double)dependency[w]);
				}

				if (w != s) {
					l_bet_cen[w] += dependency[w] * g.getReach(s);
				}
			}
			delete [] check_map;
	}

	// divide by 2 for undirected graph
	// if (!g.directed_) {
	// 	for (auto &it : l_bet_cen) {
	// 		it /= 2;
	// 	}
	// }
	//return l_bet_cen;
}

//DFS travel
bool dfsTravelBridge(Graph &g, std::set<std::pair<uint32_t, uint32_t> > &p_bridge, 	std::vector<Graph> &p_graph_shatter) {
	uint32_t g_time = 0;
	//init vector of Node in graph by size
	std::vector<Vertex> l_node(g.edge_list_.size());
	std::stack<std::pair<uint32_t, uint32_t> > m_list_edge;
	std::stack<uint32_t> m_node_component;
	uint32_t l_graph_shatter_size_begin = p_graph_shatter.size();
	uint32_t l_graph_shatter_size = p_graph_shatter.size();

	std::vector<std::vector<uint32_t> > l_graph_id(g.edge_list_.size());

	std::stack<uint32_t> S;
	//uint32_t m_num_tree_node;
	bool check = true;

	for (auto &i : g.vertex_list_) {
		//l_graph_shatter_size_start = p_graph_shatter.size();
		//l_graph_shatter_size_start = p_graph_shatter.size();
		//m_num_tree_node = g_time;
		//std::vector<Graph> l_graph;
		//GraphComponent l_graph_component;
		if (l_node[i].color == WHITE) {
			//m_num_tree_node++;
			uint32_t s = i;
			S.push(s);
			m_node_component.push(s);
			l_node[s].color = GRAY;
			l_node[s].discovery_time = ++g_time;
			l_node[s].low = l_node[s].discovery_time;

			while (!S.empty()) {
				uint32_t u = S.top();
				//S.pop();

				check = true;
				for (auto &it : g.edge_list_[u]) {
					if (l_node[it].color == WHITE) {
						//l_node[it].discovery_time = ++g_time;
						S.push(it);
						m_node_component.push(it);
						l_node[it].color = GRAY;
						l_node[it].discovery_time = ++g_time;
						l_node[it].low = l_node[it].discovery_time;
						l_node[it].predecessor = u;
						l_node[u].num_child++;
						m_list_edge.push(std::pair<uint32_t, uint32_t>(u, it));
						check = false;
						break;
					} else if (it != l_node[u].predecessor && u != l_node[it].predecessor) {
						l_node[u].low = std::min(l_node[u].low, l_node[it].discovery_time);
						m_list_edge.push(std::pair<uint32_t, uint32_t>(u, it));
					}
				}

				if (check) {
					if (l_node[u].color == GRAY) {
						l_node[u].color = BLACK;
						l_node[u].finish_time = ++g_time;

						if (l_node[u].predecessor > -1) {
							auto l_predecessor = l_node[u].predecessor;
							l_node[l_predecessor].low = std::min(l_node[l_predecessor].low, l_node[u].low);

							if (l_predecessor > -1 && l_node[u].low > l_node[l_predecessor].discovery_time) {
								p_bridge.insert(std::pair<uint32_t, uint32_t>(l_predecessor, u));

								//Init graph with current version
								Graph l_g(g.edge_list_.size(), l_graph_shatter_size);

								//add edges into graph
								auto l_pop_edge = m_list_edge.top();
								//check if edge pop not bridge
								while (l_pop_edge.first != l_predecessor || l_pop_edge.second != u) {
									//add edge to current graph
									l_g.addEdge(l_pop_edge.first, l_pop_edge.second);
									m_list_edge.pop();
									l_pop_edge = m_list_edge.top();
								}

								m_list_edge.pop();
								//add note into graph
								auto l_pop_node = m_node_component.top();
								while (l_pop_node != u) {
									l_g.addVertex(l_pop_node);
									l_g.setReach(l_pop_node, g.getReach(l_pop_node));
									//find if node belong one bridge
									if (!l_graph_id[l_pop_node].empty()) {
										for (auto &it_li : l_graph_id[l_pop_node]) {
											p_graph_shatter[it_li].setGraphConnectedId(l_graph_shatter_size);
										}
									}
									m_node_component.pop();
									l_pop_node = m_node_component.top();
								}
								//std::cout << l_pop_node << " ";
								l_graph_id[l_predecessor].push_back(l_graph_shatter_size);
								l_g.addVertex(l_pop_node);
								l_g.setReach(l_pop_node, g.getReach(l_pop_node));
								if (!l_graph_id[l_pop_node].empty()) {
									for (auto &it_li : l_graph_id[l_pop_node]) {
										p_graph_shatter[it_li].setGraphConnectedId(l_graph_shatter_size);
									}
								}
								//have yet what graph to connect
								//GraphConnected t_graph_connect(-1, l_pop_node, l_predecessor);
								l_g.setGraphConnected(-1, l_pop_node, l_predecessor);

								m_node_component.pop();
								l_g.sortAndRemoveDuplicateEdges();
								p_graph_shatter.push_back(l_g);
								//Increase version
								l_graph_shatter_size++;
							}
						}
					}
					S.pop();
				}
			}
		}

		if (!m_node_component.empty() || !m_list_edge.empty()) {
			Graph l_g(g.edge_list_.size(),l_graph_shatter_size);

			while (!m_list_edge.empty()) {
				auto l_pop_edge = m_list_edge.top();
				l_g.addEdge(l_pop_edge.first, l_pop_edge.second);
				m_list_edge.pop();
			}

			while (!m_node_component.empty()) {
				auto l_pop_node = m_node_component.top();
				l_g.addVertex(l_pop_node);
				l_g.setReach(l_pop_node, g.getReach(l_pop_node));
				if (!l_graph_id[l_pop_node].empty()) {
					for (auto &it_li : l_graph_id[l_pop_node]) {
						//std::cout << "1-LgraphId=" << it_li << " " << l_graph_shatter_size << std::endl;
						p_graph_shatter[it_li].setGraphConnectedId(l_graph_shatter_size);
					}
				}
				m_node_component.pop();
			}

			l_g.sortAndRemoveDuplicateEdges();
			p_graph_shatter.push_back(l_g);
			l_graph_shatter_size++;
		}

		//m_num_tree_node = g_time - m_num_tree_node;

		//TODO: Update reach for each graph shattered
		uint32_t f_total_reach_g1;
		uint32_t f_total_reach_g2;
		//uint32_t f_new_graph = l_graph_shatter_size - l_graph_shatter_size_begin;

		//Only recalculate reach if origin graph forms 2 graphs ( bridge exists)
		if (l_graph_shatter_size - l_graph_shatter_size_begin > 1 ) {
			for (int i = l_graph_shatter_size_begin; i < l_graph_shatter_size - 1; i++) {
				f_total_reach_g1 = p_graph_shatter[i].getTotalReach();
				f_total_reach_g2 = 0;
				for (int j = i+1; j < l_graph_shatter_size; j++) {
					f_total_reach_g2 += p_graph_shatter[j].getTotalReach();
				}

				auto u_bridge = p_graph_shatter[i].graph_connected_.u_start;
				auto v_bridge = p_graph_shatter[i].graph_connected_.v_end;
				// std::cout << "G1 " << i << " GraphId= " << p_graph_shatter[i].graph_id_ << " " << l_graph_shatter_size - l_graph_shatter_size_begin  << " " << " GraphConnected=" << p_graph_shatter[i].graph_connected_.graph_id << " size reach=" << p_graph_shatter[i].reach_.size() << std::endl;
				// std::cout << "G2 " << i << " GraphId= " << p_graph_shatter[i+1].graph_id_ << " " << l_graph_shatter_size - l_graph_shatter_size_begin << " " << " GraphConnected=" << p_graph_shatter[i+1].graph_connected_.graph_id << " size reach=" << p_graph_shatter[i+1].reach_.size() << std::endl;
				// std::cout << "UB=" << u_bridge << " VB=" << v_bridge << " IREACH=" << p_graph_shatter[i].reach_.size() << " PIREACH=" << p_graph_shatter[p_graph_shatter[i].graph_connected_.graph_id].reach_.size() << std::endl;
				auto u_new_reach = p_graph_shatter[i].getReach(u_bridge) + f_total_reach_g2;
				auto v_new_reach = p_graph_shatter[p_graph_shatter[i].graph_connected_.graph_id].getReach(v_bridge) + f_total_reach_g1;

				//Update reach bridge
				p_graph_shatter[i].setReach(u_bridge, u_new_reach);
				p_graph_shatter[p_graph_shatter[i].graph_connected_.graph_id].setReach(v_bridge, v_new_reach);

				//Update BC score
				g_bet_cen[u_bridge] += (f_total_reach_g1 - 1) * f_total_reach_g2;
				g_bet_cen[v_bridge] += (f_total_reach_g2 - 1) * f_total_reach_g1;
			}
		}

		l_graph_shatter_size_begin = l_graph_shatter_size;
	}
	return true;
}

//DFS travel
bool dfsTravelArticulation(Graph &g, std::set<uint32_t> &p_articulation_vertices) {
	uint32_t g_time = 0;
	//init vector of Node in graph by size
	std::vector<Vertex> l_node(g.edge_list_.size());

	std::stack<uint32_t> S;
	bool check = true;
	for (auto i = 0; i < g.edge_list_.size(); i++) {
		if (l_node[i].color == WHITE) {
			uint32_t s = i;
			S.push(s);
			//std::cerr << "test\n";
			l_node[s].color = GRAY;
			l_node[s].discovery_time = ++g_time;
			l_node[s].low = l_node[s].discovery_time;

			while (!S.empty()) {
				uint32_t u = S.top();
				//S.pop();

				check = true;
				for (auto &it : g.edge_list_[u]) {
					if (l_node[it].color == WHITE) {
						//l_node[it].discovery_time = ++g_time;
						S.push(it);
						l_node[it].color = GRAY;
						l_node[it].discovery_time = ++g_time;
						l_node[it].low = l_node[it].discovery_time;
						l_node[it].predecessor = u;
						l_node[u].num_child++;
						check = false;
						break;
					} else if (it != l_node[u].predecessor && u != l_node[it].predecessor) {
						l_node[u].low = std::min(l_node[u].low, l_node[it].discovery_time);
					}
				}

				if (check) {
					if (l_node[u].color == GRAY) {
						l_node[l_node[u].predecessor].low = std::min(l_node[l_node[u].predecessor].low, l_node[u].low);
						l_node[u].color = BLACK;
						l_node[u].finish_time = ++g_time;

						if (l_node[u].predecessor == -1) {
							if (l_node[u].num_child > 1) {
								p_articulation_vertices.insert(u);
							}
						} else {
							//std::cout << i << " " << m_node[i].predecessor << " " << m_node[i].low << " " << m_node[m_node[i].predecessor].discovery_time << std::endl;
							auto l_predecessor = l_node[u].predecessor;
							if (l_node[l_predecessor].predecessor > -1 && l_node[u].low >= l_node[l_predecessor].discovery_time) {
								p_articulation_vertices.insert(l_node[u].predecessor);
							}
						}
					}
					S.pop();
				}
			}
		}
	}

	return true;
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);

    Graph g;
	std::vector<Graph> m_graph_shatter;

	//Read edge list format
    if (!readEdgeList(g,"../data/wiki-Vote.txt", false)) {
        return 0;
    }
	//std::cout << "SIZE=" << g.edge_list_.size() << std::endl;


	//Read motis format
	// if (!readMetis(g, "../data/cond-mat.graph", false)) {
    //     return 0;
    // }

	g_bet_cen.resize(g.edge_list_.size(), 0);
	//std::cout << g.getDirected() << std::endl;
	//g.printGraph();

    uint16_t option = 2;

    while (option) {
        std::cout << "Please chose function\n";
        std::cout << "Press 1 to calculate articulation point\n";
        std::cout << "Press 2 to calculate bridge\n";
        std::cout << "Press 3 to calculate betweenness centrality\n";
		std::cout << "Press 4 to print graph\n";
        std::cout << "Press 0 to exit\n";
        std::cin >> option;

        switch (option) {
            case 0:
                break;

            case 1:
            {
				double start_time_articulation = getCurrentTimeMlsec();
                std::set<uint32_t> m_articulation;
				dfsTravelArticulation(g, m_articulation);
                std::cout << "Articulation vertex:\n";
                for (auto &it : m_articulation){
                    std::cout << "\t" << it << "\n";
                }
                std::cout << std::flush;
				std::cout << "Time articulation=" << getCurrentTimeMlsec() - start_time_articulation << "\n";
                break;
            }

            case 2:
            {
				double start_time_bridge = getCurrentTimeMlsec();
                std::set<std::pair<uint32_t, uint32_t> > m_bridge;
				//std::vector<Graph> m_graph_shatter;
				//m_graph_shatter.push_back(g);
				//m_graph_shatter.clear();
				dfsTravelBridge(g, m_bridge, m_graph_shatter);
                //print bridge
				// std::cout << "SIZE SHATTER=" << m_graph_shatter.size() << "\n";
				//
				// for (auto &it : m_graph_shatter) {
				// 	//std::cout << &it - &m_graph_shatter[0] << "\n";
				// 	std::cout << "Graph ID=" << it.getGraphId() <<"\n";
				// 	it.printGraph();
				// 	it.printBridgeInfo();
				// 	std::cout << "Total reach=" << it.getTotalReach() << std::endl;
				// 	std::cout << "----------\n";
				// }
				//
				// std::cout << "Bridge:\n";
            	// for (auto &it : m_bridge){
            	// 	std::cout << "\t" << it.first << " " << it.second << "\n";
            	// }
				std::cout << "Time Bridge=" << getCurrentTimeMlsec() - start_time_bridge << "\n";
				//option = 0;
				break;
            }

            case 3:
            {
				double start_time_bc = getCurrentTimeMlsec();
                //std::vector<double> m_bet_cen;
				for (auto &it : m_graph_shatter) {
					betweennessCentrality(it, g_bet_cen);
				}
				if (!g.getDirected()) {
					for (auto & it : g_bet_cen){
		            	it /= 2;
	            	}
				}

                std::cout << std::fixed << std::setprecision(6);
            	//print betweenness centrality value
            	for (auto & it : g_bet_cen){
            		std::cout << it << std::endl;
            	}
				std::cout << "Time betweenness centrality=" << getCurrentTimeMlsec() - start_time_bc<< "\n";
				//option = 0;
                break;
            }

			case 4:
			{
				g.printGraph();
				break;
			}

            default:
                std::cout << "Bad input\n";
        }
    }
    return 0;
}
