#include <iostream>
#include <vector>
#include <set>
#include <fstream>
#include <stack>
#include <queue>
#include <iomanip>

//Travel color, WHITE: unvisited, GRAY: visited but not finish, BLACK: finished
enum eColor {
	WHITE,
	GRAY,
	BLACK
};


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

class Graph{
   // uint32_t v_num_;
   // uint32_t e_num_;
    std::vector<std::vector<uint32_t> > edge_list_;
    std::set<uint32_t> vertex_list_;
    bool directed_ = false;

    public:

    // Graph(bool p_directed) {
    //     directed_ = p_directed;
    // }

    // uint32_t getNumberOfVectex() {
    //     return vertex_list_.size();
    // }
    //
    // uint32_t getNumberOfEdge() {
    //     return edge_list_.size();
    // }

    bool readEdgeList(std::string p_file_name, bool p_directed = false){
        directed_ = p_directed;
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
            vertex_list_.insert(it.first);
            vertex_list_.insert(it.second);

    		//std::cout << graph.size();
    		if (edge_list_.size() < (it.first + 1)) {
    			edge_list_.resize(it.first + 1);
    		}
    		if (edge_list_.size() < (it.second + 1)) {
    			edge_list_.resize(it.second + 1);
    		}

    		edge_list_[it.first].push_back(it.second);
    		if (!directed_) {
    			edge_list_[it.second].push_back(it.first);
    		}
    	}

    	// //sort adjacency list
    	for (auto &it : edge_list_) {
    		std::sort(it.begin(), it.end());
    	}
        return true;
    }

    //print graph
    void printGraph() {
    	for (auto &it : edge_list_){
    		for (auto &it1 : it) {
    			std::cout << &it - &edge_list_[0] << " " << it1 << std::endl;
    		}
    	}
    }

    //DFS travel
    std::vector<Vertex> dfsTravel() {
        uint32_t g_time = 0;
        //init vector of Node in graph by size
    	std::vector<Vertex> l_node(edge_list_.size());

        std::stack<uint32_t> S;
    	bool check = true;
    	for (auto i = 0; i < edge_list_.size(); i++) {
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
                    for (auto &it : edge_list_[u]) {
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
                        } else if (it != l_node[u].predecessor) {
    						l_node[u].low = std::min(l_node[u].low, l_node[it].discovery_time);
    					}
                    }

    				if (check) {
    					S.pop();
    					if (l_node[u].color == GRAY) {
    						l_node[l_node[u].predecessor].low = std::min(l_node[l_node[u].predecessor].low, l_node[u].low);
    						l_node[u].color = BLACK;
                        	l_node[u].finish_time = ++g_time;
    					}
    				}
                }
            }
    	}

        return l_node;
    }

    std::set<uint32_t> getArticulationVertex() {
        std::set<uint32_t> l_articulation;
        std::vector<Vertex> l_node = dfsTravel();

        //Find articulation vertex
    	for (auto i = 0; i < edge_list_.size(); i++){
    		if (l_node[i].predecessor == -1) {
    			if (l_node[i].num_child > 1) {
    				l_articulation.insert(i);
    			}
    		} else {
				if (l_node[l_node[i].predecessor].predecessor > -1 && l_node[i].low >= l_node[l_node[i].predecessor].discovery_time) {
					l_articulation.insert(l_node[i].predecessor);
				}
    		}
    	}
        return l_articulation;
    }

    std::set<std::pair<uint32_t, uint32_t> > getBridge() {
        std::vector<Vertex> l_node = dfsTravel();
        std::set<std::pair<uint32_t, uint32_t> > l_bridge;
        //Find bridge
    	for (auto i = 0; i < edge_list_.size(); i++) {
			if (l_node[i].predecessor == -1)
				continue;
			if (l_node[i].low > l_node[l_node[i].predecessor].discovery_time) {
				l_bridge.insert(std::pair<uint32_t, uint32_t>(l_node[i].predecessor, i));
			}
    	}
        return l_bridge;
    }

    //Compute betweenness centrality
    std::vector<double> betweennessCentrality() {
        std::vector<double> l_bet_cen(edge_list_.size());

        //set all the value of betweenness centrality of all vectices equal to 0 before computing
    	for (auto &it : l_bet_cen){
    		it = 0;
    	}

        //travel all vectice
    	for (size_t i = 0; i < edge_list_.size(); i++) {
    		//for (size_t j = 0; j < vertex_list_[i].size(); j++) {
    			uint32_t s = i;
    			//std::cerr << vertex_list_.size() << std::endl;
    			//sleep(1);

    			uint32_t *check_map = new uint32_t [edge_list_.size()/32 + 1];
    			std::memset(check_map, 0, (edge_list_.size()/32+1)*sizeof(uint32_t));

    			std::stack<uint32_t> S; //visited vertices
    			std::vector<std::vector<uint32_t> > predecessor(edge_list_.size());
    			std::queue<uint32_t> Q; //unvisited vertices
    			std::vector<int32_t> number_of_shortest_path(edge_list_.size());
    			for (auto &it : number_of_shortest_path){
    				it = 0;
    			}

    			number_of_shortest_path[s] = 1;

    			std::vector<int64_t> distance(edge_list_.size());
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
    				for (size_t k = 0; k < edge_list_[v].size(); k++) {
    					uint32_t w = edge_list_[v][k];
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

    			std::vector<double> dependency(edge_list_.size());
    			for (auto &it : dependency) {
    				it = 0;
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
    					l_bet_cen[w] += dependency[w];
    				}
    			}
    			delete [] check_map;
    	}

    	// divide by 2 for undirected graph
        if (!directed_) {
            for (auto &it : l_bet_cen) {
                it /= 2;
            }
        }
        return l_bet_cen;
    }
};

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);

    Graph g;
    if (!g.readEdgeList("../data/test2.txt", false)) {
        return 0;
    }

    uint16_t option = 1;

    while (option) {
        std::cout << "Please chose function\n";
        std::cout << "Press 1 to calculate articulation point\n";
        std::cout << "Press 2 to calculate bridge\n";
        std::cout << "Press 3 to calculate betweenness centrality\n";
        std::cout << "Press 0 to exit\n";
        std::cin >> option;

        switch (option) {
            case 0:
                break;

            case 1:
            {
                std::set<uint32_t> m_articulation = g.getArticulationVertex();
                std::cout << "Articulation vertex:\n";
                for (auto &it : m_articulation){
                    std::cout << "\t" << it << "\n";
                }
                std::cout << std::flush;
                break;
            }

            case 2:
            {
                std::set<std::pair<uint32_t, uint32_t> > m_bridge = g.getBridge();
                //print bridge
            	std::cout << "Bridge:\n";
            	for (auto &it : m_bridge){
            		std::cout << "\t" << it.first << " " << it.second << "\n";
            	}
                break;
            }

            case 3:
            {
                std::vector<double> m_bet_cen = g.betweennessCentrality();
                std::cout << std::fixed << std::setprecision(6);
            	//print betweenness centrality value
            	for (auto & it : m_bet_cen){
            		std::cout << it << std::endl;
            	}
                break;
            }

            default:
                std::cout << "Bad input\n";
        }
    }
    // std::vector<double> m_bet_cen = g.betweennessCentrality();
    // std::cout << std::fixed << std::setprecision(6);
	// //print betweenness centrality value
	// for (auto & it : m_bet_cen){
	// 	std::cout << it << std::endl;
	// }

    // std::set<uint32_t> m_articulation = g.getArticulationVertex();
    // std::cout << "Articulation vertex:\n";
    // for (auto &it : m_articulation){
    //     std::cout << "\t" << it << "\n";
    // }
    // std::cout << std::flush;

    // std::set<std::pair<uint32_t, uint32_t> > m_bridge = g.getBridge();
    // //print bridge
	// std::cout << "Bridge:\n";
	// for (auto &it : m_bridge){
	// 	std::cout << "\t" << it.first << " " << it.second << "\n";
	// }

    return 0;
}
