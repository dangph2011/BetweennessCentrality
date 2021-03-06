#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstdio>
#include <stack>
#include <queue>
#include <cstring>
#include <unistd.h>
#include <iomanip>
#include <bitset>

//print graph
void print_graph(std::vector<std::vector<unsigned> > graph) {
	for (auto &it : graph){
		for (auto &it1 : it) {
			std::cout << &it - &graph[0] << " " << it1 << std::endl;
		}
	}
}

//set bit at index in the array addr
static inline void set_bit(unsigned index, unsigned *addr)
{
	unsigned p = index >> 5, x = index & 31; //& ~ %31
	addr[p] |= 1 << x;
}

//test bit at index in the array addr
static inline int test_bit(unsigned index, const unsigned *addr)
{
	unsigned p = index >> 5, x = index & 31;
	return (addr[p] >> x) & 1;
}

//clear map at positions specified by queue array
// static inline void reset_map(unsigned *map, unsigned* queue, unsigned top)
// {
// 	for (unsigned i = 0; i < top; i++)	map[queue[i]>>5] = 0;
// 	//map[queue[0:top]>>5] = 0;
// }

int main() {
	std::ifstream f_in;
	std::string line;
	bool directed = false;
	std::vector<std::vector<unsigned> > graph;
	unsigned u, v;
	std::vector<std::pair<unsigned, unsigned> > edgesList;

	f_in.open("../data/test2.txt");
	if(f_in.fail())
	{
	    std::cout << "Error: Opening file";
		return 0;
	}
	//read data from file
    while(f_in >> u >> v){
		//put to vector
		edgesList.push_back(std::pair<unsigned, unsigned>(u,v));

    }

	//Store in adjacency list
    for (auto &it : edgesList) {
		//std::cout << graph.size();
		if (graph.size() < (it.first + 1)) {
			graph.resize(it.first + 1);
		}
		if (graph.size() < (it.second + 1)) {
			graph.resize(it.second + 1);
		}

		graph[it.first].push_back(it.second);
		if (!directed) {
			graph[it.second].push_back(it.first);
		}
	}

	// //sort adjacency list
	for (auto &it : graph) {
		std::sort(it.begin(), it.end());
	}

	//test print graph
	//print_graph(graph);

	//Compute betweenness centrality
	std::vector<double> betCen(graph.size());
	//set all the value of betweenness centrality of all vectices equal to 0
	for (auto &it : betCen){
		it = 0;
	}

	//travel all vectice
	for (size_t i = 0; i < graph.size(); i++) {
		//for (size_t j = 0; j < graph[i].size(); j++) {
			unsigned s = i;
			//std::cerr << graph.size() << std::endl;
			//sleep(1);

			unsigned *checkMap = new unsigned [graph.size()/32 + 1];
			std::memset(checkMap, 0, (graph.size()/32+1)*sizeof(unsigned));

			std::stack<unsigned> S; //visited vertices
			std::vector<std::vector<unsigned> > predecessor(graph.size());
			std::queue<unsigned> Q; //unvisited vertices
			std::vector<int> numberOfShortestPath(graph.size());
			for (auto &it : numberOfShortestPath){
				it = 0;
			}

			numberOfShortestPath[s] = 1;

			std::vector<long int> distance(graph.size());
			for (auto &it : distance){
				it = -1;
			}

			distance[s] = 0;
			Q.push(s);
			set_bit(s, checkMap);
			while (!Q.empty()) {
				unsigned v = Q.front();
				Q.pop();
				S.push(v);
				//find all neighbor w of v
				for (size_t k = 0; k < graph[v].size(); k++) {
					unsigned w = graph[v][k];
					//check if w is visited vertex

					if (!test_bit(w, checkMap)) {
						set_bit(w, checkMap);	//set w is visited vertex
						if (distance[w] < 0) {
							Q.push(w);
							distance[w] = distance[v] + 1;
						}
					}

					if (distance[w] == distance[v] + 1) {
						numberOfShortestPath[w] += numberOfShortestPath[v];
						predecessor[w].push_back(v);
 					}
				}

			}

			std::vector<double> dependency(graph.size());
			for (auto &it : dependency) {
				it = 0;
			}
			//S returns vertices in order of non-increasing distance from s
			while (!S.empty()) {
				unsigned w = S.top();
				S.pop();
				for (auto &it : predecessor[w]) {
					int v = it;
					dependency[v] += ((double)numberOfShortestPath[v] / (double)numberOfShortestPath[w]) * (1 + (double)dependency[w]);
				}

				if (w != s) {
					betCen[w] += dependency[w];
				}
			}
			delete [] checkMap;
	}

	std::cout << std::fixed << std::setprecision(6);
    if (!directed) {
        for (auto &it : betCen) {
            it /= 2;
        }
    }
	//print betweenness centrality value
	for (auto & it : betCen){
		std::cout << it << std::endl;
	}

	return 0;
}
