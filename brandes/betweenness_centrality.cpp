#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstdio>
#include <stack>
#include <queue>
#include <cstring>

//print graph
void print_graph(std::vector<std::vector<unsigned> > graph) {
	for (std::vector<std::vector<unsigned> >::iterator it = graph.begin(); it !=graph.end(); it++){
		for (std::vector<unsigned>::iterator it1 = it->begin(); it1 != it->end(); it1++) {
			std::cout << it - graph.begin() << " " << (*it1) << std::endl;
		}
	}
}

//set bit at index in the array addr
static inline void set_bit(unsigned index, unsigned *addr)
{
	unsigned p = index >> 5, x = index & 31;
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
	std::vector<std::vector<unsigned> > graph;
	f_in.open("../data/facebook_combined.txt");
	if(f_in.fail())
	{
	    std::cout << "Error: Opening file";
		return 0;
	}
	unsigned u, v;
	std::vector<std::pair<unsigned, unsigned> > edgesList;
	//read data from file
    while(f_in >> u >> v){
		//put to vector
		edgesList.push_back(std::pair<unsigned, unsigned>(u,v));
    }

	//Store in adjacency list
    for (std::vector<std::pair<unsigned, unsigned> >::iterator it = edgesList.begin(); it!= edgesList.end(); it++) {
		//std::cout << graph.size();
		if (graph.size() < (*it).first + 1) {
			graph.resize((*it).first+1);
		}
		graph[(*it).first].push_back((*it).second);
	}

	//sort adjacency list
	for (std::vector<std::vector<unsigned> >::iterator it = graph.begin(); it !=graph.end(); it++) {
		std::sort(it->begin(), it->end());
	}

	//test print graph
	//print_graph(graph);

	//Compute betweenness centrality
	std::vector<double> betCen(graph.size());
	//set all the value of betweenness centrality of all vectices equal to 0
	for (std::vector<double>::iterator it = betCen.begin(); it != betCen.end(); it++){
		(*it) = 0;
	}
	for (size_t i = 0; i < graph.size(); i++) {
		for (size_t j = 0; j < graph[i].size(); j++) {
			unsigned *checkMap = new unsigned [graph.size()/4 + 4];
			std::memset(checkMap, 0, (graph.size()/4 + 4)*sizeof(unsigned));

			std::stack<unsigned> S; //visited vertices
			std::vector<std::vector<unsigned> > predecessor(graph.size());
			std::queue<unsigned> Q; //unvisited vertices
			std::vector<int> numberOfShortestPath(graph.size());
			for (std::vector<int>::iterator it = numberOfShortestPath.begin(); it != numberOfShortestPath.end(); it++){
				(*it) = 0;
			}

			int s = graph[i][j];
			numberOfShortestPath[s] = 1;

			std::vector<int> distance(graph.size());
			for (std::vector<int>::iterator it = numberOfShortestPath.begin(); it != numberOfShortestPath.end(); it++){
				(*it) = -1;
			}

			distance[s] = 0;
			Q.push(s);

			while (!Q.empty()) {
				unsigned v = Q.front();
				Q.pop();
				set_bit(v, checkMap);
				S.push(v);
				//find all neighbor w of v
				for (size_t k = 0; k < graph[v].size(); k++) {
					unsigned w = graph[v][k];
					if (!test_bit(w, checkMap)) {
						if (distance[w] < 0) {
							Q.push(w);
							distance[w] = distance[v] + 1;
						}

						if (distance[w] == distance[v] + 1) {
							numberOfShortestPath[w] += numberOfShortestPath[v];
							predecessor[w].push_back(v);
	 					}
					}
				}
			}

			std::vector<int> dependency(graph.size());
			for (std::vector<int>::iterator it = dependency.begin(); it != dependency.end(); it++) {
				(*it) = 0;
			}
			//S returns vertices in order of non-increasing distance from s
			while (!S.empty()) {
				unsigned w = S.top();
				S.pop();
				for (std::vector<unsigned>::iterator it = predecessor[w].begin(); it != predecessor[w].end(); it++) {
					int v = (*it);
					dependency[v] += (1 / predecessor[w].size()) * (1 + dependency[w]);
				}

				if (w != s) {
					betCen[w] += dependency[w];
				}
			}
			delete [] checkMap;
		}
	}

	return 0;
}
