#include <vector>
#include <iostream>
#include "../src/BNLabel.h"

int main(int argc, char** argv){	
	BNLabel<0> pll;
	//char* filename = "/Users/jagsly/Data/prjs/2_hop_shortest_path/BrandNewLabeling/data/slashdot";
	char* filename = argv[1];
	std::cout<<filename<<std::endl;
		
	std::vector<std::pair<int, int> > es;
	std::vector<std::vector<int> > adj,adj_1;
	std::vector<int> inv, inv_1;

	//pll
	pll.LoadGraph(filename, es);
	pll.ConstructAdj(es, adj);
	pll.DegreeOrdering(adj, inv);
	
	pll.ConstructPLLIndex(adj, inv);		
	char* labelName = strcat(argv[1],".label.pll");
	pll.StoreIndex(labelName);
//	std::cout<<pll.QueryDistance(0,100)<<std::endl;
}
