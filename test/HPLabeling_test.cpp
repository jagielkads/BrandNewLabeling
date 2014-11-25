#include <vector>
#include <iostream>
#include "../src/BNLabel.h"

int main(int argc, char** argv){	
	BNLabel<0> hdl;
	//char* filename = "/Users/jagsly/Data/prjs/2_hop_shortest_path/BrandNewLabeling/data/slashdot";
	char* filename = argv[1];
	std::cout<<filename<<std::endl;
		
	std::vector<std::pair<int, int> > es;
	std::vector<std::vector<int> > adj,adj_1;
	std::vector<int> inv, inv_1;

	//pll
	hdl.LoadGraph(filename, es);
	hdl.ConstructAdj(es, adj);
	hdl.DegreeOrdering(adj, inv);
	
	hdl.ConstructHDIndex(adj, inv, argv[1]);		
	char* labelName = strcat(argv[1],".label");
	hdl.StoreIndex(labelName);
	std::cout<<hdl.QueryDistance(0,100)<<std::endl;
}
