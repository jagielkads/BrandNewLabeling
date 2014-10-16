#include <vector>
#include <iostream>
#include "../src/BNLabel.h"

int main(int argc, char** argv){	
	BNLabel<0> hdl;
	char* filename = "../data/slashdot";
	std::vector<std::pair<int, int> > es;
	std::vector<std::vector<int> > adj,adj_1;
	std::vector<int> inv, inv_1;

	//pll
	hdl.LoadGraph(filename, es);
	hdl.ConstructAdj(es, adj);
	hdl.DegreeOrdering(adj, inv);
	
	hdl.ConstructHDIndex(adj, inv);		
}
