#include <vector>
#include <iostream>
#include "../src/BNLabel.h"

int main(int argc, char** argv){	
	BNLabel<0> hdl,pll;
	//char* filename = "/Users/jagsly/Data/prjs/2_hop_shortest_path/BrandNewLabeling/data/slashdot";
	char* filename = argv[1];
	if(argc!=4){
		std::cout<<"a.out [srcPath] [outPath1] [outPath2]"<<std::endl;
		return 0;
	}
		
	std::vector<std::pair<int, int> > es;
	std::vector<std::vector<int> > adj,adj_1;
	std::vector<int> inv, inv_1;

	
	//pll
	pll.LoadGraph(filename, es);
	pll.ConstructAdj(es, adj);
	pll.DegreeOrdering(adj, inv);
	
	pll.ConstructPPLIndex(adj, inv);		
	//pll.StoreIndex(argv[2]);

	//hd
	hdl.LoadGraph(filename,es);
	hdl.ConstructAdj(es, adj);
	hdl.DegreeOrdering(adj, inv);
	
	hdl.ConstructHDIndex(adj, inv);
	//hdl.StoreIndex(argv[3]);	

	srand(pll.GetCurrentTimeSec());
	int N = 100;
	std::vector<std::pair<double, double> > dt(N);	
	for(int n = 0; n < N; n++){
		int v = rand()%pll.GetNumVertices();
		int w = rand()%pll.GetNumVertices();
		double t1 = -pll.GetCurrentTimeSec();
		int d1 = pll.QueryDistance(v, w);
		t1 += pll.GetCurrentTimeSec();
		t1 = -pll.GetCurrentTimeSec();
		d1 = pll.QueryDistance(v, w);
		t1 += pll.GetCurrentTimeSec();
		double t2 = -hdl.GetCurrentTimeSec();
		int d2 = hdl.QueryDistance(v, w);
		t2 += hdl.GetCurrentTimeSec();
		t2 = -hdl.GetCurrentTimeSec();
		d2 = hdl.QueryDistance(v, w);
		t2 += hdl.GetCurrentTimeSec();
		//std::cout<<(d1==d2)<<"\t"<<(t2-t1)/t2<<std::endl;
		std::cout<<t1<<"\t"<<t2<<std::endl;
	}	
}
