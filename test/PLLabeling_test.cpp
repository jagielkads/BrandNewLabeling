#include <vector>
#include <iostream>
#include "../src/BNLabel.h"

int main(int argc, char** argv){	
	BNLabel<0> pll,pll_1;
	char* filename = "slashdot";
	std::vector<std::pair<int, int> > es;
	std::vector<std::vector<int> > adj,adj_1;
	std::vector<int> inv, inv_1;

	//pll
	pll.LoadGraph(filename, es);
	pll.ConstructAdj(es, adj);
	pll.DegreeOrdering(adj, inv);
	
	double indextime = -pll.GetCurrentTimeSec();
	pll.ConstructPPLIndex(adj, inv);		
	std::cout<<indextime+pll.GetCurrentTimeSec()<<std::endl;
	std::cout<<pll.StoreInfo("R/info1")<<std::endl;

	//pll_1
	pll_1.LoadGraph(filename, es);
	pll_1.ConstructAdj(es, adj_1);
	inv_1.resize(adj_1.size());
	for(int v = 0; v < inv_1.size(); v++) inv_1[v] = v;	
	indextime = -pll_1.GetCurrentTimeSec();
	pll_1.ConstructPPLIndex(adj_1, inv_1);		
	std::cout<<indextime+pll_1.GetCurrentTimeSec()<<std::endl;
	std::cout<<pll_1.StoreInfo("R/info2")<<std::endl;

	srand(pll.GetCurrentTimeSec());
	int N = 3000;
	std::vector<std::pair<double, double> > dt(N);	
	for(int n = 0; n < N; n++){
		int v = rand()%pll.GetNumVertices();
		int w = rand()%pll.GetNumVertices();
		double t1 = -pll.GetCurrentTimeSec();
		int d1 = pll.QueryDistance(v, w);
		t1 += pll.GetCurrentTimeSec();
		double t2 = -pll_1.GetCurrentTimeSec();
		int d2 = pll_1.QueryDistance(v, w);
		t2 += pll_1.GetCurrentTimeSec();
		dt[n] = std::make_pair(t1,t2);
		std::cout<<(d1==d2)<<"\t"<<(t2-t1)/t2<<std::endl;
	}	
	std::ofstream ofs("R/querytime");	
	for(int n = 0; n < dt.size(); n++){
		ofs << dt[n].first << "\t" << dt[n].second << "\n";
	}
	ofs.close();
}
