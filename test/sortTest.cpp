#include "../src/BNLabel.h"

int main(){
	BNLabel<0> bnl;
	std::vector<int> numv(0);
	//for(int i = 0; i < 20; i=i+2){
	//	numv.push_back(i);
	//}
	for(int i = 0; i < numv.size(); i++){
		std::cout<<numv[i]<<",";
	}
	std::cout<<std::endl;
	int n = 19;
	int i1 = bnl.binary_search(numv,n,0,numv.size());
	std::vector<int>::iterator it = numv.begin();
	numv.insert(it+i1,n);

	for(int i = 0; i < numv.size(); i++){
		std::cout<<numv[i]<<",";
	}
	std::cout<<std::endl;

	std::pair<std::vector<int>, std::vector<uint8_t> > arr;
	for(int i = 0 ; i < 30; i ++){
		arr.first.push_back(29-i);
		arr.second.push_back(i);
	}
	bnl.quick_sort(arr,0,arr.first.size());
	for(int i = 0; i < arr.first.size(); i++){
		std::cout<<arr.first[i]<<",";
	}
	std::cout<<std::endl;
	for(int i = 0; i < arr.second.size(); i++){
		std::cout<<int(arr.second[i])<<",";
	}
	std::cout<<std::endl;
}	
