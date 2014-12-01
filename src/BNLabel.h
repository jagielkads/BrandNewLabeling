#ifndef PLLABELING_H
#define PLLABELING_H

#include <fstream>
#include <vector>
#include <algorithm>
#include <climits>
#include <sys/time.h>
#include <xmmintrin.h>
#include <utility>
#include <iostream>

template<int kNumBitParallelRoots = 50>
class BNLabel {
	public:
	bool LoadGraph(const char *filename, std::vector<std::pair<int, int> > &es);
	bool ConstructAdj(const std::vector<std::pair<int, int> > &es, std::vector<std::vector<int> > &adj);
	bool ConstructPLLIndex(const std::vector<std::vector<int> > &adj, const std::vector<int> &inv);
	bool ConstructHDIndex(const std::vector<std::vector<int> > &adj, const std::vector<int> &inv);
	bool DegreeOrdering(std::vector<std::vector<int> > &adj, std::vector<int> &inv);
	bool GraphUpdate(std::vector<std::pair<int, int> > &es, const int num_upd);

	bool StoreInfo(const char *filename);
	bool StoreIndex(const char *filename);
	bool LoadIndex(const char *filename);
	
	inline int QueryDistance(int v, int w);

	int GetNumVertices() { return num_v_; }	

	static const uint8_t INF8;
	
	void Free();
	BNLabel()
		: num_v_(0), index_(NULL), time_load_(0), time_indexing_(0) {}
	virtual ~BNLabel() {
		Free();
	}
	
	double GetCurrentTimeSec() {
		struct timeval tv;
		gettimeofday(&tv, NULL);
		return tv.tv_sec + tv.tv_usec * 1e-6;
	}
	double time_load_, time_indexing_;

//	private:	
	int num_v_;


	struct index_t{
		uint8_t bpspt_v[kNumBitParallelRoots];
		uint64_t bpspt_d[kNumBitParallelRoots];
		uint32_t *spt_v;
		uint8_t *spt_d;
	} __attribute__((aligned(64)));
	
	index_t *index_;
	
	public:	
	void quick_sort(std::pair<std::vector<int>, std::vector<uint8_t> > &pair_array, int first, int last);
	int binary_search(const std::vector<int> &array, const int &key, const int &low, const int & high);
};

template<int kNumBitParallelRoots>
const uint8_t BNLabel<kNumBitParallelRoots>::INF8 = 100;

template<int kNumBitParallelRoots>
bool BNLabel<kNumBitParallelRoots>
::LoadGraph(const char *filename, std::vector<std::pair<int, int> > &es){
	std::ifstream ifs(filename);
	if(ifs == false)
		return false;
	for (int v,w; ifs>> v >> w; ) {
		es.push_back(std::make_pair(v,w));
	}	
	if (ifs.bad()) return false;
	int &V = this->num_v_;
	for( int i = 0; i < es.size(); i++){
		V = std::max(V, std::max(es[i].first, es[i].second) + 1 );
	}
	return true;
}

template<int kNumBitParallelRoots>
bool BNLabel<kNumBitParallelRoots>
::ConstructAdj(const std::vector<std::pair<int, int> > &es, std::vector<std::vector<int> > &adj){
	int &V = this->num_v_;
	if(adj.size() != V) adj.resize(V);
	for(int i = 0; i < es.size(); i++){
		int v = es[i].first, w = es[i].second;
		if (v != w){
			adj[v].push_back(w);
			//adj[w].push_back(v);
		}
	}
	
	index_ = (index_t*) malloc( V * sizeof(index_t));

	if(index_ == NULL){
		num_v_ = 0;
		return false;
	}
	for(int v = 0; v < V; v++){
		index_[v].spt_v = NULL;
		index_[v].spt_d = NULL;
	}

	return true;
}

template<int kNumBitParallelRoots>
bool BNLabel<kNumBitParallelRoots>
::GraphUpdate(std::vector<std::pair<int, int> > &es, const int num_upd){	
	int &V = this->num_v_;
	if(num_upd >= V || num_upd <= 0) return false;
	int c_count = num_upd;
	srand(this->GetCurrentTimeSec());
	while(c_count--){
		int flag = rand()%2;
		if (flag == 0){//Insert Edge
			int v = rand()%V;
			int u = rand()%V;
			while(v==u) u = rand()%V;
			for(int i = 0; i < es.size(); i++){
				if((es[i].first==v && es[i].second==u) || (es[i].first==u&&es[i].second==v)){
					v = rand()%V;
					u = rand()%V;	
					while(v==i) u = rand()%V;
					i = 0;       		
				}
			}
		} else if (flag == 1){//Delete Edge
			int del_idx = rand()%es.size();
			int v = es[del_idx].first;
			int u = es[del_idx].second;
			es.erase(es.begin()+del_idx);
			for(int i = 0; i < es.size(); i++){
				if(es[i].first==u && es[i].second==v){
					es.erase(es.begin()+i);
					break;
				}
			}
		} //else if (flag == 2){//Insert Edge from a random point to another
		//} else { //if (flag == 3){ //Delete Edge
		//}
	}
	return true;
}

template<int kNumBitParallelRoots>
bool BNLabel<kNumBitParallelRoots>
::DegreeOrdering(std::vector<std::vector<int> > &adj, std::vector<int> &inv){	
	int &V = this->num_v_;
	std::vector<std::pair<int, int> > deg(this->num_v_);
	for(int i = 0; i < V; i++){
		deg[i] = std::make_pair(adj[i].size(), i);	
	}
	std::sort(deg.rbegin(), deg.rend());
	inv.resize(V);
	for (int i = 0; i < V; ++i) inv[i] = deg[i].second;
	// Relabel the vertex IDs
	std::vector<int> rank(V);
	for (int i = 0; i < V; ++i) rank[deg[i].second] = i;
	std::vector<std::vector<int> > new_adj(V);
	for (int v = 0; v < V; ++v) {
		for (size_t i = 0; i < adj[v].size(); ++i) {
			new_adj[rank[v]].push_back(rank[adj[v][i]]);
		}   
	}   
        adj.swap(new_adj);	
	return true;
}

template<int kNumBitParallelRoots>
bool BNLabel<kNumBitParallelRoots>
::ConstructPLLIndex(const std::vector<std::vector<int> > &adj, const std::vector<int> &inv){
	int &V = num_v_;
	std::vector<std::pair<std::vector<int>, std::vector<uint8_t> > >
		tmp_idx(V, make_pair(std::vector<int>(1, V), std::vector<uint8_t>(1, INF8))); //Labels
	std::vector<bool> vis(V);
	std::vector<int> que(V);
	std::vector<uint8_t> dst_r(V +1, INF8);
	
	std::vector<bool> usd(V,false);	
	for(int r = 0; r < V; r++){ //outter loop ,the r-th pruned
		if(usd[r]) continue;
		index_t &idx_r = index_[inv[r]]; //T
		const std::pair<std::vector<int>, std::vector<uint8_t> >
			&tmp_idx_r = tmp_idx[r]; //T
		for (size_t i = 0; i < tmp_idx_r.first.size(); i++){
			dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i]; //T
		}	
		
		int que_t0 = 0, que_t1 = 0, que_h = 0;
		que[que_h++] = r;
		vis[r] = true;
		que_t1 = que_h;
		
		double etime = -GetCurrentTimeSec();
		for(int d = 0; que_t0 < que_h; d++){
			for(int que_i = que_t0; que_i < que_t1; que_i++){
				int v = que[que_i];
				std::pair<std::vector<int>, std::vector<uint8_t> >	&tmp_idx_v = tmp_idx[v];
				//index_t &idx_v = index_[inv[v]];	//this one only used in bit parallel
				
				_mm_prefetch(&tmp_idx_v.first[0], _MM_HINT_T0);
				_mm_prefetch(&tmp_idx_v.second[0], _MM_HINT_T0);	
				if (usd[v]) continue;
				for(size_t i = 0; i < tmp_idx_v.first.size(); i++){
					int w = tmp_idx_v.first[i];
					int td = tmp_idx_v.second[i] + dst_r[w];
					if (td <= d) goto pruned;
				}
				tmp_idx_v.first.back() = r;
				tmp_idx_v.second.back() = d;
				tmp_idx_v.first.push_back(V);
				tmp_idx_v.second.push_back(INF8);
				for(size_t i = 0; i < adj[v].size(); i++){

					int w = adj[v][i];
					if(!vis[w]) {
						que[que_h++] = w;
						vis[w] = true;
					}
				}

				pruned:
					{}
			}
			que_t0 = que_t1;
			que_t1 = que_h;
		}
		//std::cout<<etime+GetCurrentTimeSec()<<std::endl;
		for(int i = 0; i < que_h; i++) vis[que[i]] = false;
		for(size_t i = 0; i < tmp_idx_r.first.size(); i++){
			dst_r[tmp_idx_r.first[i]] = INF8;
		}
		usd[r] = true;
	}
	//Appending Labels
	for(int v = 0; v < V; v++){
		int k = tmp_idx[v].first.size();
		
		index_[inv[v]].spt_v = (uint32_t*) malloc( k * sizeof(uint32_t));
		index_[inv[v]].spt_d = (uint8_t*) malloc( k * sizeof(uint8_t));
		
		if(!index_[inv[v]].spt_v || !index_[inv[v]].spt_d ){
			Free();
			return false;
		}
		for(int i = 0; i < k; i++) index_[inv[v]].spt_v[i] = tmp_idx[v].first[i];
		for(int i = 0; i < k; i++) index_[inv[v]].spt_d[i] = tmp_idx[v].second[i];
		tmp_idx[v].first.clear();
		tmp_idx[v].second.clear();
	}
	return true;
}

template<int kNumBitParallelRoots>
void BNLabel<kNumBitParallelRoots>
::Free(){
	for(int v = 0; v < num_v_; v++){
		free(index_[v].spt_v);
		free(index_[v].spt_d);
	}
	free(index_);
	index_ = NULL;
	num_v_ = 0;
}

template<int kNumBitParallelRoots>
bool BNLabel<kNumBitParallelRoots>
::StoreIndex(const char *filename){
	//|V|
	//v0 |v0| l1 d1 l2 d2
	int &V = this->num_v_;
	std::ofstream ofs(filename);
	if(ofs == false)
		return false;
	ofs << V << std::endl;
	std::cout << V << std::endl;
	int cnt = 0;
 	for(int v = 0; v < V; v++){
		int count = 0;
		while(1){
			int label_v = index_[v].spt_v[count];
			if(label_v == V)
				break;
			int label_d = index_[v].spt_d[count++];
			cnt++;

		}	
		ofs << count << std::endl;
		//if(count != 0)	
		//	ave_d = (ave_d+0.5)/count;
		//else
		//	ave_d = 0;
			
		//ofs << v<< " " << count<< " " << ave_d<<std::endl;
		//ofs << v<< " " << count<< " " << max_d<<std::endl;
		count = 0;
		//while(1){
		//	int label = index_[v].spt_v[count];
		//	int label_d = index_[v].spt_d[count++];
		//	ofs << " " << label;
		//	ofs << " " << label_d;
		//	if(label == V)
		//		break;
		//}
		//ofs << std::endl;
	}	
	ofs << cnt << std::endl;
	ofs.close();
	if (ofs.bad()) return false;
	return true;
}

template<int kNumBitParallelRoots>
bool BNLabel<kNumBitParallelRoots>
::StoreInfo(const char *filename){
	int &V = this->num_v_;
	std::ofstream ofs(filename);
	if(ofs == false)
		return false;
 	for(int v = 0; v < V; v++){
		int count = 0;
		while(1){
			if(index_[v].spt_v[count++] == V)
				break;
		}		
		ofs << count << std::endl;
	}	
	ofs.close();
	if (ofs.bad()) return false;
	return true;
}

template<int kNumBitParallelRoots>
int BNLabel<kNumBitParallelRoots>
::QueryDistance(int v, int w) {
	if(v >= num_v_ || w >= num_v_) return v == w ? 0 : INT_MAX;

	const index_t &idx_v = index_[v];
	const index_t &idx_w = index_[w];
	int d = INF8;

	_mm_prefetch(&idx_v.spt_v[0], _MM_HINT_T0);
	_mm_prefetch(&idx_v.spt_d[0], _MM_HINT_T0);
	_mm_prefetch(&idx_w.spt_v[0], _MM_HINT_T0);
	_mm_prefetch(&idx_w.spt_d[0], _MM_HINT_T0);

	for (int i1 = 0, i2 = 0 ; ; ) {
		int v1 = idx_v.spt_v[i1], v2 = idx_w.spt_v[i2];
		if(v1 == v2) {
			if(v1 == num_v_) break;
			int td = idx_v.spt_d[i1] + idx_w.spt_d[i2];
			if(td < d) d = td;
			i1++;
			i2++;
		} else {
			i1 += v1 < v2 ? 1 : 0;
			i2 += v1 > v2 ? 1 : 0;
		}
	}

	if ( d >= INF8 - 2) d = INT_MAX;
	return d;
}

template<int kNumBitParallelRoots>
bool BNLabel<kNumBitParallelRoots>
::ConstructHDIndex(const std::vector<std::vector<int> > &adj, const std::vector<int> &inv){
	int &V = num_v_;
	std::vector<std::pair<std::vector<int>, std::vector<uint8_t> > >
		all_idx(V); //all labels
	std::vector<std::pair<std::vector<int>, std::vector<uint8_t> > >
		prev_idx(V); //prev labels
	std::vector<std::pair<std::vector<int>, std::vector<uint8_t> > >
		cur_idx(V); //cur labels
	int pnum;
	//std::vector<bool> vis(V);
	//std::vector<int> que(V);
	//std::vector<uint8_t> dst_r(V +1, INF8);
	
	//to output info of each iteration
	std::vector<std::vector<int> > itInfo;	

	//Initialization: Generating labels based on directed linking	
	pnum = 0;
	for(int v = 0; v < V; v++){
		std::pair<std::vector<int>, std::vector<uint8_t> >
			&all_idx_v = all_idx[v];
		std::pair<std::vector<int>, std::vector<uint8_t> >
			&prev_idx_v = prev_idx[v];
		
		//all_idx_v.first.push_back(v);
		//all_idx_v.second.push_back(0);
		//prev_idx_v.first.push_back(v);
		//prev_idx_v.second.push_back(0);
		//pnum++;	
		for(int r = 0; r < adj[v].size(); r++){
			int w = adj[v][r];
			if(w < v){//only push higher-rank label into V
				all_idx_v.first.push_back(w);
				all_idx_v.second.push_back(1);
				prev_idx_v.first.push_back(w);
				prev_idx_v.second.push_back(1);
				pnum++;
			}
		}
	}	


	int iteration = 0;
	//iteration till the end
	while(pnum!=0){
		std::cout<<++iteration<<" iteration go."<<std::endl;
		pnum = 0;


		double itert,sortt, gent, prunt;
		itert = -GetCurrentTimeSec();
		sortt = -GetCurrentTimeSec();
		//quicksort for all idx
		//for(int u = 0; u < V; u++)
		//	quick_sort(all_idx[u], 0, all_idx[u].first.size()-1);
		for(int v = 0; v < V; v++){
			std::pair<std::vector<int>, std::vector<uint8_t> >	
				&all_idx_v = all_idx[v];
			std::vector<std::pair<int, uint8_t> > all_idx_v_pair(all_idx_v.first.size());
			for(int r = 0; r < all_idx_v.first.size(); r++){
				all_idx_v_pair[r] = std::make_pair(all_idx_v.first[r], all_idx_v.second[r]);
			}
			std::sort(all_idx_v_pair.begin(),all_idx_v_pair.end());
			for(int r = 0; r < all_idx_v.first.size(); r++){
				all_idx_v.first[r] = all_idx_v_pair[r].first;
				all_idx_v.second[r] = all_idx_v_pair[r].second;
			}	
		}
		//inverse label
		std::vector<std::pair<std::vector<int>, std::vector<uint8_t> > >
			inv_all_idx(V);
		for(int u = 0; u < V; u++){
			for(int v = 0; v < all_idx[u].first.size(); v++){
				inv_all_idx[all_idx[u].first[v]].first.push_back(u);
				inv_all_idx[all_idx[u].first[v]].second.push_back(all_idx[u].second[v]);
			}
		}
		std::cout<<"Label sorting time:"<<sortt+GetCurrentTimeSec()<<std::endl;

		gent = -GetCurrentTimeSec();
		//Label Generation	
		for(int u = 0; u < V; u++){/*{{{*/
			std::pair<std::vector<int>, std::vector<uint8_t> >
				&prev_idx_u = prev_idx[u];		
			std::pair<std::vector<int>, std::vector<uint8_t> >
				&all_idx_u = all_idx[u];
			std::pair<std::vector<int>, std::vector<uint8_t> >
				&inv_all_idx_u = inv_all_idx[u];
			
			for(int i1 = 0; i1 < prev_idx_u.first.size(); i1++){
				int v1 = prev_idx_u.first[i1]; // L(u): u <-> v1
				uint8_t d1 = prev_idx_u.second[i1];
				
				//Generating by Rule 1
				double rule1time = -GetCurrentTimeSec();
				for(int i2 = 0; i2 < all_idx_u.first.size(); i2++){
					int v2 = all_idx_u.first[i2];//L(u): u <-> v2
					uint8_t d2 = all_idx_u.second[i2];
					if(v1 != v2 && v1 < v2 && d2 == 1){// generating label of v1 <-> v2, r(v1) > r(v2)
						int vidx, vlb;
						uint8_t dlb = d1 + d2;
						vidx = v1 > v2? v1: v2;
						vlb = v1 < v2? v1: v2;
						std::pair<std::vector<int>, std::vector<uint8_t> >
							&cur_idx_vidx = cur_idx[vidx], &all_idx_vidx = all_idx[vidx];
						int loc = binary_search(all_idx_vidx.first, vlb, 0, all_idx_vidx.first.size());
						if(all_idx_vidx.first.size()==0 || vlb != all_idx_vidx.first[loc]){//avoid duplications in all labels
							loc = binary_search(cur_idx_vidx.first, vlb, 0, cur_idx_vidx.first.size());
							if(cur_idx_vidx.first.size()==0 || vlb != cur_idx_vidx.first[loc]){//avoid duplications in cur labels
								std::vector<int>::iterator fit = cur_idx_vidx.first.begin();
								std::vector<uint8_t>::iterator sit = cur_idx_vidx.second.begin(); 
								cur_idx_vidx.first.insert(fit+loc,vlb);
								cur_idx_vidx.second.insert(sit+loc,dlb);
							}
						}
					}
				}		
				rule1time += GetCurrentTimeSec();

				//Generating by Rule 2
				double rule2time = -GetCurrentTimeSec();
				for(int i2 = 0; i2 < inv_all_idx_u.first.size(); i2++){
					int v2 = inv_all_idx_u.first[i2];//L[v2]: v2 <-> u
					uint8_t d2 = inv_all_idx_u.second[i2];
					if(v1 != v2 && v1 < v2 && d2 == 1){//generating label of v1 <-> v2, r(v1) > r(v2)
						int vidx, vlb;
						uint8_t dlb = d1 + d2;
						vidx = v1 > v2? v1: v2;
						vlb = v1 < v2? v1: v2;
						std::pair<std::vector<int>, std::vector<uint8_t> >
							&cur_idx_vidx = cur_idx[vidx], &all_idx_vidx = all_idx[vidx];
						int loc = binary_search(all_idx_vidx.first, vlb, 0, all_idx_vidx.first.size());
						if(all_idx_vidx.first.size()==0 || vlb != all_idx_vidx.first[loc]){//avoid duplications in all labels
							loc = binary_search(cur_idx_vidx.first, vlb, 0, cur_idx_vidx.first.size());
							if(cur_idx_vidx.first.size()==0 || vlb != cur_idx_vidx.first[loc]){//avoid duplications in cur labels
								std::vector<int>::iterator fit = cur_idx_vidx.first.begin();
								std::vector<uint8_t>::iterator sit = cur_idx_vidx.second.begin(); 
								cur_idx_vidx.first.insert(fit+loc,vlb);
								cur_idx_vidx.second.insert(sit+loc,dlb);
							}
						}
					}
				}	
				rule2time += GetCurrentTimeSec();
				//std::cout<<"1:"<<rule1time/(rule1time+rule2time)<<std::endl;
				//std::cout<<"2:"<<rule2time/(rule1time+rule2time)<<std::endl;
			}
		}/*}}}*/
		std::cout<<"Label Generaton time:"<<gent+GetCurrentTimeSec()<<std::endl;
		//Clearing prev label
		for(int u = 0; u < V; u++){
			prev_idx[u].first.clear();
			prev_idx[u].second.clear();
		}
		
		prunt = -GetCurrentTimeSec();	

		//Label Pruning
		for(int u = 0; u < V; u++){/*{{{*/
			std::pair<std::vector<int>, std::vector<uint8_t> >
				&all_idx_u = all_idx[u], &prev_idx_u = prev_idx[u], &cur_idx_u = cur_idx[u];
			for(int i = 0; i < cur_idx_u.first.size(); i++){
				int tv = cur_idx_u.first[i];		
				int td = cur_idx_u.second[i];
				
				//std::pair<std::vector<int>, std::vector<uint8_t> >
				//	&all_idx_v = all_idx[i];
				std::pair<std::vector<int>, std::vector<uint8_t> >
					&all_idx_v = all_idx[tv];

				int dmin = INT_MAX;
				for(int i1 = 0, i2 = 0; ;){
					if(i1==all_idx_u.first.size() || i2==all_idx_v.first.size() || dmin <= td) break;
					int v1 = all_idx_u.first[i1], v2 = all_idx_v.first[i2];
					if( v1 == v2 ){
						int dd = all_idx_u.second[i1] + all_idx_v.second[i2];
						if(dd < dmin) dmin = dd;
						i1++;
						i2++;
					} else{
						i1 += v1 < v2 ? 1 : 0;
						i2 += v1 > v2 ? 1 : 0;
					}
				}
				if(dmin > td){//not pruned
					//if(dmin!=INT_MAX) std::cout<<"dmin:"<<dmin<<"+td:"<<td<<",INTMAX:"<<INT_MAX<<std::endl;
					prev_idx_u.first.push_back(tv);
					prev_idx_u.second.push_back(td);
					all_idx_u.first.push_back(tv);
					all_idx_u.second.push_back(td);
					pnum++;
				}
			}
			cur_idx_u.first.clear();
			cur_idx_u.second.clear();
		}/*}}}*/
		std::cout<<"Label Pruning time:"<<prunt+GetCurrentTimeSec()<<std::endl;
		std::cout<<"Time in this iteration:"<<itert+GetCurrentTimeSec()<<std::endl;
		std::cout<<"prev num:"<<pnum<<std::endl;
		//for(int v = 0; v < V; v++){
		//	std::cout<<v<<":";
		//	for(int u = 0; u < all_idx[v].first.size();u++){
		//		std::cout<<all_idx[v].first[u]<<",";
		//	}
		//	std::cout<<std::endl;
		//}
		std::vector<int> thisItInfo(V);
		for(int v = 0 ; v < V; v++){
			thisItInfo[v] = prev_idx[v].first.size();
		}
		itInfo.push_back(thisItInfo);
	}	
	/*
	char* itName = std::strcat(dName, ".info");
	std::ofstream ofs(itName);
	for(int u = 0; u < V; u++){
		for(int v = 0 ; v < itInfo.size(); v++){
			ofs << itInfo[v][u] << ",";
		}
		ofs << std::endl;
	}
	ofs.close();
	*/
	for(int v = 0 ; v < V; v++){
		/*
		int k = all_idx[v].first.size() + 2;
		index_[inv[v]].spt_v = (uint32_t*) malloc( k * sizeof(uint32_t));
		index_[inv[v]].spt_d = (uint8_t*) malloc( k * sizeof(uint8_t));
		if(!index_[inv[v]].spt_v || !index_[inv[v]].spt_d ){
			Free();
			return false;
		}
		

		//We restore the ranking id here
		for(int i = 0; i < k - 2; i++) index_[inv[v]].spt_v[i] = all_idx[v].first[i];
		for(int i = 0; i < k - 2; i++) index_[inv[v]].spt_d[i] = all_idx[v].second[i];
		index_[inv[v]].spt_v[k-2] = inv[v];
		index_[inv[v]].spt_d[k-2] = 0;
		index_[inv[v]].spt_v[k-1] = V;
		index_[inv[v]].spt_d[k-1] = INF8;
		
		all_idx[v].first.clear();
		all_idx[v].second.clear();
		prev_idx[v].first.clear();
		prev_idx[v].second.clear();
		*/
		
		//We keep the ranking as the vertex id here
		int k = all_idx[v].first.size() + 2;
		index_[v].spt_v = (uint32_t*) malloc( k * sizeof(uint32_t));
		index_[v].spt_d = (uint8_t*) malloc( k * sizeof(uint8_t));
		if(!index_[v].spt_v || !index_[v].spt_d ){
			Free();
			return false;
		}
		for(int i = 0; i < k - 2; i++) index_[v].spt_v[i] = all_idx[v].first[i];
		for(int i = 0; i < k - 2; i++) index_[v].spt_d[i] = all_idx[v].second[i];
		index_[v].spt_v[k-2] = v;
		index_[v].spt_d[k-2] = 0;
		index_[v].spt_v[k-1] = V;
		index_[v].spt_d[k-1] = INF8;
		
		all_idx[v].first.clear();
		all_idx[v].second.clear();
		prev_idx[v].first.clear();
		prev_idx[v].second.clear();
	}

	return true;
}

template<int kNumBitParallelRoots>
void BNLabel<kNumBitParallelRoots>
::quick_sort(std::pair<std::vector<int>, std::vector<uint8_t> > &pair_array, int first, int last){
	if(first < 0 || last < 0 || first >= pair_array.first.size() || first >= pair_array.second.size() || 
			pair_array.first.size() != pair_array.second.size() ) return;
	if(first != last){
		int left = first;
		int right = last;
		int pivot = left++;
		while(left != right){
			if(pair_array.first[left] < pair_array.first[pivot]){
				++left;
			} else {
				while(left != right && pair_array.first[pivot] < pair_array.first[right]) --right;
				int tmp_first = pair_array.first[left];
				uint8_t tmp_second = pair_array.second[left];
				pair_array.first[left] = pair_array.first[right];
				pair_array.second[left] = pair_array.second[right];
				pair_array.first[right] = tmp_first;
				pair_array.second[right] = tmp_second;
			}
		}

		--left;
		int tmp_first = pair_array.first[left];
		uint8_t tmp_second = pair_array.second[left];
		pair_array.first[left] = pair_array.first[pivot];
		pair_array.second[left] = pair_array.second[pivot];
		pair_array.first[pivot] = tmp_first;
		pair_array.second[pivot] = tmp_second;

		quick_sort(pair_array, first, left);
		quick_sort(pair_array, right, last);
	}	
}

template<int kNumBitParallelRoots>
int BNLabel<kNumBitParallelRoots>
::binary_search(const std::vector<int> &array, const int &key, const int &low, const int & high){
	int imax = high-1, imin = low;
	while(imax >= imin){
		int imid = imin + (imax-imin)/2;
		if(array[imid] == key)
			return imid;
		else if (array[imid] < key)
			imin = imid + 1;
		else
			imax = imid - 1;
	}
	return imin;		
}
#endif //PLLABELING_H
