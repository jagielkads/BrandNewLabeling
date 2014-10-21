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
	bool ConstructPPLIndex(const std::vector<std::vector<int> > &adj, const std::vector<int> &inv);
	bool ConstructHDIndex(const std::vector<std::vector<int> > &adj, const std::vector<int> &inv);
	bool DegreeOrdering(std::vector<std::vector<int> > &adj, std::vector<int> &inv);

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

	private:	
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
::ConstructPPLIndex(const std::vector<std::vector<int> > &adj, const std::vector<int> &inv){
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
::StoreInfo(const char *filename){
	std::ofstream ofs(filename);
	if(ofs == false)
		return false;
	if(index_==NULL)
		return false;
	int &V = num_v_;
	for(int v = 0; v < V; v++){
		const index_t &idx = index_[v];
		int counter = 0;
		for(int i = 0 ; ; i++){
			if(idx.spt_v[i] == V) break;
			
			counter++;
		}
		ofs << counter << std::endl;
	}
	ofs.close();
	return true;
}

template<int kNumBitParallelRoots>
bool BNLabel<kNumBitParallelRoots>
::StoreIndex(const char *filename){
  std::ofstream ofs(filename);
  uint32_t num_v = num_v_;//, num_bpr = kNumBitParallelRoots;
  ofs.write((const char*)&num_v,   sizeof(num_v));
  //ofs.write((const char*)&num_bpr, sizeof(num_bpr));

  for (int v = 0; v < num_v_; ++v) {
		index_t &idx = index_[v];

		//for (int i = 0; i < kNumBitParallelRoots; ++i) {
		//int8_t d = idx.bpspt_d[i];
		//uint64_t a = idx.bpspt_s[i][0];
		//uint64_t b = idx.bpspt_s[i][1];
		//ofs.write((const char*)&d, sizeof(d));
		//ofs.write((const char*)&a, sizeof(a));
		//ofs.write((const char*)&b, sizeof(b));
		//}

		int32_t s;
		for (s = 1; idx.spt_v[s - 1] != num_v; ++s) continue;  // Find the sentinel
		ofs.write((const char*)&s, sizeof(s));
		for (int i = 0; i < s; ++i) {
		int32_t l = idx.spt_v[i];
		int8_t  d = idx.spt_d[i];
		ofs.write((const char*)&l, sizeof(l));
		ofs.write((const char*)&d, sizeof(d));
		}
  
	}
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
	

	//Initialization: Generating labels based on directed linking	
	pnum = 0;
	for(int v = 0; v < V; v++){
		std::pair<std::vector<int>, std::vector<uint8_t> >
			&all_idx_v = all_idx[v];
		std::pair<std::vector<int>, std::vector<uint8_t> >
			&prev_idx_v = prev_idx[v];
		
		all_idx_v.first.push_back(v);
		all_idx_v.second.push_back(0);
		prev_idx_v.first.push_back(v);
		prev_idx_v.second.push_back(0);
		pnum++;	
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
		///quicksort for all idx
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
				for(int i2 = 0; i2 < all_idx_u.first.size(); i2++){
					int v2 = all_idx_u.first[i2];//L(u): u <-> v2
					uint8_t d2 = all_idx_u.second[i2];
					if(v1 != v2){// generating label of v1 <-> v2
						int vidx, vlb;
						uint8_t dlb = d1 + d2;
						vidx = v1 > v2? v1: v2;
						vlb = v1 < v2? v1: v2;
						std::pair<std::vector<int>, std::vector<uint8_t> >
							&cur_idx_vidx = cur_idx[vidx];
						int loc = binary_search(cur_idx_vidx.first, vlb, 0, cur_idx_vidx.first.size());
						if(cur_idx_vidx.first.size()==0 || vlb != cur_idx_vidx.first[loc]){//not duplicates,generating new label
							std::vector<int>::iterator fit = cur_idx_vidx.first.begin();
							std::vector<uint8_t>::iterator sit = cur_idx_vidx.second.begin(); 
							cur_idx_vidx.first.insert(fit+loc,vlb);
							cur_idx_vidx.second.insert(sit+loc,dlb);
						}
					}
				}		

				//Generating by Rule 2
				for(int i2 = 0; i2 < inv_all_idx_u.first.size(); i2++){
					int v2 = inv_all_idx_u.first[i2];//L[v2]: v2 <-> u
					uint8_t d2 = inv_all_idx_u.second[i2];
					if(v1 != v2){//generating label of v1 <-> v2
						int vidx, vlb;
						uint8_t dlb = d1 + d2;
						vidx = v1 > v2? v1: v2;
						vlb = v1 < v2? v1: v2;
						std::pair<std::vector<int>, std::vector<uint8_t> >
							&cur_idx_vidx = cur_idx[vidx];
						int loc = binary_search(cur_idx_vidx.first, vlb, 0, cur_idx_vidx.first.size());
						if(cur_idx_vidx.first.size()==0 || vlb != cur_idx_vidx.first[loc]){//not duplicates,generating new label
							std::vector<int>::iterator fit = cur_idx_vidx.first.begin();
							std::vector<uint8_t>::iterator sit = cur_idx_vidx.second.begin(); 
							cur_idx_vidx.first.insert(fit+loc,vlb);
							cur_idx_vidx.second.insert(sit+loc,dlb);
						}
					}
				}	
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
				
				std::pair<std::vector<int>, std::vector<uint8_t> >
					&all_idx_v = all_idx[i];
				int dmin = INT_MAX;
				for(int i1 = 0, i2 = 0; ;){
					if(i1==all_idx_u.first.size() || i2==all_idx_v.first.size() || dmin < td) break;
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
