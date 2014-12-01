// 2_hop_cover.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include<vector>
#include<fstream>
#include<utility>
#include<algorithm>
#include<iostream>

int num_v_ = 0;
struct path{
	int s, t;
	int w;
	std::vector<int> mp;

	path(){}
	path(int ss, int tt, int ww, std::vector<int> mpp){
		s = ss; t = tt; w = ww;  mp = mpp;
	}
	bool inpath(int v){
		if (s == v || t == v)
			return true;
		for (int i = 0; i < mp.size(); i++){
			if (mp[i] == v)
				return true;
		}
		return false;
	}
};

bool LoadGraph(char* sName, std::vector<std::pair<int, int> > &es, std::vector<std::vector<int> > &adj){
	/*
	The graph is stored as directed graph for simplicity. For any u,v is connected, there is only u v or v u in the file.
	*/
	std::ifstream ifs(sName, std::ifstream::in);
	if (ifs.bad())
		return false;
	for (int u, v; ifs >> u >> v;) {
		es.push_back(std::make_pair(u, v));
	}

	if (ifs.bad())
		return false;
	ifs.close();

	int &V = num_v_;

	for (int i = 0; i < es.size(); i++)
		V = std::max(V, std::max(es[i].first, es[i].second));
	V++;
	if (adj.size() != V) adj.resize(V);
	for (int i = 0; i < es.size(); i++){
		if (es[i].first != es[i].second){
			//std::cout << es[i].first << ":" << es[i].second << std::endl;
			adj[es[i].first].push_back(es[i].second);
			//adj[es[i].second].push_back(es[i].first);
		}
	}
	return true;
}

bool Dijkstra(int s, const std::vector<std::vector<int> > &adj, std::vector<path> &ssp){
	int V = num_v_;
	if (s < 0 || s >= V)
		return false;
	std::vector<int> d(V);
	//std::vector<int> pv(V);
	std::vector<std::vector<std::pair<int, int> > > pv(V);
	std::vector<int> Q(V);
	int Qsize = 0;

	d[s] = 0;
	for (int v = 0; v < V; v++){
		if (v != s){
			d[v] = INT_MAX;
		}
		Q[v] = 0;
		Qsize++;
	}

	while (Qsize != 0){
		int min_d = INT_MAX;
		int u = -1;

		for (int v = 0; v < Q.size(); v++){
			if (Q[v] == 0){
				if (d[v] < min_d){
					min_d = d[v];
					u = v;
				}
			}
		}

		Q[u] = 1;
		Qsize--;

		for (int c = 0; c < adj[u].size(); c++){
			int v = adj[u][c];
			int alt = d[u] + 1;
			if (alt < d[v]){
				d[v] = alt;
				for (int i = 0; i < pv[v].size(); i++){
					if (pv[v][i].second > d[v])
						pv[v].erase(pv[v].begin() + (i--));
				}
			}
			if (alt == d[v]) pv[v].push_back(std::make_pair(u, alt));
		}
	}

	if (ssp.size() != 0) ssp.clear();
	path self(s, s, 0, std::vector<int>(0));
	ssp.push_back(self);

	std::vector<bool> visited(V, false);
	visited[s] = true;
	int vSize = V - 1;
	while (vSize != 0){
		int min_d = INT_MAX;
		int u = -1;
		for (int v = 0; v < V; v++){
			if (v != s &&visited[v] == false){
				if (min_d > d[v]){
					min_d = d[v];
					u = v;
				}
			}
		}

		visited[u] = true;
		vSize--;
		for (int i = 0; i < pv[u].size(); i++){
			int tpv = pv[u][i].first;
			for (int j = 0; j < ssp.size(); j++){
				if (ssp[j].s == s && ssp[j].t == tpv){
					std::vector<int> tp(ssp[j].mp);
					if (d[u] != 0 && d[u] != 1){
						tp.push_back(tpv);
					}
					path tpath(s, u, d[u], tp);
					ssp.push_back(tpath);
				}
			}
		}

	}

}
int _tmain(int argc, _TCHAR* argv[]){
	std::vector<std::pair<int, int> > es(0);
	std::vector<std::vector<int> > adj(0);
	char* sName = "G:\\prjs\\git\\BrandNewLabeling\\data\\sample.txt";
	LoadGraph(sName, es, adj);

	std::vector<std::vector<path> > ssps(num_v_);
	for (int v = 0; v < num_v_; v++){
		Dijkstra(v, adj, ssps[v]);
	}

	int c0 = 0, c1 = 0, c2 = 0;
	int u0 = 0, u1 = 1, u2 = 2;
	for (int v = 0; v < ssps.size(); v++){
		for (int i = 0; i < ssps[v].size(); i++){
			if (ssps[v][i].inpath(u0)) c0++;
			if (ssps[v][i].inpath(u1)) c1++;
			if (ssps[v][i].inpath(u2)) c2++;
			std::cout << ssps[v][i].s << "->";
			for (int j = 0; j < ssps[v][i].mp.size(); j++){
				std::cout << ssps[v][i].mp[j] << "->";
			}
			std::cout << ssps[v][i].t << std::endl;
		}
	}



	return 0;
}

