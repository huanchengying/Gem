#include "gem/all.h"
#include <list>

#include <stdio.h>
#include <math.h>
#include <unordered_map>
#include <queue>
#include <string>
#include <iostream>
#include <fstream>
#include <set>

using namespace std;
using namespace gem;

typedef float VertexType;
typedef float EdgeType;

struct AdjCellP {
  IndexType start_vertex, end_vertex, start_dst_vertex, end_dst_vertex;
  uint32_t edge_cnt;
  EdgeType min_edge;
  vector<uint32_t> index;
  vector<DirectEdge<EdgeType> > edges;

  void FromEdgeTriplets(vector<EdgeTriplet<EdgeType> >& vec, size_t start, size_t cnt) {
    index.clear(); edges.clear();
    if (cnt == 0) { start_vertex = end_vertex = start_dst_vertex = end_dst_vertex = edge_cnt = min_edge = 0; return;}
    sort(vec.begin() + start, vec.begin() + start + cnt, EdgeTripletGreaterByFrom<EdgeType>());
    start_vertex = vec[start].from;
    end_vertex = vec[start + cnt - 1].from;
    IndexType vertex_cnt = 1 + end_vertex - start_vertex;

    start_dst_vertex = end_dst_vertex = vec[start].to;
    min_edge = vec[start].value;
    index.push_back(0);
    IndexType cur_vertex = start_vertex;
    for (size_t i = start; i < start + cnt; ++i) {
      while (cur_vertex < vec[i].from) {
        index.push_back(i - start);
        ++cur_vertex;
      }
      start_dst_vertex = min(start_dst_vertex, vec[i].to);
      end_dst_vertex = max(end_dst_vertex, vec[i].to);
      min_edge = min(min_edge, vec[i].value);
    }
    index.push_back(cnt);
    if (cur_vertex != end_vertex || index.size() != vertex_cnt + 1) {
      cerr << "[ERROR] FromEdgeTriplets Index\n";
    }

    edge_cnt = cnt;
    if (edge_cnt != cnt) {
      cerr << "[ERROR] FromEdgeTriplets Edge Cnt\n";
    }
    for (size_t i = start; i < start + cnt; ++i) {
      IndexType ver = vec[i].to;
      EdgeType val = vec[i].value;
      edges.emplace_back(ver, val);
    }
  }
};

double *pagerank;
double *deg;
double get_pagerank(EdgeTriplet<EdgeType> a){
    return min(pagerank[a.from],pagerank[a.to])/(a.value);
  //  return min(deg[a.from],deg[a.to])/(a.value);
}
bool cmp1 (EdgeTriplet<EdgeType> a, EdgeTriplet<EdgeType> b){
    return get_pagerank(a)<get_pagerank(b);
}
set<pair<int,double>>ss;
typedef pair<int,double> ID;
ID* gg[2];
int binary_search(ID* nn,int num, double p){
    int le=0,ri=num-1;
    while(le<=ri){
        int mid=(le+ri)/2;
        auto now=nn+mid;
        if(now->second<=p) le=mid+1;
        else ri=mid-1;
    }
    return le;
}
int cmp2(pair<int,double>a, pair<int,double>b){
    return a.second <b.second;
}
int main(int argc, char ** argv) {
  FILE *fp = fopen (argv[1], "r");
  cerr << "Input:\t" << argv[1] << "\n";
  string output_prefix(argv[2]);
  
  uint64_t grid_size = atoll(argv[3]);
  IndexType grid_width = atoi(argv[4]);
  uint64_t abstract_size = atoll(argv[5]);

  // *******************
  //  Vertex
  // *******************
  uint64_t tot_vertex_cnt;
  if (fread(&tot_vertex_cnt , sizeof(tot_vertex_cnt), 1, fp) != 1) {
    cerr << "Read tot_vertex_cnt Error!\n";
    exit(1);
  }
  cerr << "Vertex Cnt:\t" << tot_vertex_cnt << endl;

  string vertex_file_path(output_prefix + "-W" + std::to_string(grid_width) + ".vertex");
  cerr << "Vertex File: " << vertex_file_path << "\n";
  size_t vertices_size = sizeof(VertexType) * tot_vertex_cnt;
  ofstream vertex_file(vertex_file_path.c_str(), ios::out | ios::binary);
  char zero = '\0';
  for (size_t i = 0; i < vertices_size; ++i) vertex_file.write(&zero, sizeof(zero));
  vertex_file.close();

  // *******************
  //  Edge
  // *******************
  string info_file_path(output_prefix + "-W" + std::to_string(grid_width) + ".rgraph_info");
  cerr << "Info File: " << info_file_path << "\n";
  ofstream info_file(info_file_path.c_str(), ios::out);
  info_file << tot_vertex_cnt << " -- Tot Vertex Cnt\n";
  

  uint64_t tot_edge_cnt;
  if (fread(&tot_edge_cnt , sizeof(tot_edge_cnt), 1, fp) != 1) {
    cerr << "Read tot_edge_cnt Error!\n";
    exit(1);
  }
  cerr << "Edge Cnt:\t" << tot_edge_cnt << endl;
  info_file << tot_edge_cnt << " -- Tot Edge Cnt\n";
  
  info_file << sizeof(VertexType) << " -- Size of Vertex\n";
  info_file << sizeof(EdgeType) << " -- Size of Edge\n";

  vector<EdgeTriplet<EdgeType> > edges(tot_edge_cnt);
  size_t edge_size = sizeof(IndexType) * 2 + sizeof(EdgeType);
  size_t edges_size = edge_size * tot_edge_cnt;
  if (fread((void *)(&edges[0]) , sizeof(uint8_t), edges_size, fp) != edges_size) {
    cerr << "Read Edge Error!\n";
    exit(1);
  }
  for(uint64_t i=0;i<tot_edge_cnt;i++) edges[i].value=1.0*i/10+1;
  pagerank=(double *)malloc((tot_vertex_cnt+10)*sizeof(double));
  deg=(double *)malloc((tot_vertex_cnt+10)*sizeof(double));
  //for(int i=0;i<edges.size();i++) deg[edges[i].from]+=edges[i].value;
  for(int i=0;i<edges.size();i++) deg[edges[i].from]++;
  auto cur_t=chrono::steady_clock::now();
  for(int i=0;i<10;i++){
    #pragma omp parallel num_threads(16) 
    for(int j=0;j<edges.size();j++){
        pagerank[edges[j].to]+=pagerank[edges[j].from]/deg[edges[j].from];
    }
  }
  vector<EdgeTriplet<EdgeType> > abstract(abstract_size);
  int no_num=0;
  int id[20];
  pair<int,int> *id_;
  ID *tmp;
  id_=(pair<int,int> *)malloc((abstract_size+10)*sizeof(pair<int,int>));
  tmp=(ID *)malloc((abstract_size+10)*sizeof(ID));
  gg[0]=(ID *)malloc((abstract_size+10)*sizeof(ID));
  gg[1]=(ID *)malloc((abstract_size+10)*sizeof(ID));
  vector<vector<ID>> gk(abstract_size+10);
  int cur_gg=0;
  for(uint64_t i=0;i<tot_edge_cnt;i+=abstract_size){
    if(i==0){
       for(int j=0;j<abstract_size;j++){
            tmp[j]=make_pair(i+j,get_pagerank(edges[i+j]));
       } 
       sort(tmp,tmp+abstract_size,cmp2);
       for(int j=0;j<abstract_size;j++) gg[cur_gg][j]=tmp[j];
       continue;
    }
    for(int j=0;j<20;j++) cerr<<gg[0][j].second<<" ";

    int num=min(abstract_size,(tot_edge_cnt-i));
    {
        #pragma omp parallel num_threads(16)
        for(int j=0;j<num;j++){
            id[j]=binary_search(gg[cur_gg],abstract_size,get_pagerank(edges[i+j]));
       //     cerr<<id[j]<<"\n";
        }
    }
    cerr<<"FUCK\n";
    for(int j=0;j<num;j++){
        gk[id[j]].push_back(make_pair(i+j,get_pagerank(edges[i+j])));
    }
    no_num=0;
    cerr<<"FUCK\n";
    for(int j=0;j<abstract_size;j++){
       sort(gk[j].begin(),gk[j].end(),cmp2);
       for(auto it=gk[j].begin();it<gk[j].end();it++){
            gg[cur_gg^1][no_num++]=(*it);
            if(no_num>=abstract_size) break;
       }
       gk[j].clear();
       if(no_num<abstract_size){
            gg[cur_gg^1][no_num++]=gg[cur_gg][j];
       }
       else break;
    }
    cur_gg=cur_gg^1;
  }
  cerr<<"Preprocess time: "<< chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - cur_t).count()/1000.0<<"\    n";
  /*
  for (uint64_t i = 0; i < abstract_size; ++i) {
    abstract[i] = edges[edges.size() - 1 - i];
  }
  edges.resize(edges.size() - abstract_size);
  cerr<<"Preprocessing time: "<<chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - cur_t).count()/1000.0<<"\n";
  vector<uint64_t> barrier{0};
  while (barrier.back() < tot_vertex_cnt) barrier.push_back(min(barrier.back() + grid_width, tot_vertex_cnt));
  cerr << "Grid Size: " << grid_size << "\n";
  info_file << grid_size << " -- Grid Size\n";
  cerr << "Grid Width: " << grid_width << "\n";
  info_file << grid_width << " -- Grid Width\n";

  cerr << "Barrier\n";
  for (auto b: barrier) cerr << "\t" << b;
  cerr << "\n";
  // info_file << barrier.size() << " -- Barrier Size\n";
  // for (auto b: barrier) info_file << "\t" << b;
  // info_file << "\n";

  vector<vector<vector<EdgeTriplet<EdgeType> > > > edge_grid
    = vector<vector<vector<EdgeTriplet<EdgeType> > > >(
        barrier.size() - 1,
        vector<vector<EdgeTriplet<EdgeType> > >(barrier.size() - 1));
  for (auto e: edges) {
    edge_grid[e.from/grid_width][e.to/grid_width].push_back(e);
  }
  cerr << "Grid Size\n";
  for (int i = 0; i < barrier.size() - 1; ++i) {
    for (int j = 0; j < barrier.size() - 1; ++j) {
      cerr << "\t" << edge_grid[i][j].size();//descriptor.GridCnt(i, j);
    }
    cerr << endl;
  }

  vector<tuple<EdgeType, int, int, size_t, size_t> > cells;
  for (int i = 0; i < barrier.size() - 1; ++i) {
    for (int j = 0; j < barrier.size() - 1; ++j) {
      if (edge_grid[i][j].size() == 0) continue;
      sort(edge_grid[i][j].begin(), edge_grid[i][j].end(), EdgeTripletLesserByValue<EdgeType>());

      size_t step = edge_grid[i][j].size();
      if (step > grid_size) {
        int cnt = edge_grid[i][j].size()/grid_size + (edge_grid[i][j].size() % grid_size == 0 ? 0 : 1);
        step = edge_grid[i][j].size()/cnt + (edge_grid[i][j].size() % cnt == 0 ? 0 : 1);
      }
      for (size_t start = 0; start < edge_grid[i][j].size(); start += step) {
        size_t end = min(edge_grid[i][j].size(), start + step);
        EdgeType min_velue = edge_grid[i][j][end-1].value;
        cells.emplace_back(min_velue, i, j, start, end);
      }
    }
  }
  sort(cells.begin(), cells.end());
  

  string edge_file_path(output_prefix + "-W" + std::to_string(grid_width) + ".rgraph_edge");
  cerr << "Edge File: " << edge_file_path << "\n";
  ofstream edge_file(edge_file_path.c_str(), ios::out | ios::binary);

  uint64_t cursor = 0;
  AdjCellP adj_cell;
  
  // Reverse Abstract
  auto rabstract = abstract;
  for (auto& e: rabstract) {
    auto f = e.from;
    e.from = e.to;
    e.to = f;
  }
  adj_cell.FromEdgeTriplets(rabstract, 0, rabstract.size());
  cerr 
      << "\t" << adj_cell.start_vertex << "\t" << adj_cell.end_vertex 
      << "\t" << adj_cell.start_dst_vertex << "\t" << adj_cell.end_dst_vertex
      << "\t" << adj_cell.edge_cnt << "\t" << adj_cell.min_edge << "\t" << cursor;
  info_file 
      << "\t" << adj_cell.start_vertex << "\t" << adj_cell.end_vertex 
      << "\t" << adj_cell.start_dst_vertex << "\t" << adj_cell.end_dst_vertex
      << "\t" << adj_cell.edge_cnt << "\t" << adj_cell.min_edge << "\t" << cursor;
  edge_file.write((char *)&adj_cell.index[0], sizeof(adj_cell.index[0]) * adj_cell.index.size());
  cursor += sizeof(adj_cell.index[0]) * adj_cell.index.size();
  edge_file.write((char *)&adj_cell.edges[0], sizeof(adj_cell.edges[0]) * adj_cell.edges.size());
  cursor += sizeof(adj_cell.edges[0]) * adj_cell.edges.size();
  cerr << "\t" << cursor << "\n";
  info_file << "\t" << cursor << "\n";

  cerr << "Cell\n";
  info_file << cells.size() + 1 << " -- Cell Cnt\n";

  // Abstract
  adj_cell.FromEdgeTriplets(abstract, 0, abstract.size());
  cerr 
      << "\t" << adj_cell.start_vertex << "\t" << adj_cell.end_vertex 
      << "\t" << adj_cell.start_dst_vertex << "\t" << adj_cell.end_dst_vertex
      << "\t" << adj_cell.edge_cnt << "\t" << adj_cell.min_edge << "\t" << cursor;
  info_file 
      << "\t" << adj_cell.start_vertex << "\t" << adj_cell.end_vertex 
      << "\t" << adj_cell.start_dst_vertex << "\t" << adj_cell.end_dst_vertex
      << "\t" << adj_cell.edge_cnt << "\t" << adj_cell.min_edge << "\t" << cursor;
  edge_file.write((char *)&adj_cell.index[0], sizeof(adj_cell.index[0]) * adj_cell.index.size());
  cursor += sizeof(adj_cell.index[0]) * adj_cell.index.size();
  edge_file.write((char *)&adj_cell.edges[0], sizeof(adj_cell.edges[0]) * adj_cell.edges.size());
  cursor += sizeof(adj_cell.edges[0]) * adj_cell.edges.size();
  cerr << "\t" << cursor << "\n";
  info_file << "\t" << cursor << "\n";

  for (auto& c: cells) {
    EdgeType min_value; int i, j; size_t start, end;
    std::tie(min_value, i, j, start, end) = c;

    // One Cell
    adj_cell.FromEdgeTriplets(edge_grid[i][j], start, end - start);
    cerr 
      << "\t" << adj_cell.start_vertex << "\t" << adj_cell.end_vertex 
      << "\t" << adj_cell.start_dst_vertex << "\t" << adj_cell.end_dst_vertex
      << "\t" << adj_cell.edge_cnt << "\t" << adj_cell.min_edge << "\t" << cursor;
    info_file 
      << "\t" << adj_cell.start_vertex << "\t" << adj_cell.end_vertex 
      << "\t" << adj_cell.start_dst_vertex << "\t" << adj_cell.end_dst_vertex
      << "\t" << adj_cell.edge_cnt << "\t" << adj_cell.min_edge << "\t" << cursor;
    edge_file.write((char *)&adj_cell.index[0], sizeof(adj_cell.index[0]) * adj_cell.index.size());
    cursor += sizeof(adj_cell.index[0]) * adj_cell.index.size();
    edge_file.write((char *)&adj_cell.edges[0], sizeof(adj_cell.edges[0]) * adj_cell.edges.size());
    cursor += sizeof(adj_cell.edges[0]) * adj_cell.edges.size();
    cerr << "\t" << cursor << "\n";
    info_file << "\t" << cursor << "\n";
  }
  */
//  edge_file.close();
 // info_file.close();
}
