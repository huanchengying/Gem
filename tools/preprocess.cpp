#include "gem/all.h"

#include <stdio.h>
#include <unordered_map>
#include <queue>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;
using namespace gem;

typedef int VertexType;
typedef int EdgeType;

typedef int EdgeType;
int main(int argc, char ** argv) {
  FILE *fp = fopen (argv[1], "r");
  cerr << "Input:\t" << argv[1] << "\n";
  string output_prefix(argv[2]);

  // *******************
  //  Vertex
  // *******************
  uint32_t tot_vertex_cnt;
  if (fread(&tot_vertex_cnt , sizeof(tot_vertex_cnt), 1, fp) != 1) {
    cerr << "Read tot_vertex_cnt Error!\n";
    exit(1);
  }
  cerr << "Vertex Cnt:\t" << tot_vertex_cnt << endl;

  // string vertex_file_path(output_prefix + ".vertex");
  // cerr << "Vertex File: " << vertex_file_path << "\n";
  // size_t vertices_size = sizeof(VertexType) * tot_vertex_cnt;
  // ofstream vertex_file(vertex_file_path.c_str(), ios::out | ios::binary);
  // char zero = '\0';
  // for (size_t i = 0; i < vertices_size; ++i) vertex_file.write(&zero, sizeof(zero));
  // vertex_file.close();

  // *******************
  //  Edge
  // *******************
  uint32_t tot_edge_cnt;
  if (fread(&tot_edge_cnt , sizeof(tot_edge_cnt), 1, fp) != 1) {
    cerr << "Read tot_edge_cnt Error!\n";
    exit(1);
  }
  cerr << "Edge Cnt:\t" << tot_edge_cnt << endl;

  vector<EdgeTriplet<EdgeType> > edges(tot_edge_cnt);
  size_t edge_size = sizeof(IndexType) * 2 + sizeof(EdgeType);
  size_t edges_size = edge_size * tot_edge_cnt;
  if (fread((void *)(&edges[0]) , sizeof(uint8_t), edges_size, fp) != edges_size) {
    cerr << "Read Edge Error!\n";
    exit(1);
  }
  sort(edges.begin(), edges.end(), EdgeTripletGreaterByFrom<EdgeType>());

  // string edge_file_path(output_prefix + ".edge");
  // cerr << "Edge File: " << edge_file_path << "\n";
  // ofstream edge_file(edge_file_path.c_str(), ios::out | ios::binary);
  // edge_file.write((char *)&edges[0], edges_size);
  // edge_file.close();

  // *******************
  //  Info
  // *******************
  string info_file_path(output_prefix + ".info");
  cerr << "Info File: " << info_file_path << "\n";
  // ofstream info_file(info_file_path.c_str(), ios::out);
  // info_file << tot_vertex_cnt << " -- Tot Vertex Cnt\n";
  // info_file << tot_edge_cnt << " -- Tot Edge Cnt\n";
  // info_file << sizeof(VertexType) << " -- Size of Vertex\n";
  // info_file << sizeof(EdgeType) << " -- Size of Edge\n";
  // info_file.close();


  // *******************
  //  Index
  // *******************
  vector<size_t> index{0};
  IndexType cur = 0;
  for (size_t i = 0; i < edges.size(); ++i) {
    auto& e = edges[i];
    while (cur < e.from) {
      ++cur;
      index.push_back(i);
    }
  }
  for (; cur < tot_vertex_cnt; ++cur) {
    index.push_back(edges.size());
  }
  cerr << index.size() << "==\n";

  string index_file_path(output_prefix + ".index");
  cerr << "Index File: " << index_file_path << "\n";
  ofstream index_file(index_file_path.c_str(), ios::out | ios::binary);
  index_file.write((char *)&index[0], index.size() * sizeof(size_t));
  index_file.close();
}
