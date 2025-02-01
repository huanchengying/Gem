#ifndef GEM_GEM
#define GEM_GEM

#include "common.h"
#include "bit_map.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <fstream>
#include <string>

using namespace std;

#define CHAR_BUF_SIZE 256

namespace gem
{

template<typename VertexType, typename EdgeType>
class gem {
public:
  gem(char *input_prefix);
  ~gem();

  uint64_t TotVertexCnt() { return _tot_vertex_cnt; }
  uint64_t TotEdgeCnt() { return _tot_edge_cnt; }
  inline VertexType& Vertex(IndexType u) { return vertices[u]; }
  inline EdgeTriplet<EdgeType>* EdgePtr(size_t i) { return edges + i; }
#ifdef USE_INDEX
  inline size_t Index(IndexType u) { return indexes[u]; }
#endif
  inline void InitVertex(VertexType x) {
    std::fill(vertices, vertices + _tot_vertex_cnt, x);
  }

#ifdef USE_ACTIVITY
  BitMap *in_activity, *out_activity;
#endif

private:
  void analysisInfoFile(string file_name);

  uint64_t _tot_vertex_cnt, _tot_edge_cnt;
  char *mmap_edges, *mmap_vertices;
  VertexType *vertices;
  EdgeTriplet<EdgeType> *edges;

#ifdef USE_INDEX
  char *mmap_index;
  size_t *indexes;
#endif
};

template<typename VertexType, typename EdgeType>
gem<VertexType, EdgeType>::gem(char *input_prefix) {
  string prefix(input_prefix);
  analysisInfoFile(prefix + ".info");

  {
    string vertex_file_name(prefix + ".vertex");
    int vertex_fd = open(vertex_file_name.c_str(), O_RDWR);
    mmap_vertices = (char *)mmap(0, sizeof(VertexType) * _tot_vertex_cnt, PROT_WRITE, MAP_PRIVATE, vertex_fd, 0);
    vertices = (VertexType *)mmap_vertices;
    close(vertex_fd);
    madvise(mmap_vertices, sizeof(VertexType) * _tot_vertex_cnt, MADV_WILLNEED);
  }

  {
    string edge_file_name(prefix + ".edge");
    int edge_fd = open(edge_file_name.c_str(), O_RDONLY);
    size_t edges_size = 2 * sizeof(IndexType) + sizeof(EdgeType);
    mmap_edges = (char *)mmap(0, edges_size * _tot_edge_cnt, PROT_READ, MAP_PRIVATE, edge_fd, 0);
    edges = (EdgeTriplet<EdgeType> *)mmap_edges;
    close(edge_fd);
  }

#ifdef USE_ACTIVITY
  in_activity = new BitMap(_tot_vertex_cnt);
  out_activity = new BitMap(_tot_vertex_cnt);
#endif

#ifdef USE_INDEX
  {
    string index_file_name(prefix + ".index");
    int index_fd = open(index_file_name.c_str(), O_RDONLY);
    mmap_index = (char *)mmap(0, sizeof(size_t) * _tot_vertex_cnt, PROT_READ, MAP_PRIVATE, index_fd, 0);
    indexes = (size_t *)mmap_index;
    close(index_fd);
    madvise(mmap_index, sizeof(size_t) * _tot_vertex_cnt, MADV_WILLNEED);
  }
#endif
}

template<typename VertexType, typename EdgeType>
gem<VertexType, EdgeType>::~gem() {
  munmap(mmap_vertices, sizeof(VertexType) * _tot_vertex_cnt);

  size_t edges_size = 2 * sizeof(IndexType) + sizeof(EdgeType);
  munmap(mmap_edges, edges_size * _tot_edge_cnt);

#ifdef USE_ACTIVITY
  free(in_activity);
  free(out_activity);
#endif

#ifdef USE_INDEX
  munmap(mmap_index, sizeof(size_t) * _tot_vertex_cnt);
#endif
}

template<typename VertexType, typename EdgeType>
void gem<VertexType, EdgeType>::analysisInfoFile(string file_name) {
  char buf[CHAR_BUF_SIZE];
  ifstream f(file_name.c_str());
  
  f >> _tot_vertex_cnt; f.getline(buf, CHAR_BUF_SIZE);
  cerr << "Vertex Cnt: " << _tot_vertex_cnt << "\n";

  f >> _tot_edge_cnt; f.getline(buf, CHAR_BUF_SIZE);
  cerr << "Edge Cnt: " << _tot_edge_cnt << "\n";

  size_t sz;
  f >> sz; f.getline(buf, CHAR_BUF_SIZE);
  if (sz != sizeof(VertexType)) {
    cerr << "[ERROR] Vertex Size Mismatch! " << sz << "\n";
    exit(1);
  }
  f >> sz; f.getline(buf, CHAR_BUF_SIZE);
  if (sz != sizeof(EdgeType)) {
    cerr << "[ERROR] Edge Size Mismatch! " << sz << "\n";
    exit(1);
  }

  f.close();
}


} // namespace gem
#endif // gem_gem
