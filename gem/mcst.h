#ifndef GEM_MCST_GEM
#define GEM_MCST_GEM

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
class Mcstgem {
public:
  Mcstgem(char *input_prefix);
  ~Mcstgem();

  uint64_t TotVertexCnt() { return _tot_vertex_cnt; }
  uint64_t TotEdgeCnt() { return _tot_edge_cnt; }

  inline IndexType  boot(IndexType u){
	if(root[u]==u) return u;
	return root[u]=boot(root[u]);
  }
  inline void Initroot(IndexType cnt){
	for(int i=0;i<=cnt;i++){
		root[i]=i;
	}
  }
  inline void Setroot(IndexType u,IndexType x){
	root[u]=x;
  }
  BitMap *in_activity, *out_activity;
  std::vector<AdjCell<EdgeType> > cells;

  AdjCell<EdgeType> rabstract;
  std::vector<VertexType> min_expect;
private:
  void analysisInfoFile(string file_name);

  uint64_t _tot_vertex_cnt, _tot_edge_cnt, mmap_edges_size;
  char *mmap_edges, *mmap_vertices;
  IndexType *root;
};

template<typename VertexType, typename EdgeType>
Mcstgem<VertexType, EdgeType>::Mcstgem(char *input_prefix) {
  string prefix(input_prefix);

  {
    string edge_file_name(prefix + ".rgraph_edge");
    int edge_fd = open(edge_file_name.c_str(), O_RDONLY);
    struct stat sb;
    fstat(edge_fd, &sb);
    mmap_edges_size = sb.st_size;
    mmap_edges = (char *)mmap(0, mmap_edges_size, PROT_READ, MAP_PRIVATE, edge_fd, 0);
    close(edge_fd);
  }

  analysisInfoFile(prefix + ".rgraph_info");

  {
    string vertex_file_name(prefix + ".vertex");
    int vertex_fd = open(vertex_file_name.c_str(), O_RDWR);
    mmap_vertices = (char *)mmap(0, sizeof(IndexType) * _tot_vertex_cnt, PROT_WRITE, MAP_PRIVATE, vertex_fd, 0);
    root = (IndexType *)mmap_vertices;
    close(vertex_fd);
    madvise(mmap_vertices, sizeof(IndexType) * _tot_vertex_cnt, MADV_WILLNEED);
  }

  in_activity = new BitMap(_tot_vertex_cnt);
  out_activity = new BitMap(_tot_vertex_cnt);
}

template<typename VertexType, typename EdgeType>
Mcstgem<VertexType, EdgeType>::~Mcstgem() {
  munmap(mmap_vertices, sizeof(VertexType) * _tot_vertex_cnt);

  size_t edges_size = 2 * sizeof(IndexType) + sizeof(EdgeType);
  munmap(mmap_edges, mmap_edges_size);

  free(in_activity);
  free(out_activity);
}

template<typename VertexType, typename EdgeType>
void Mcstgem<VertexType, EdgeType>::analysisInfoFile(string file_name) {
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

  uint64_t grid_size, grid_width;
  f >> grid_size; f.getline(buf, CHAR_BUF_SIZE);
  cerr << "Grid Size: " << grid_size << "\n";
  f >> grid_width; f.getline(buf, CHAR_BUF_SIZE);
  cerr << "Grid Width: " << grid_width << "\n";

  f >> rabstract.start_vertex >> rabstract.end_vertex >> rabstract.start_dst_vertex >> rabstract.end_dst_vertex >> rabstract.edge_cnt >> rabstract.min_edge;
  uint64_t cursor1, cursor2;
  f >> cursor1 >> cursor2;
  rabstract.index = (uint32_t *)(mmap_edges + cursor1);
  rabstract.edges = (DirectEdge<EdgeType> *)(mmap_edges + cursor1 + sizeof(uint32_t) * (rabstract.end_vertex - rabstract.start_vertex + 2));
  if (cursor2 - cursor1 != sizeof(uint32_t) * (rabstract.end_vertex - rabstract.start_vertex + 2) + sizeof(DirectEdge<EdgeType>) * rabstract.edge_cnt) {
    cerr << "[ERROR] analysisInfoFile cursor\n";
  }

  int cell_cnt;
  f >> cell_cnt; f.getline(buf, CHAR_BUF_SIZE);
  cerr << "Cell Cnt: " << cell_cnt << "\n";
  while (cell_cnt--) {
    AdjCell<EdgeType> c;
    f >> c.start_vertex >> c.end_vertex >> c.start_dst_vertex >> c.end_dst_vertex >> c.edge_cnt >> c.min_edge;
    uint64_t cursor1, cursor2;
    f >> cursor1 >> cursor2;
    c.index = (uint32_t *)(mmap_edges + cursor1);
    c.edges = (DirectEdge<EdgeType> *)(mmap_edges + cursor1 + sizeof(uint32_t) * (c.end_vertex - c.start_vertex + 2));
    if (cursor2 - cursor1 != sizeof(uint32_t) * (c.end_vertex - c.start_vertex + 2) + sizeof(DirectEdge<EdgeType>) * c.edge_cnt) {
      cerr << "[ERROR] analysisInfoFile cursor\n";
    }
    cells.push_back(c);
  }
  f.close();
  min_expect.resize(cells.size());
}


} // namespace gem
#endif // gem_RGRAPH_gem
