#ifndef GEM_COMMON
#define GEM_COMMON

#include <limits.h>
#include <algorithm>
#include <sys/mman.h>

namespace gem
{

typedef int IndexType;
static IndexType INDEX_BOTTOM_BOUND = -1;
static IndexType INDEX_UPPER_BOUND = INT_MAX;

struct __attribute__((packed)) EmptyProperty { int dummy[]; };

template<typename EdgeType>
struct __attribute__((packed)) DirectEdge {
  IndexType ver;
  EdgeType val;
  DirectEdge() {}
  DirectEdge(IndexType x, EdgeType y) : ver(x), val(y) {}
};

template<typename EdgeType>
struct __attribute__((packed)) EdgeTriplet {
  IndexType from, to;
  EdgeType value;

  EdgeTriplet() {}
  EdgeTriplet(IndexType f, IndexType t, EdgeType v)
  : from(f), to(t), value(v) {}
};

template<typename EdgeType>
struct EdgeTripletGreaterByFrom {
  bool operator()(const EdgeTriplet<EdgeType>& lx, const EdgeTriplet<EdgeType>& rx ) const {
    if (lx.from != rx.from) return lx.from < rx.from;
    return lx.to < rx.to;
  }
};


template<typename EdgeType>
struct EdgeTripletLesserByValue {
  bool operator()(const EdgeTriplet<EdgeType>& lx, const EdgeTriplet<EdgeType>& rx ) const {
    return lx.value > rx.value;
  }
};

template<typename EdgeType>
struct EdgeTripletGreaterByValue {
  bool operator()(const EdgeTriplet<EdgeType>& lx, const EdgeTriplet<EdgeType>& rx ) const {
    return lx.value < rx.value;
  }
};

template<typename EdgeType>
struct EdgeTripletPartitionerByValue {
  EdgeType barrier;
  EdgeTripletPartitionerByValue(EdgeType b) : barrier(b) {}
  bool operator()(const EdgeTriplet<EdgeType>& x) const {
    return x.value < barrier;
  }
};

template<typename EdgeType>
struct EdgeTripletGreaterByGrid {
  uint64_t grid_width;

  EdgeTripletGreaterByGrid(uint64_t w) : grid_width(w) {}

  bool operator()(const EdgeTriplet<EdgeType>& lx, const EdgeTriplet<EdgeType>& rx ) const {
    uint32_t lgx = lx.from/grid_width, rgx = rx.from/grid_width;
    if (lgx != rgx) return lgx < rgx;

    uint32_t lgy = lx.to/grid_width, rgy = rx.to/grid_width;
    if (lgy != rgy) return lgx < rgx;

    if (lx.from != rx.from) return lx.from < rx.from;
    return lx.to < rx.to;
  }
};



template<typename EdgeType>
struct Cell {
  IndexType start_vertex, end_vertex, start_dst_vertex, end_dst_vertex;
  EdgeType min_edge;
  size_t start_edge_index, end_edge_index;
};

template<typename EdgeType>
struct AdjCell {
  IndexType start_vertex, end_vertex, start_dst_vertex, end_dst_vertex;
  uint32_t edge_cnt;
  EdgeType min_edge;
  
  uint32_t *index;
  DirectEdge<EdgeType> *edges;

  inline DirectEdge<EdgeType> *EdgeBegin(IndexType u) { return edges + index[u - start_vertex]; }
  inline DirectEdge<EdgeType> *EdgeEnd(IndexType u) { return edges + index[1 + u - start_vertex]; }
};

} // namespace gem
#endif // gem_COMMON
