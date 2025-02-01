#include "gem/all.h"

#include <stdio.h>
#include <unordered_map>
#include <queue>
#include <numeric>
#include <chrono>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <sys/resource.h>


using namespace std;
using namespace gem;

typedef float EdgeType;
typedef float VertexType;

vector<IndexType> TEST_CASE;
vector<EdgeType> gridgraph_dist;
vector<float> gridgraph_reach, gridgraph_final;

struct WorkerController {
  CVBarrier start_barrier, end_barrier;
  atomic<uint64_t> start_cursor;
  vector<uint64_t> stage_barrier;
  uint64_t end_cursor, task_block;
  int task_type;
  bool last_work;
  VertexType upper_bound;

  WorkerController(int num_threads) : start_barrier(num_threads + 1), end_barrier(num_threads + 1), last_work(false) {}
  inline void StartAllWorkers() { start_barrier.Wait(); }
  inline void WaitAllWorkers() { end_barrier.Wait(); }
  
  inline void StopAllWorkers() {
    task_type = 0; 
    start_cursor = end_cursor = 0;
    task_block = 1;
    StartAllWorkers();
    last_work = true;
    WaitAllWorkers();
  }

  // Task Type 1: Relax Edge
  vector<AdjCell<EdgeType> *> cells;
  bool has_relaxed, this_relaxed;

  inline void InitRelaxEdgeTask(vector<AdjCell<EdgeType> *>& c, uint64_t block) {
    task_type = 1; task_block = block;
    cells = c;

    start_cursor = 0;
    stage_barrier.clear();
    for (auto& ptr : cells) {
      uint64_t last = (stage_barrier.size() == 0 ? 0 : stage_barrier.back());
      stage_barrier.push_back(last + 1 + ptr->end_vertex - ptr->start_vertex);
    }
    end_cursor = stage_barrier.back();
  }

  // Task Type 2: Skip Task
  AdjCell<EdgeType> *cell;
  VertexType bound;
  bool skip;

  inline void InitSkipTask(AdjCell<EdgeType> *c, uint64_t block, VertexType b) {
    task_type = 2; task_block = block; bound = b;
    cell = c;
    start_cursor = 0; end_cursor = 1 + cell->end_vertex - cell->start_vertex;
    skip = true;
  }

  // Task Type 3: Min Vertex
  vector<VertexType> min_vertex;
  inline void InitMinVertexTask(uint64_t n, VertexType init) {
    task_type = 3; task_block = 1;
    start_cursor = 0; end_cursor = n;
    min_vertex.resize(n);
    fill(min_vertex.begin(), min_vertex.end(), init);
  }

  // Task Type 4: RAbstract
  AdjCell<EdgeType> *rabstract;
  inline void InitRAbstractTask(RGraphgem<VertexType, EdgeType>& gem, uint64_t block) {
    task_type = 4; rabstract = &(gem.rabstract); task_block = block;
    start_cursor = 0; end_cursor = 1 + rabstract->end_vertex - rabstract->start_vertex;
  }
  // Task 5:
  inline void InitRCellsTask(RGraphgem<VertexType, EdgeType>& gem) {
    task_type = 5; task_block = 1;
    start_cursor = 1; end_cursor = gem.cells.size();
  }
};

bool intersect(uint64_t ll, uint64_t lr, uint64_t rl, uint64_t rr) {
  if (lr <= rl) return false;
  if (ll >= rr) return false;
  return true;
}

struct Worker {
  int id;
  WorkerController& controller;
  RGraphgem<VertexType, EdgeType>& gem;

  Worker(int id_, WorkerController& controller_, RGraphgem<VertexType, EdgeType>& gem_)
    : id(id_), controller(controller_), gem(gem_) {}

  void operator() () const {
    while (!controller.last_work) {
      controller.start_barrier.Wait();
      switch (controller.task_type) {
      case 0: break; // End Task;
      case 1: RelaxEdgeTask(); break;
      case 2: SkipTask(); break;
      case 3: MinVertexTask(); break;
      case 4: RAbstractTask(); break;
      case 5: RCellsTask(); break;
      }
      controller.end_barrier.Wait();
    }
  }

  void RelaxEdgeTask() const {
    bool local_relaxed = false;
    VertexType local_upper_bound = controller.upper_bound;

    while (true) {
      uint64_t cur = controller.start_cursor.fetch_add(controller.task_block);
      if (cur >= controller.end_cursor) break;

      uint64_t stage = 0, block = controller.task_block;
      while (controller.stage_barrier[stage] < cur) ++stage;

      AdjCell<EdgeType> *c = controller.cells[stage];
      uint64_t barrier = controller.stage_barrier[stage];
      uint64_t low_barrier = (stage == 0 ? 0 : controller.stage_barrier[stage-1]);

      for (uint64_t i = 0; i < block; ++i, ++cur) {
        if (cur == barrier) {
          ++stage;
          if (stage == controller.stage_barrier.size()) break;
          c = controller.cells[stage];
          low_barrier = barrier;
          barrier = controller.stage_barrier[stage];
        }

        IndexType from = c->start_vertex + cur - low_barrier;
        if (!gem.in_activity->get_bit(from)) continue;
        VertexType src_dist = gem.Vertex(from);
        if (src_dist >= local_upper_bound) continue;

        for (auto ptr = c->EdgeBegin(from), end = c->EdgeEnd(from); ptr < end; ++ ptr) {
          VertexType new_dist = src_dist + ptr->val;
          VertexType& dst_dist = gem.Vertex(ptr->ver);
          if (dst_dist > new_dist) {
            local_relaxed = true;
            dst_dist = new_dist;
            
            gem.in_activity->set_bit(ptr->ver);
            gem.out_activity->set_bit(ptr->ver);
          }
        }
      }
    }
    if (local_relaxed) {
      controller.this_relaxed = true;
      controller.has_relaxed = true;
    }
  }

  void RAbstractTask() const {
    bool local_relaxed = false;
    auto& rabstract = gem.rabstract;
    while (true) {
      uint64_t cur = controller.start_cursor.fetch_add(controller.task_block);
      if (cur >= controller.end_cursor) break;

      IndexType from = rabstract.start_vertex + cur;
      for (IndexType 
        from = rabstract.start_vertex + cur, 
        end = min(from + (IndexType)controller.task_block, rabstract.end_vertex); 
        from < end; ++from) {

        VertexType src_dist = gem.Vertex(from);
        for (auto ptr = rabstract.EdgeBegin(from), end = rabstract.EdgeEnd(from); ptr < end; ++ ptr) {
          VertexType new_dist = src_dist + ptr->val;
          VertexType& dst_dist = gem.Vertex(ptr->ver);
          if (dst_dist > new_dist) {
            local_relaxed = true;
            dst_dist = new_dist;
          }
        }
        for (int i = 1; i < gem.cells.size(); ++i) {
          auto& c = gem.cells[i];
          if (from < c.start_dst_vertex || from >= c.end_dst_vertex) continue;
          VertexType new_dist = src_dist + c.min_edge;
          VertexType& dst_dist = gem.min_expect[i];
          if (dst_dist > new_dist) {
            local_relaxed = true;
            dst_dist = new_dist;
          }
        }
      }
    }
    if (local_relaxed) {
      controller.has_relaxed = true;
    }
  }

  void RCellsTask() const {
    bool local_relaxed = false;
    auto& rabstract = gem.rabstract;
    while (true) {
      uint64_t cur = controller.start_cursor.fetch_add(controller.task_block);
      if (cur >= controller.end_cursor) break;

      auto& this_c = gem.cells[cur];
      VertexType src_dist = gem.min_expect[cur];
      for (IndexType u = this_c.start_vertex, end = this_c.end_vertex; u < end; ++u) {
        VertexType& dst_dist = gem.Vertex(u);
        if (dst_dist > src_dist) {
          local_relaxed = true;
          dst_dist = src_dist;
        }
      }
      for (int i = 1; i < gem.cells.size(); ++i) {
        if (i == cur) continue;
        auto& c = gem.cells[i];
        if (!intersect(c.start_dst_vertex, c.end_dst_vertex, this_c.start_vertex, this_c.end_vertex)) continue;
        VertexType new_dist = src_dist + c.min_edge;
        VertexType& dst_dist = gem.min_expect[i];
        if (dst_dist > new_dist) {
          local_relaxed = true;
          dst_dist = new_dist;
        }
      }
    }
    if (local_relaxed) {
      controller.has_relaxed = true;
    }
  }

  void SkipTask() const {
    VertexType bound = controller.bound;
    while (controller.skip) {
      uint64_t cur = controller.start_cursor.fetch_add(controller.task_block);
      if (!controller.skip || cur >= controller.end_cursor) break;

      for (IndexType 
        start_vertex = controller.cell->start_vertex + cur,
        end_vertex = min((uint64_t)controller.cell->end_vertex, start_vertex + controller.task_block);
        start_vertex < end_vertex; ++start_vertex) {

        if (!gem.in_activity->get_bit(start_vertex)) continue;

        VertexType src_dist = gem.Vertex(start_vertex);
        if (src_dist < bound) {
          controller.skip = false;
          break;
        }
      }
    }
  }

  void MinVertexTask() const {
    while (true) {
      uint64_t cur = controller.start_cursor.fetch_add(controller.task_block);
      if (cur >= controller.end_cursor) break;
      auto& c = gem.cells[cur + 1];
      VertexType min_value = controller.min_vertex[cur];
      for (IndexType start_vertex = c.start_vertex, end_vertex = c.end_vertex; start_vertex < end_vertex; ++start_vertex) {
        if (!gem.in_activity->get_bit(start_vertex)) continue;
        min_value = min(min_value, gem.Vertex(start_vertex));
      }
      controller.min_vertex[cur] = min_value;
    }
  }
};



int main(int argc, char ** argv) {
  cerr.precision(4);

  RGraphgem<VertexType, EdgeType> gem(argv[1]);
  VertexType UPPER_BOUND = gem.TotVertexCnt() * 2.0;
  if (UPPER_BOUND < 0) exit(1);

  cerr << "A --" << gem.rabstract.edge_cnt << "\n";
  if (gem.cells.size() >= 2) {
    cerr << "W --" << gem.cells[1].end_vertex - gem.cells[1].start_vertex << "\n";
    cerr << "B --" << gem.cells[1].edge_cnt << "\n";
  }

  int num_threads = atoi(argv[2]);
  cerr << "num_threads: " << num_threads << "\n";
  
  WorkerController controller(num_threads);
  controller.upper_bound = UPPER_BOUND;
  vector<thread> workers;
  for (int i = 0; i < num_threads; ++i) {
    workers.emplace_back(Worker(i, controller, gem));
  }

  uint64_t task_block = atoll(argv[3]);
  cerr << "task_block: " << task_block << "\n";

  int mrt = atoi(argv[4]);

  {
    ifstream input(argv[5]);
    IndexType src, dst;
    EdgeType dist;
    float reach, final;
    while (input >> src >> dst >> dist >> reach >> final) {
      TEST_CASE.push_back(src); TEST_CASE.push_back(dst);
      gridgraph_dist.push_back(dist);
      gridgraph_reach.push_back(reach);
      gridgraph_final.push_back(final);
    }
  }
  

  vector<double> rgraph_t, normal_t, sorted_t, reach_t;
  printf("src\tdst\tdis\treach\tfinal\n");
  for (int _ = 0; _ < TEST_CASE.size(); _ += 2) {
    IndexType qsrc = TEST_CASE[_], qdst = TEST_CASE[_+1];
    cerr << qsrc << "\t" << qdst << "\t" 
      << gridgraph_dist[_/2] << "\t"
      << gridgraph_reach[_/2] << "\t"
      << gridgraph_final[_/2] << "\t";

    {
      // auto start_t = chrono::steady_clock::now();
      // gem.InitVertex(UPPER_BOUND);
      // fill(gem.min_expect.begin(), gem.min_expect.end(), UPPER_BOUND);
      // gem.Vertex(qdst) = 0;
      // for (int i = 1; i < gem.cells.size(); ++i) {
      //   auto& c = gem.cells[i];
      //   if (qdst < c.start_dst_vertex || qdst >= c.end_dst_vertex) continue;
      //   gem.min_expect[i] = c.min_edge;
      // }

      // controller.has_relaxed = true;
      // while (controller.has_relaxed) {
      //   controller.has_relaxed = false;
      //   controller.InitRAbstractTask(gem, task_block);
      //   controller.StartAllWorkers(); // Start Work
      //   controller.WaitAllWorkers(); // End Work

      //   controller.InitRCellsTask(gem);
      //   controller.StartAllWorkers(); // Start Work
      //   controller.WaitAllWorkers(); // End Work
      // }

      // rgraph_t.push_back(chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - start_t).count()/1000.0);
      // cerr << rgraph_t.back() << "\t";
      // // cerr << "\n";
      // // for (int i = 1; i < gem.cells.size(); ++i) {
      // //   cerr << "\t" << gem.min_expect[i] << "," << gem.cells[i].min_edge;
      // // }
      for (int i = 1; i < gem.cells.size(); ++i) {
        gem.min_expect[i] = gem.cells[i].min_edge;
      }
    }

    // {
    //   auto start_t = chrono::steady_clock::now();
    //   gem.InitVertex(UPPER_BOUND);
    //   gem.Vertex(qsrc) = 0;
      
    //   gem.in_activity->clear();
    //   gem.in_activity->set_bit(qsrc);

    //   AdjCell<EdgeType> *abstract = &gem.cells[0];

    //   controller.has_relaxed = true;
    //   {
    //     vector<AdjCell<EdgeType> *> task{abstract};
    //     controller.InitRelaxEdgeTask(task, task_block);
    //     controller.StartAllWorkers(); // Start Work
    //     controller.WaitAllWorkers(); // End Work
    //   }
      
    //   while (controller.has_relaxed) {
    //     controller.has_relaxed = false;
    //     gem.out_activity->clear();
    //     controller.upper_bound = gem.Vertex(qdst);

    //     for (int i = 1; i < gem.cells.size(); ++i) {
    //       auto& c = gem.cells[i];
    //       // if (c.min_edge >= controller.upper_bound) continue;
    //       if (gem.min_expect[i] >= controller.upper_bound) continue;//HERE
    //       if (!gem.in_activity->get_bits(c.start_vertex, (1 + c.end_vertex) - c.start_vertex)) continue;
    //       // if (controller.upper_bound != UPPER_BOUND && c.min_edge != 1) {
    //       //   controller.InitSkipTask(&c, 1024, controller.upper_bound - c.min_edge);
    //       if (controller.upper_bound != UPPER_BOUND && gem.min_expect[i] != 1) {
    //         controller.InitSkipTask(&c, 1024, controller.upper_bound - gem.min_expect[i]);
    //         controller.StartAllWorkers(); // Start Work
    //         controller.WaitAllWorkers(); // End Work
    //         if (controller.skip) continue;
    //       }

    //       controller.this_relaxed = false;
    //       vector<AdjCell<EdgeType> *> task{&c};
    //       controller.InitRelaxEdgeTask(task, task_block);
    //       controller.StartAllWorkers(); // Start Work
    //       controller.WaitAllWorkers(); // End Work
          
    //       if (controller.this_relaxed) {
    //         vector<AdjCell<EdgeType> *> task{abstract};
    //         if (intersect(c.start_vertex, c.end_vertex, c.start_dst_vertex, c.end_dst_vertex)) {
    //           for (int i = 1; i < mrt; ++i) { task.push_back(&c); /*task.push_back(abstract);*/ };
    //         }
    //         controller.InitRelaxEdgeTask(task, task_block);
    //         controller.StartAllWorkers(); // Start Work
    //         controller.WaitAllWorkers(); // End Work
    //       }
    //     }
    //     swap(gem.in_activity, gem.out_activity);
    //   }

    //   normal_t.push_back(rgraph_t.back() + chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - start_t).count()/1000.0);
    //   cerr << gem.Vertex(qdst) << "\t" << normal_t.back() << "\t";
    // }

    {
      auto start_t = chrono::steady_clock::now();
      gem.InitVertex(UPPER_BOUND);
      gem.Vertex(qsrc) = 0;
  
      gem.in_activity->clear();
      gem.in_activity->set_bit(qsrc);

      AdjCell<EdgeType> *abstract = &gem.cells[0];

      controller.has_relaxed = true;
      {
        vector<AdjCell<EdgeType> *> task{abstract};
        controller.InitRelaxEdgeTask(task, task_block);
        controller.StartAllWorkers(); // Start Work
        controller.WaitAllWorkers(); // End Work
      }
      EdgeType cur_min = gem.Vertex(qdst);
      auto last_t = chrono::steady_clock::now();
      while (controller.has_relaxed) {
        controller.has_relaxed = false;
        gem.out_activity->clear();
        controller.upper_bound = gem.Vertex(qdst);
        
        controller.InitMinVertexTask(gem.cells.size() - 1, controller.upper_bound);
        controller.StartAllWorkers(); // Start Work
        controller.WaitAllWorkers(); // End Work
        vector<pair<VertexType, int> > min_expect;
        for (int i = 1; i < gem.cells.size(); ++i) {
          auto& c = gem.cells[i];
          min_expect.emplace_back(gem.min_expect[i] + controller.min_vertex[i-1], i);
          // if (qdst < c.start_dst_vertex || qdst >= c.end_dst_vertex) {
          //   min_expect.back().first += 1;
          // }
        }
        sort(min_expect.begin(), min_expect.end());

        for (int i = 0; i < min_expect.size(); ++i) {
          if (min_expect[i].first >= controller.upper_bound) break;
          auto& c = gem.cells[min_expect[i].second];

          controller.this_relaxed = false;
          vector<AdjCell<EdgeType> *> task{&c};
          controller.InitRelaxEdgeTask(task, task_block);
          controller.StartAllWorkers(); // Start Work
          controller.WaitAllWorkers(); // End Work
          
          if (controller.this_relaxed) {
            vector<AdjCell<EdgeType> *> task{abstract};
            if (intersect(c.start_vertex, c.end_vertex, c.start_dst_vertex, c.end_dst_vertex)) {
              for (int i = 1; i < mrt; ++i) { task.push_back(&c); /*task.push_back(abstract);*/ };
            }
            controller.InitRelaxEdgeTask(task, task_block);
            controller.StartAllWorkers(); // Start Work
            controller.WaitAllWorkers(); // End Work
          }
        }
        if (cur_min > gem.Vertex(qdst)){
          cur_min = gem.Vertex(qdst);
          last_t = chrono::steady_clock::now();
        }

        if (!controller.has_relaxed) {
          vector<AdjCell<EdgeType> *> task{abstract};
          controller.InitRelaxEdgeTask(task, task_block);
          controller.StartAllWorkers(); // Start Work
          controller.WaitAllWorkers(); // End Work
        }
        swap(gem.in_activity, gem.out_activity);
      }
      if (cur_min > gem.Vertex(qdst)){
        cur_min = gem.Vertex(qdst);
        last_t = chrono::steady_clock::now();
      }
      reach_t.push_back(/*rgraph_t.back() + */chrono::duration_cast<chrono::milliseconds>(last_t - start_t).count()/1000.0);
      sorted_t.push_back(/*rgraph_t.back() + */chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - start_t).count()/1000.0);
      cerr << gem.Vertex(qdst) << "\t"  << reach_t.back() << "\t" << sorted_t.back() << "\t";
    }
    cerr << "\n";
  }

  if (normal_t.size()) {
    sort(normal_t.begin(), normal_t.end());
    cerr << "Normal:\n"
         << "\tmin\t" << normal_t[0] << "\n"
         << "\tavg\t" << accumulate(normal_t.begin(), normal_t.end(), 0.0) / (double)normal_t.size() << "\n"
         << "\tmed\t" << normal_t[normal_t.size()/2] << "\n"
         << "\tmax\t" << normal_t.back() << "\n";
  }
  if (reach_t.size()) {
    sort(gridgraph_reach.begin(), gridgraph_reach.end());
    sort(reach_t.begin(), reach_t.end());
    cerr << "Reach:\n"
         << "\tmin\t" << reach_t[0] << "\n"
         << "\tavg\t" << accumulate(reach_t.begin(), reach_t.end(), 0.0) / (double)reach_t.size() 
         << "\t" <<  accumulate(gridgraph_reach.begin(), gridgraph_reach.end(), 0.0)/accumulate(reach_t.begin(), reach_t.end(), 0.0) << "\n"
         << "\tmed\t" << reach_t[reach_t.size()/2] << "\t"
         << "\t" << gridgraph_reach[gridgraph_reach.size()/2]/reach_t[reach_t.size()/2] << "\n"
         << "\tmax\t" << reach_t.back() << "\n";
  }
  if (sorted_t.size()) {
    sort(sorted_t.begin(), sorted_t.end());
    cerr << "Sorted:\n"
         << "\tmin\t" << sorted_t[0] << "\n"
         << "\tavg\t" << accumulate(sorted_t.begin(), sorted_t.end(), 0.0) / (double)sorted_t.size() 
         << "\t" <<  accumulate(gridgraph_final.begin(), gridgraph_final.end(), 0.0)/accumulate(sorted_t.begin(), sorted_t.end(), 0.0) << "\n"
         << "\tmed\t" << sorted_t[sorted_t.size()/2] << "\t"
         << "\t" << gridgraph_final[gridgraph_final.size()/2]/sorted_t[sorted_t.size()/2] << "\n"
         << "\tmax\t" << sorted_t.back() << "\n";
  }
  print_reads();
  controller.StopAllWorkers();
  for (int i = 0; i < num_threads; ++i) workers[i].join();
}
