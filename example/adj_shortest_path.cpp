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

typedef int EdgeType;
typedef int VertexType;

vector<int> TEST_CASE{
794578,  3971329, 152,// 11.56 
1545711, 663326,  74,//  9.695 
3107147, 4450106, 140,// 10.75 
4450744, 4207456, 159,// 12.32 
3743533, 3403077, 95,//  11.29 
2004504, 3453695, 62,//  10.82 
2614488, 4370158, 53,//  10.22 
1871967, 2291172, 68,//  11.28 
232250,  2668443, 96,//  9.192 
2982831, 2179521, 37//  11.22
};

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
};

struct Worker {
  int id;
  WorkerController& controller;
  Adjgem<VertexType, EdgeType>& gem;

  Worker(int id_, WorkerController& controller_, Adjgem<VertexType, EdgeType>& gem_)
    : id(id_), controller(controller_), gem(gem_) {}

  void operator() () const {
    while (!controller.last_work) {
      controller.start_barrier.Wait();
      switch (controller.task_type) {
      case 0: break; // End Task;
      case 1: RelaxEdgeTask(); break;
      case 2: SkipTask(); break;
      case 3: MinVertexTask(); break;
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
    controller.this_relaxed |= local_relaxed;
    controller.has_relaxed |= local_relaxed;
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

bool intersect(uint64_t ll, uint64_t lr, uint64_t rl, uint64_t rr) {
  if (lr <= rl) return false;
  if (ll >= rr) return false;
  return true;
}

int main(int argc, char ** argv) {
  cerr.precision(4);

  Adjgem<VertexType, EdgeType> gem(argv[1]);
  VertexType UPPER_BOUND = gem.TotVertexCnt() * 200;
  if (UPPER_BOUND < 0) exit(1);

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

  vector<double> normal_t, sorted_t;
  printf("src\tdst\tdis\n");
  for (int _ = 0; _ < TEST_CASE.size(); _ += 3) {
    IndexType qsrc = TEST_CASE[_], qdst = TEST_CASE[_+1];
    cerr << qsrc << "\t" << qdst << "\t" << TEST_CASE[_+2] << "\t";

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
      
      while (controller.has_relaxed) {
        controller.has_relaxed = false;
        gem.out_activity->clear();
        controller.upper_bound = gem.Vertex(qdst);

        for (int i = 1; i < gem.cells.size(); ++i) {
          auto& c = gem.cells[i];
          if (c.min_edge >= controller.upper_bound) continue;
          if (!gem.in_activity->get_bits(c.start_vertex, (1 + c.end_vertex) - c.start_vertex)) continue;
          if (controller.upper_bound != UPPER_BOUND && c.min_edge != 1) {
            controller.InitSkipTask(&c, 1024, controller.upper_bound - c.min_edge);
            controller.StartAllWorkers(); // Start Work
            controller.WaitAllWorkers(); // End Work
            if (controller.skip) continue;
          }

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
        swap(gem.in_activity, gem.out_activity);
      }

      normal_t.push_back(chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - start_t).count()/1000.0);
      cerr << gem.Vertex(qdst) << "\t" << normal_t.back() << "\t";
    }

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
          min_expect.emplace_back(c.min_edge + controller.min_vertex[i-1], i);
          if (qdst < c.start_dst_vertex || qdst >= c.end_dst_vertex) {
            min_expect.back().first += 1;
          }
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
        if (!controller.has_relaxed) {
          vector<AdjCell<EdgeType> *> task{abstract};
          controller.InitRelaxEdgeTask(task, task_block);
          controller.StartAllWorkers(); // Start Work
          controller.WaitAllWorkers(); // End Work
        }
        swap(gem.in_activity, gem.out_activity);
      }

      sorted_t.push_back(chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - start_t).count()/1000.0);
      cerr << gem.Vertex(qdst) << "\t" << sorted_t.back() << "\t";
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
  if (sorted_t.size()) {
    sort(sorted_t.begin(), sorted_t.end());
    cerr << "Sorted:\n"
         << "\tmin\t" << sorted_t[0] << "\n"
         << "\tavg\t" << accumulate(sorted_t.begin(), sorted_t.end(), 0.0) / (double)sorted_t.size() << "\n"
         << "\tmed\t" << sorted_t[sorted_t.size()/2] << "\n"
         << "\tmax\t" << sorted_t.back() << "\n";
  }
  controller.StopAllWorkers();
  for (int i = 0; i < num_threads; ++i) workers[i].join();
}
