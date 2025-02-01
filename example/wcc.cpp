#include "gem/all.h"

#include <stdio.h>
#include <unordered_map>
#include <queue>
#include <numeric>
#include <chrono>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <map>


using namespace std;
using namespace gem;

typedef EmptyProperty EdgeType;
typedef IndexType VertexType;

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

  inline void InitWCCTask(vector<AdjCell<EdgeType> *>& c, uint64_t block) {
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
};

bool intersect(uint64_t ll, uint64_t lr, uint64_t rl, uint64_t rr) {
  if (lr <= rl) return false;
  if (ll >= rr) return false;
  return true;
}

struct Worker {
  int id;
  WorkerController& controller;
  Reachgem<VertexType, EdgeType>& gem;

  Worker(int id_, WorkerController& controller_, Reachgem<VertexType, EdgeType>& gem_)
    : id(id_), controller(controller_), gem(gem_) {}

  void operator() () const {
    while (!controller.last_work) {
      controller.start_barrier.Wait();
      switch (controller.task_type) {
      case 0: break; // End Task;
      case 1: WCCTask(); break;
      }
      controller.end_barrier.Wait();
    }
  }

  void WCCTask() const {
    bool local_relaxed = false;

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
        VertexType& src = gem.Vertex(from);
        
        for (auto ptr = c->EdgeBegin(from), end = c->EdgeEnd(from); ptr < end; ++ ptr) {
          VertexType& dst = gem.Vertex(ptr->ver);
          if (dst > src) {
            dst = src;
            local_relaxed = true;
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
};



int main(int argc, char ** argv) {
  cerr.precision(4);

  Reachgem<VertexType, EdgeType> gem(argv[1]);
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
  int mrt_abs = atoi(argv[5]);

  gem.in_activity->clear();
  for (int i = 0; i < gem.TotVertexCnt(); ++i) {
    gem.Vertex(i) = i;
    gem.in_activity->set_bit(i);
  }
  auto start_t = chrono::steady_clock::now();
      
  AdjCell<EdgeType> *abstract = &gem.cells[0];
  controller.has_relaxed = true;
  {
    // vector<AdjCell<EdgeType> *> task;
    // for (int i = 0; i < mrt_abs; ++i) task.push_back(abstract);
    // controller.InitWCCTask(task, task_block);
    // controller.StartAllWorkers(); // Start Work
    // controller.WaitAllWorkers(); // End Work
    queue<IndexType> q;
    for (IndexType i = 0; i < gem.TotVertexCnt(); ++i) {
      VertexType& src = gem.Vertex(i);
      for (auto ptr = abstract->EdgeBegin(i), end = abstract->EdgeEnd(i); ptr < end; ++ ptr) {
        VertexType& dst = gem.Vertex(ptr->ver);
        if (dst > src) {
          dst = src;
          gem.in_activity->set_bit(ptr->ver);
        }
      }
      if (i) q.push(i);
    }
    while (!q.empty()) {
      IndexType x = q.front(); q.pop();
      VertexType& src = gem.Vertex(x);
      for (auto ptr = abstract->EdgeBegin(x), end = abstract->EdgeEnd(x); ptr < end; ++ ptr) {
        VertexType& dst = gem.Vertex(ptr->ver);
        if (dst > src) {
          dst = src;
          q.push(ptr->ver);
          gem.in_activity->set_bit(ptr->ver);
        }
      }
    }
  }
  cerr << chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - start_t).count()/1000.0 << "==\n";

  while (controller.has_relaxed) {
    controller.has_relaxed = false;
    gem.out_activity->clear();
    for (int i = 1; i < gem.cells.size(); ++i) {
      auto& c = gem.cells[i];
      if (!gem.in_activity->get_bits(c.start_vertex, c.end_vertex-c.start_vertex)) continue;
      controller.this_relaxed = false;
      vector<AdjCell<EdgeType> *> task{&c};
      controller.InitWCCTask(task, task_block);
      controller.StartAllWorkers(); // Start Work
      controller.WaitAllWorkers(); // End Work

      if (controller.this_relaxed) {
        vector<AdjCell<EdgeType> *> task{abstract};
        controller.InitWCCTask(task, task_block);
        controller.StartAllWorkers(); // Start Work
        controller.WaitAllWorkers(); // End Work
      }
    }
    swap(gem.in_activity, gem.out_activity);
  }
      auto xx = chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - start_t).count()/1000.0;
      map<IndexType, uint64_t> wcc_cnt;
      for (int i = 0; i < gem.TotVertexCnt(); ++i) { 
        wcc_cnt[gem.Vertex(i)] += 1; 
      }
      cerr << wcc_cnt.size() << "\t"  << xx << "\t";
      // for (auto i: wcc_cnt) {
      //   cerr << "\t" << i.first << " " << i.second << "\n";
      // }

    cerr << "\n";
 print_reads();
  controller.StopAllWorkers();
  for (int i = 0; i < num_threads; ++i) workers[i].join();

}
