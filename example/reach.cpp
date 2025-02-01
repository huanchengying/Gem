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

typedef EmptyProperty EdgeType;
typedef bool VertexType;

vector<IndexType> TEST_CASE;
vector<float> gridgraph_reach;

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

  inline void InitReachTask(vector<AdjCell<EdgeType> *>& c, uint64_t block) {
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
      case 1: ReachTask(); break;
      }
      controller.end_barrier.Wait();
    }
  }

  void ReachTask() const {
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
        if (gem.Vertex(from) != true) continue;
        
        for (auto ptr = c->EdgeBegin(from), end = c->EdgeEnd(from); ptr < end; ++ ptr) {
          VertexType& dst = gem.Vertex(ptr->ver);
          if (!dst) {
            dst = local_relaxed = true;
            gem.in_activity->set_bit(ptr->ver);
            gem.out_activity->set_bit(ptr->ver);
          }
        }
      }
    }
    controller.this_relaxed |= local_relaxed;
    controller.has_relaxed |= local_relaxed;
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

  {
    ifstream input(argv[6]);
    IndexType src, dst;
    float reach;
    while (input >> src >> dst >> reach) {
      TEST_CASE.push_back(src); TEST_CASE.push_back(dst);
      gridgraph_reach.push_back(reach);
    }
  }
  

  vector<double> reach_t, reach2_t;
  printf("src\tdst\tdis\treach\n");
  for (int _ = 0; _ < TEST_CASE.size(); _ += 2) {
    IndexType qsrc = TEST_CASE[_], qdst = TEST_CASE[_+1];
    cerr << qsrc << "\t" << qdst << "\t" 
      << gridgraph_reach[_/2] << "\t";

    {
      auto start_t = chrono::steady_clock::now();
      gem.InitVertex(false);
      gem.Vertex(qsrc) = true;
      gem.in_activity->clear();
      gem.in_activity->set_bit(qsrc);
      
      AdjCell<EdgeType> *abstract = &gem.cells[0];

      controller.has_relaxed = true;
      {
        // vector<AdjCell<EdgeType> *> task;
        // for (int i = 0; i < mrt_abs; ++i) task.push_back(abstract);
        // controller.InitReachTask(task, task_block);
        // controller.StartAllWorkers(); // Start Work
        // controller.WaitAllWorkers(); // End Work

        queue<IndexType> q; q.push(qsrc);
        while (!q.empty()) {
          IndexType x = q.front(); q.pop();
          for (auto ptr = abstract->EdgeBegin(x), end = abstract->EdgeEnd(x); ptr < end; ++ ptr) {
            VertexType& dst_dist = gem.Vertex(ptr->ver);
            if (!dst_dist) {
              dst_dist = true;            
              gem.in_activity->set_bit(ptr->ver);
              q.push(ptr->ver);
            }
          }
        }
      }
      VertexType cur_reach = gem.Vertex(qdst);
      auto last_t = chrono::steady_clock::now();

      while (!cur_reach && controller.has_relaxed) {
        controller.has_relaxed = false;
        gem.out_activity->clear();
        for (int i = 1; i < gem.cells.size(); ++i) {
          auto& c = gem.cells[i];
          if (!gem.in_activity->get_bits(c.start_vertex, c.end_vertex-c.start_vertex)) continue;

          controller.this_relaxed = false;
          vector<AdjCell<EdgeType> *> task{&c};
          controller.InitReachTask(task, task_block);
          controller.StartAllWorkers(); // Start Work
          controller.WaitAllWorkers(); // End Work
            
          cur_reach = gem.Vertex(qdst);
          last_t = chrono::steady_clock::now();
          if (cur_reach) break;

          if (controller.this_relaxed) {
            vector<AdjCell<EdgeType> *> task{abstract};
            if (intersect(c.start_vertex, c.end_vertex, c.start_dst_vertex, c.end_dst_vertex)) {
              for (int i = 1; i < mrt; ++i) { task.push_back(&c); /*task.push_back(abstract);*/ };
            }
            controller.InitReachTask(task, task_block);
            controller.StartAllWorkers(); // Start Work
            controller.WaitAllWorkers(); // End Work

            cur_reach = gem.Vertex(qdst);
            if (cur_reach) break;
            last_t = chrono::steady_clock::now();
          }
        }
        swap(gem.in_activity, gem.out_activity);
      }
      if (!cur_reach){
        cur_reach = gem.Vertex(qdst);
        last_t = chrono::steady_clock::now();
      }
      reach_t.push_back(chrono::duration_cast<chrono::milliseconds>(last_t - start_t).count()/1000.0);
      cerr << gem.Vertex(qdst) << "\t"  << reach_t.back() << "\t";
    }


    {
      auto start_t = chrono::steady_clock::now();
      gem.InitVertex(false);
      gem.Vertex(qsrc) = true;
      gem.in_activity->clear();
      gem.in_activity->set_bit(qsrc);
      
      AdjCell<EdgeType> *abstract = &gem.cells[0];

      controller.has_relaxed = true;
      {
        queue<IndexType> q; q.push(qsrc);
        while (!q.empty()) {
          IndexType x = q.front(); q.pop();
          for (auto ptr = abstract->EdgeBegin(x), end = abstract->EdgeEnd(x); ptr < end; ++ ptr) {
            VertexType& dst = gem.Vertex(ptr->ver);
            if (!dst) {
              dst = true;
              if (ptr->ver == qdst) break;
              q.push(ptr->ver);
              gem.in_activity->set_bit(ptr->ver);
            }
          }
        }
      }
      VertexType cur_reach = gem.Vertex(qdst);
      auto last_t = chrono::steady_clock::now();

      while (!cur_reach && controller.has_relaxed) {
        controller.has_relaxed = false;
        gem.out_activity->clear();
        for (int i = 1; i < gem.cells.size(); ++i) {
          auto& c = gem.cells[i];
          if (!gem.in_activity->get_bits(c.start_vertex, c.end_vertex-c.start_vertex)) continue;

          controller.this_relaxed = false;
          vector<AdjCell<EdgeType> *> task{&c};
          controller.InitReachTask(task, task_block);
          controller.StartAllWorkers(); // Start Work
          controller.WaitAllWorkers(); // End Work
            
          cur_reach = gem.Vertex(qdst);
          last_t = chrono::steady_clock::now();
          if (cur_reach) break;

          if (controller.this_relaxed) {
            vector<AdjCell<EdgeType> *> task{abstract};
            if (intersect(c.start_vertex, c.end_vertex, c.start_dst_vertex, c.end_dst_vertex)) {
              for (int i = 1; i < mrt; ++i) { task.push_back(&c); /*task.push_back(abstract);*/ };
            }
            controller.InitReachTask(task, task_block);
            controller.StartAllWorkers(); // Start Work
            controller.WaitAllWorkers(); // End Work

            cur_reach = gem.Vertex(qdst);
            if (cur_reach) break;
            last_t = chrono::steady_clock::now();
          }
        }
        swap(gem.in_activity, gem.out_activity);
      }
      if (!cur_reach){
        cur_reach = gem.Vertex(qdst);
        last_t = chrono::steady_clock::now();
      }
      reach2_t.push_back(chrono::duration_cast<chrono::milliseconds>(last_t - start_t).count()/1000.0);
      cerr << gem.Vertex(qdst) << "\t"  << reach2_t.back() << "\t";
    }
    cerr << "\n";
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
  if (reach2_t.size()) {
    sort(gridgraph_reach.begin(), gridgraph_reach.end());
    sort(reach2_t.begin(), reach2_t.end());
    cerr << "Reach:\n"
         << "\tmin\t" << reach2_t[0] << "\n"
         << "\tavg\t" << accumulate(reach2_t.begin(), reach2_t.end(), 0.0) / (double)reach2_t.size() 
         << "\t" <<  accumulate(gridgraph_reach.begin(), gridgraph_reach.end(), 0.0)/accumulate(reach2_t.begin(), reach2_t.end(), 0.0) << "\n"
         << "\tmed\t" << reach2_t[reach2_t.size()/2] << "\t"
         << "\t" << gridgraph_reach[gridgraph_reach.size()/2]/reach2_t[reach2_t.size()/2] << "\n"
         << "\tmax\t" << reach2_t.back() << "\n";
  }
  print_reads();
  controller.StopAllWorkers();
  for (int i = 0; i < num_threads; ++i) workers[i].join();
}
