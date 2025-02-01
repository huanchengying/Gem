#define USE_ACTIVITY
#define USE_UPPER_BOUND
// #define USE_INDEX
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
  bool last_work, updated;
  atomic<uint64_t> start_cursor;
  uint64_t end_cursor, start_task, end_task, task_block, task_size;
  CVBarrier start_barrier, end_barrier;
  int task_round;

#ifdef USE_UPPER_BOUND
  VertexType upper_bound;
#endif

  WorkerController(int num_threads) : start_barrier(num_threads + 1), end_barrier(num_threads + 1), last_work(false) {}
  inline void StartAllWorkers() { start_barrier.Wait(); }
  inline void WaitAllWorkers() { end_barrier.Wait(); }
  inline void InitTask(uint64_t start, uint64_t end, uint64_t block, int max_round = 1) {
    start_cursor = start_task = start; 
    end_task = end; task_block = block; 
    task_size = end_task - start_task;
    task_round = max_round;
    end_cursor = start_cursor + task_size * task_round;
  }
  inline void StopAllWorkers() {
    InitTask(0, 0, 1);
    StartAllWorkers();
    last_work = true;
    WaitAllWorkers();
  }
};

uint64_t IOBUFFER;
struct Worker {
  int id;
  WorkerController& controller;
  gem<VertexType, EdgeType>& gem;

  Worker(int id_, WorkerController& controller_, gem<VertexType, EdgeType>& gem_)
    : id(id_), controller(controller_), gem(gem_) {}

  void operator() () const {
    vector<char> buf(IOBUFFER);
    while (!controller.last_work) {
      controller.start_barrier.Wait();
      bool local_updated = false;

#ifdef USE_UPPER_BOUND
      VertexType local_upper_bound = controller.upper_bound;
#endif

      while (true) {
        uint64_t cur_start = controller.start_cursor.fetch_add(controller.task_block);
        if (cur_start >= controller.end_cursor) break;
        uint64_t cur_end = min(cur_start + controller.task_block, controller.end_cursor);

        uint64_t edge_start = (cur_start - controller.start_task) % controller.task_size + controller.start_task;
        uint64_t edge_end = (cur_end - controller.start_task) % controller.task_size + controller.start_task;

#ifdef USE_ACTIVITY
        if (edge_start < edge_end) {
          IndexType start_v = gem.EdgePtr(edge_start)->from;
          IndexType end_v = gem.EdgePtr(edge_end - 1)->from;
          if (!gem.in_activity->get_bits(start_v, (1 + end_v) - start_v)) {
            continue;
          }
        }
#endif
        for (auto 
          e = gem.EdgePtr(edge_start), 
          end = gem.EdgePtr(edge_end), 
          turn_point = gem.EdgePtr(controller.end_task); 
          e != end; ++e) 
        {
          if (e == turn_point) {
            e = gem.EdgePtr(controller.start_task);
            if (e == end) break;
          }
#ifdef USE_ACTIVITY
          if (!gem.in_activity->get_bit(e->from)) {

#ifdef USE_INDEX
            if (edge_start < edge_end) {
              IndexType next_v = e->from + 1;
              size_t next_index = gem.Index(e->from + 1);
              while (next_v < gem.TotVertexCnt() && next_index < edge_end) {
                if (gem.in_activity->get_bit(next_v)) break;
                next_index = gem.Index(++next_v);
              }
              if (next_v >= gem.TotVertexCnt() || next_index >= edge_end) break;
              e = gem.EdgePtr(next_index);
            } else {
              continue;
            }
#else
            continue;
#endif
          }
#endif
          VertexType& src_dist = gem.Vertex(e->from);
          VertexType new_dist = src_dist + e->value;
#ifdef USE_UPPER_BOUND
          if (new_dist >= local_upper_bound) continue;
#endif
          VertexType& dst_dist = gem.Vertex(e->to);
          if (dst_dist > new_dist) {
            local_updated = true;
            dst_dist = new_dist;
#ifdef USE_ACTIVITY
            gem.in_activity->set_bit(e->to);
            gem.out_activity->set_bit(e->to);
#endif
          }
        }
      }
      controller.updated |= local_updated;
      controller.end_barrier.Wait();
    }
  }
};


int main(int argc, char ** argv) {
  gem<VertexType, EdgeType> gem(argv[1]);

  int num_threads = atoi(argv[2]);
  cerr << "num_threads: " << num_threads << "\n";
  
  WorkerController controller(num_threads);
  vector<thread> workers;
  for (int i = 0; i < num_threads; ++i) {
    workers.emplace_back(Worker(i, controller, gem));
  }

  uint64_t available_mem = atoll(argv[3]);
  cerr << "available_mem: " << available_mem << "\n";
  uint32_t partition_cnt = ((2 * sizeof(IndexType) + sizeof(EdgeType)) * gem.TotEdgeCnt()) / available_mem + 1;
  uint64_t partition_step = gem.TotEdgeCnt()/partition_cnt + (gem.TotEdgeCnt() % partition_cnt == 0 ? 0 : 1);
  cerr << "partition_step: " << partition_step << " size: " << (2 * sizeof(IndexType) + sizeof(EdgeType)) * partition_step << "\n";

  uint64_t task_block = atoll(argv[4]);
  cerr << "task_block: " << task_block << "\n";

  int mrt = atoi(argv[5]);

  VertexType UPPER_BOUND = gem.TotVertexCnt() * 1.0;
  if (UPPER_BOUND < 0) exit(1);

  cerr.precision(4);
  vector<double> orig_t, rgraph_t, reach_t, final_t;
  printf("src\tdst\tdis\n");
  // for (int _ = 0; _ < TEST_CASE.size(); _ += 3) {
  //   IndexType qsrc = TEST_CASE[_], qdst = TEST_CASE[_+1];
  //   cerr << qsrc << "\t" << qdst << "\t" << TEST_CASE[_+2] << "\t" << 0.0 << "\t";
  //   orig_t.push_back(0);
  for (int _ = 10; _ > 0; --_) {
    IndexType qsrc = rand()*rand()%gem.TotVertexCnt();
    IndexType qdst = rand()*rand()%gem.TotVertexCnt();
    if (qsrc == qdst) continue;

    // Original Answer
    {
      auto start_t = chrono::steady_clock::now();
      gem.InitVertex(UPPER_BOUND);
      gem.Vertex(qsrc) = 0;
      VertexType upper_bound = UPPER_BOUND;

      for (bool updated = true; updated; ) {
        updated = false;
        for (auto e = gem.EdgePtr(0), end = gem.EdgePtr(gem.TotEdgeCnt()); e < end; ++e) {
          VertexType& src_dist = gem.Vertex(e->from);
          VertexType new_dist = src_dist + e->value;
          VertexType& dst_dist = gem.Vertex(e->to);
          if (dst_dist > new_dist) {
            updated = true;
            dst_dist = new_dist;
          }
        }
      }

      VertexType original_dist = gem.Vertex(qdst);
      if (original_dist == UPPER_BOUND) {++_; continue;}
      orig_t.push_back(chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - start_t).count()/1000.0);
      cerr << qsrc << "\t" << qdst << "\t" << original_dist << "\t" << orig_t.back() << "\t";
    }

    {
      auto start_t = chrono::steady_clock::now();
      gem.InitVertex(UPPER_BOUND);
      gem.Vertex(qsrc) = 0;
      madvise(gem.EdgePtr(0), 1024*1024*32, MADV_WILLNEED);

#ifdef USE_ACTIVITY
      gem.in_activity->clear();
      gem.in_activity->set_bit(qsrc);
#endif
      controller.updated = true;
      while (controller.updated) {
        controller.updated = false;

#ifdef USE_ACTIVITY
        gem.out_activity->clear();
#endif
        for (uint64_t start = 0; start < gem.TotEdgeCnt(); start += partition_step) {
          uint64_t end = min(gem.TotEdgeCnt(), start + partition_step);

#ifdef USE_ACTIVITY
          IndexType start_v = gem.EdgePtr(start)->from;
          IndexType end_v = gem.EdgePtr(end - 1)->from;
          if (!gem.in_activity->get_bits(start_v, (1 + end_v) - start_v)) {
            // cerr << "Skip One Block\n";
            continue;
          }
#endif
          if (end < gem.TotEdgeCnt()) {
            madvise(gem.EdgePtr(end), 1024*1024*32, MADV_WILLNEED);
          }
#ifdef USE_UPPER_BOUND
          controller.upper_bound = gem.Vertex(qdst);
#endif
          controller.InitTask(start, end, task_block, mrt);
          controller.StartAllWorkers(); // Start Work
          controller.WaitAllWorkers(); // End Work
        }

#ifdef USE_ACTIVITY
        swap(gem.in_activity, gem.out_activity);
#endif
      }

      VertexType this_dist = gem.Vertex(qdst);
      final_t.push_back(chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - start_t).count()/1000.0);
      cerr << this_dist << "\t" << final_t.back() << "\t";
    }

    cerr << "\n";
  }

  if (orig_t.size()) {
    sort(orig_t.begin(), orig_t.end());
    cerr << "On Original Graph:\n"
         << "\tmin\t" << orig_t[0] << "\n"
         << "\tavg\t" << accumulate(orig_t.begin(), orig_t.end(), 0.0) / (double)orig_t.size() << "\n"
         << "\tmed\t" << orig_t[orig_t.size()/2] << "\n"
         << "\tmax\t" << orig_t.back() << "\n";

    sort(final_t.begin(), final_t.end());
    cerr << "final:\n"
         << "\tmin\t" << final_t[0] << "\n"
         << "\tavg\t" << accumulate(final_t.begin(), final_t.end(), 0.0) / (double)final_t.size() << "\n"
         << "\tmed\t" << final_t[final_t.size()/2] << "\n"
         << "\tmax\t" << final_t.back() << "\n";
  }
  print_reads();
  controller.StopAllWorkers();
  for (int i = 0; i < num_threads; ++i) workers[i].join();

}
