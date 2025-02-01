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

struct pi{
	uint32_t x,y,z;
};
pi *ne;
int main(int argc, char ** argv) {
  cerr.precision(4);

  Mcstgem<VertexType, EdgeType> gem(argv[1]);
  VertexType UPPER_BOUND = gem.TotVertexCnt() * 2.0;
  if (UPPER_BOUND < 0) exit(1);
  if (gem.cells.size() >= 2) {
    cerr << "W --" << gem.cells[1].end_vertex - gem.cells[1].start_vertex << "\n";
    cerr << "B --" << gem.cells[1].edge_cnt << "\n";
  }
  uint64_t w=atoll(argv[2]);
  ne=(pi *)malloc(sizeof(pi)*w);
  AdjCell<EdgeType> abstract=gem.cells[0];
  uint64_t cnt=abstract.edge_cnt;
  gem.Initroot(gem.TotVertexCnt());
  auto start = chrono::steady_clock::now();
  uint64_t no=1;
  int all=0;
  for(uint64_t i=0;i<cnt;i++){
	while(abstract.index[no]<=i){
		no++;
	}
	int x=no+abstract.start_vertex-1;
	int a=gem.boot(x);
	int y=abstract.edges[i].ver;
	int b=gem.boot(y);
	if(a!=b){
		all++;
		gem.Setroot(b,a);
	}
  }
 while(1){
	int x=0;
	for(int i=1;i<gem.cells.size();i++){
		AdjCell<EdgeType> c=gem.cells[i];
		uint64_t no=1;
		for(uint64_t j=0;j<c.edge_cnt;j++){
			while(c.index[no]<=j){
				no++;
			}
			int x1=no+c.start_vertex-1;
			int a=gem.boot(x1);
			int y1=c.edges[j].ver;
			int b=gem.boot(y1);
			if(a!=b){
				if(x==w){
					if(c.edges[j].val>=ne[w-1].z){
						continue;
					}
					x--;
				}
				int f;
				for(f=x-1;f>=0;f--){
					if(ne[f].z>c.edges[j].val){
						continue;
					}
					break;
				}
				f++;
				for(int q=x;q>f;q--){
					ne[q]=ne[q-1];
				}
				pi pp;
				pp.x=x1;
				pp.y=y1;
				pp.z=c.edges[j].val;
				ne[f]=pp;
				x++;
			}
		}
	}
	if(x==0) break;
	for(int i=0;i<x;i++){
		int a=gem.boot(ne[i].x);
		int b=gem.boot(ne[i].y);
		if(a!=b){
			all++;
			gem.Setroot(b,a);
		}
	}
 }
   auto last_t = chrono::steady_clock::now();
   FILE *fp=fopen("out.txt","a+");
   fprintf(fp,"%.6lf\n",chrono::duration_cast<chrono::milliseconds>(last_t - start).count()/1000.0);
   fclose(fp);
   print_reads();
}
