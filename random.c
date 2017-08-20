#include "random.h"

__thread int __id;  
int isSymmetric = 0;
int createUndir = 0;
int not_processed = 1; 
int load_mode = 2;
int numa = 0;
algo_t algo_phase = ALGO;  
volatile int parallels_done; 
volatile int waiting;       
struct edge* edge_array_out;
struct edge* edge_array_in;

pthread_spinlock_t* locks;
int silent;    // Do not print debug info messages if 1
char* in_frontier_next;
char* in_frontier;
uint32_t items_in_frontier = 0;
uint32_t front_degree = 0;
algo_mode mode= PUSH;

int switch_mode = 1;
int a_mode = 0;

#if defined CREATE_WEIGHT
uint32_t* weights;
uint32_t* weights_in;
#endif


static void usage(void) {
	printf("Usage: ./random -f <graph file> -n <nb_nodes> [-m -u -r -L -w -b -a [mode] -p [bfs_root] -s -N ]\n");
	printf("\t-u: create undirected on load (For example for WCC running with load mode 2 or 3)\n");
	printf("\t-r [switch]: switch between PUSH and PULL where applicable (default 1)\n");
	printf("\t-n: NB_NODES, e.g.,\n");
	printf("\t\trmat20 1048576\n");
	printf("\t\trmat23 8388608\n");
	printf("\t\trmat25 33554432\n");
	printf("\t\trmat27 134217728\n");
	printf("\t-a [mode] : Algo mode, 0 PUSH, 1 PULL, (default 0)\n");
	printf("\t-w: weighted input graph\n");
	printf("\t-r [bfs_root]: BFS & SSSP root\n");
	printf("\t-s:  symmetrict graph, if not given set of incoming edges will be created \n"); 
	printf("\t-N: when running pr_numa or bfs_numa, run the NUMA aware version \n"); 
	printf("\t-m [0 1 2 3 4 5] : default = 2 0 grid sorted; 1 grid nosort; 2 adj sorted; 3 adj nosort; 4 adjacency created without sort; 5 grid fully sorted\n"); 
	_exit(-1);
}

int main(int argc, char **argv) {
	uint64_t start, stop;
	uint64_t load_start, load_stop;
	setlocale(LC_NUMERIC, "");
	int update = 0;

	int c;
	while ((c = getopt (argc, argv, "f:n:m:uwLa:p:srNh")) != -1) {
		switch (c) {
			case 'f':
				filename = optarg;
				break;
			case 'm':
				load_mode = atoi(optarg);
				printf("Load mode = %d\n", load_mode);
				break;
			case 'n':
				NB_NODES = atol(optarg);
				//					#if GRID	
				P = (NB_NODES * 8 / 1024) / 20;
				if(P> 2000) P = 256;
				if(P == 0 ) P = 4;

				printf("P = %d\n", P);
				in_frontier = (char*) malloc(NB_NODES * sizeof(char));
				in_frontier_next = (char* ) malloc(NB_NODES * sizeof(char));
				memset(in_frontier, 0, NB_NODES * sizeof(char));
				memset(in_frontier_next, 0 , NB_NODES * sizeof(char)); 

				break;
			case 'u':
				createUndir = 1;
				printf("#Creating Undir\n");
				break;
			case 'r':
				switch_mode = 0;
				break;
			case 's':
				isSymmetric = 1;
				printf(" no in-edges created\n");
				break;
			case 'a':
				a_mode = atoi(optarg);
				if(a_mode == 0) mode = PUSH;
				else mode = PULL;
				break;
			case 'w':
				weighted_graph = 1;
				break;
			case 'p':
				BFS_ROOT= atol(optarg);
				break;
			case 'N':
				numa = 1;
				break;
			case 'h':
				usage();
			case '?':
				if (optopt == 'f' || optopt == 'n' || optopt == 'a')
					fprintf (stderr, "Option -%c requires an argument.\n", optopt);
				else
					fprintf (stderr, "Unknown option `-%c'.\n", optopt);
				usage();
			default:
				usage();
		}
	}
	if(!NB_NODES || !filename)
		usage();

	printf("#Running with %d update threads and %d algo threads\n", NB_CONCURRENCY, ALGO_NB_THREADS);
	rdtscll(start);

	setWorkers(ALGO_NB_THREADS);
	init(update);
	rdtscll(stop);
	printf("#Total loading time %lu ( %fs )\n", stop - start, ((float)(stop - start))/(float)get_cpu_freq());

	current_algo.construct();
	rdtscll(start);
	{
		rdtscll(load_start);

		current_algo.reset(nodes);
		rdtscll(load_stop);
		printf ("#Total reset  time  %lu ( %fs )\n", load_stop - load_start, ((float)(load_stop - load_start))/(float)get_cpu_freq());

		current_algo.main(nodes);
	}
	rdtscll(stop);
	printf("#Total algo time %lu ( %fs )\n", stop - start, ((float)(stop - start))/(float)get_cpu_freq());

end:
	current_algo.destruct();
	return 0;
}
