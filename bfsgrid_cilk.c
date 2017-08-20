/***Run BFS over a grid ***/

#include "random.h"
#include <cilk/reducer_opadd.h>
#include "parallel_ligra.h"

#define bfs_ROOT_NODE 0

struct node* node_list;
uint32_t* dist;
int switched =0;

static void print_stats(void);

static struct thread_stats {
	uint64_t tasks, updates;
} thread_stats[64];


static inline void bfsgrid_algo(); 
static short* active;
static pthread_t threads[ALGO_NB_THREADS];
static struct thread_buffer thread_buffers[ALGO_NB_THREADS];
static int iterations = 0;
inline std::pair<size_t, size_t> get_partition_range(const size_t vertices, const size_t partitions, const size_t partition_id) {
	const size_t split_partition = vertices % partitions;
	const size_t partition_size = vertices / partitions + 1;
	if (partition_id < split_partition) {
		const size_t begin = partition_id * partition_size;
		const size_t end = (partition_id + 1) * partition_size;
		return std::make_pair(begin, end);
	}
	const size_t split_point = split_partition * partition_size;
	const size_t begin = split_point + (partition_id - split_partition) * (partition_size - 1);
	const size_t end = split_point + (partition_id - split_partition + 1) * (partition_size - 1);
	return std::make_pair(begin, end);
}


inline int writeMin(long* curr, long newV) {
	long c; int r =0;
	do c = *curr;
	while (c > newV && !(r = __sync_bool_compare_and_swap(curr,c,newV)) );
	return r; 

}

static inline void get_bfs(int i, int j)
{
	uint32_t stop_q = offsets[i][j];
	uint32_t tid;
	{
		uint32_t stop =  (j == P - 1 ? (i == P - 1? nb_edges : row_offsets[i+1] ) : row_offsets[i] +  offsets[i][j+1] ); 
		for(uint32_t start = row_offsets[i] + offsets[i][j]; start < stop; start++) {
			struct edge_t* e = &memblock[start];
			uint32_t src = e->src;
			uint32_t dst = e->dst;
			if(in_frontier[src] == 1 && dist[dst] == 0){  
				dist[dst] = dist[src] + 1;
				in_frontier_next[dst] = 1;
			}
		}

	}
}

inline int compute(uint32_t s, uint32_t stop, uint32_t i, uint32_t j) {

	for(uint32_t start=s ; start < stop; start++) {
		struct edge_t* e = &blocks[i][j][start];
		uint32_t src = e->src; 
		uint32_t dst = e->dst;
		if(in_frontier[src] == 1 && dist[dst] == 0) {
			dist[dst] = dist[src] + 1;
			in_frontier_next[dst] = 1;
		}	

	}

}
static inline void bfsgrid_algo() {
	for(uint32_t i = 0; i < P; i++) {
		if(active[i] != 0 )	
			parallel_for(uint32_t j = 0; j < P; j++) {

				if(load_mode == 0 || load_mode == 6)	get_bfs(i,j); // grid is in memblock sorted
				else compute(0,offsets[i][j], i, j); //grid is in blocks accross memory created on load
			}
	}	

}

void bfsgrid_construct(void) {
	if(load_mode != 0 && load_mode != 1 && load_mode !=6) {
		printf("To run this, you need to give grid as the data layout, otherwise work with pagerank_simple\n");
		exit(1);
	}
	uint64_t start,stop;
	rdtscll(start);
	dist = (uint32_t*) malloc(NB_NODES *sizeof(uint32_t));
	parallel_for(uint32_t i = 0; i < NB_NODES; i++) {
		dist[i] = 0;
		in_frontier_next[i] = 0;
		in_frontier[i] = 0;
	}
	rdtscll(stop);
	printf ("#Init time for state array %f\n", (float)(stop - start)/(float)get_cpu_freq());
	memset(thread_stats, 0, 64*sizeof(struct thread_stats));
}

void bfsgrid_destruct(void) {
	uint64_t nodes_discovered=0;

	for(uint64_t i = 0; i < NB_NODES; i++)
		if(dist[i] != 0) nodes_discovered++;

	printf("Total nodes discovered:: %lu\n", nodes_discovered);
}


static used void iterator(struct node *nodes) {
	float compute_time, switch_time = 0;
	uint64_t edges_to_stream =0;
	uint64_t p_active = 0;
	char* tmp;
	uint64_t iter_start, iter_stop;
	uint64_t prev_nodes=0;
	uint64_t  prev_mode = mode;
	while ( items_in_frontier != 0) {
		prev_nodes = items_in_frontier;
		prev_mode = mode;
		rdtscll(iter_start);


		bfsgrid_algo();
		uint32_t active_partitions = 0;
		uint32_t total_edges_to_stream = 0;
		tmp = in_frontier_next;
		in_frontier_next = in_frontier;
		in_frontier = tmp;
		uint32_t degree = 0;
		items_in_frontier = 0;
		rdtscll(iter_stop);
		parallel_for(uint32_t i  = 0; i < P; i++) {
			active[i] = 0 ;
		}
		cilk::reducer_opadd<unsigned long> accum(0);
		cilk::reducer_opadd<unsigned long> temp(0);

		parallel_for(uint32_t i = 0; i < NB_NODES; i++) {
			if(in_frontier[i] == 1) {
				*temp+=1 ;
				int pid = get_partition_id(i);
				active[pid] = 1;
			}
			parallel_for(uint32_t i = 0; i < NB_NODES; i++) {
				in_frontier_next[i] = 0; 
			}
			items_in_frontier = temp.get_value();
			printf("#Iter %d, items %d , time %f\n", iterations, items_in_frontier, (float)(iter_stop - iter_start)/(float)get_cpu_freq());	
			edges_to_stream = 0; p_active = 0; prev_nodes = 0; switch_time = 0;
			iterations++;
		}

	}

}
void bfsgrid_reset(struct node *nodes) {
}

/*
 * Default function that launches a bfs from node 0
 */
void bfsgrid(struct node *nodes) {
	uint64_t construct_start, construct_stop;

	rdtscll(construct_start);	
	active = (short*) malloc( P * sizeof(short));
	memset(active, 0, P * sizeof(short));
	rdtscll(construct_stop);
	printf ("#Time to init active array %lu, ( %.3f sec) \n", construct_stop - construct_start, ((float)(construct_stop - construct_start) / (float)(get_cpu_freq())) );     			


	rdtscll(construct_start);

	node_list = nodes;
	dist[BFS_ROOT] = 1;
	in_frontier[BFS_ROOT] = 1;

	items_in_frontier = 1;
	uint32_t proot = get_partition_id(BFS_ROOT);
	active[proot] = 1;

	rdtscll(construct_stop);
	printf ("#Task list time %lu, ( %.3f sec) \n", construct_stop - construct_start, ((float)(construct_stop - construct_start) / (float)(get_cpu_freq())) );
	rdtscll(construct_start);
	iterator(nodes);


	rdtscll(construct_stop);
	printf ("#Algo time %lu, ( %.3f sec) \n", construct_stop - construct_start, ((float)(construct_stop - construct_start) / (float)(get_cpu_freq())) );

}




static void print_stats(void) {
	/* Stats */
	uint64_t tasks = 0, updates = 0;
	for(size_t i = 0; i < ALGO_NB_THREADS; i++) {
		tasks += thread_stats[i].tasks;
		updates += thread_stats[i].updates;
	}
	printf("\t[bfs - TOTAL] %lu tasks done %lu updates pushed\n", tasks, updates);
}
void *bfsgrid_parallels(void *data){ 

}
struct algo_func current_algo = {
	.reset = bfsgrid_reset, .main = bfsgrid,  .construct = bfsgrid_construct, .destruct = bfsgrid_destruct, 
};
