/***** NUMA aware BFS and a version not optimized for NUMA
  Should be used to benchmark the effectivness of NUMA aware data placement
  and computation. Used the partitioning schemes of Gemini and Polymer
  The programming model is also different form what is found in systems and 
  in the paper. The computation now goes through all vertices by all threads
  which update only the neighbours on their NUMA node. 
  This improved performance on high diameter graphs with a balanced
  number of neighbours. 
 *******/
#include "random.h" 
#include <numa.h>
#include "parallel_ligra.h"
#include "bitmap.h"
#include <cilk/reducer_opadd.h>
#include "barrier.h"
#define ROOT 0
#define NUMA_PART 2

uint32_t* parent;
uint32_t* active;
uint32_t* active_next;
uint32_t no_active = 0;
uint32_t activated = 0;
struct  per_thread_task {
	uint32_t* task_array;
	uint32_t curr;
};
uint32_t part_offsets[NUMA_PART];
uint32_t* edge_part_degree[NUMA_PART];
uint64_t* edge_part_offsets[NUMA_PART];
struct edge* edge_array_numa[NUMA_PART];
struct per_thread_task per_thread_tasks[ALGO_NB_THREADS];
void bfsnuma_algo() {


}

void bfsnuma_reset(struct node* nodes) {

}


int getNuma_node(uint32_t v_id) {
	for(int i = 0; i < NUMA_PART; i++) { 
		int end = (i == NUMA_PART -1 ? NB_NODES - 1 : part_offsets[i+1]);	
		if(v_id >= part_offsets[i] && v_id <= end) return i;

	}
}
template <class E>
struct getNumaPartDst { uint32_t operator() (E a) {
	getNuma_node(a.dst);
}
};

template <class E>
struct getNumaPartSrc { uint32_t operator() (E a) {
	getNuma_node(a.src);
}
};
uint32_t edges_per_node[NUMA_PART];
uint32_t* curr_work;
static struct x_barrier x_sync;
uint32_t work_done = 0;
uint32_t total =0;
int iterations = 0;

void init_local(){
	parallel_for(uint32_t i = 0; i < ALGO_NB_THREADS;i++) {
		int node =numa_node_of_cpu(i); 
		if(numa)
			per_thread_tasks[i].task_array = (uint32_t*) numa_alloc_onnode(1024 * sizeof(uint32_t), node);
		else 
			per_thread_tasks[i].task_array = (uint32_t*) malloc(1024 * sizeof(uint32_t));	
		per_thread_tasks[i].curr = 0;
	}

	parallel_for(uint32_t i = 0; i < NB_NODES; i++) {
		parent[i] = NB_NODES;
		active[i] = NB_NODES; 
		active_next[i] = NB_NODES;
	}
	parent[ROOT] = ROOT;
	active[0] = ROOT;
	no_active = 1;
	activated = 0;


}
void bfsnuma_construct(){
	init_barrier(&x_sync, ALGO_NB_THREADS);	

	if(!numa){
		parent = (uint32_t*) malloc(NB_NODES * sizeof(uint32_t));
		active = (uint32_t*) malloc(NB_NODES * sizeof(uint32_t));
		active_next = (uint32_t*) malloc(NB_NODES * sizeof(uint32_t));

		init_local();	
		return;
	}

	for(int i = 0; i < NUMA_PART; i++) {
		part_offsets[i] = i * NB_NODES / NUMA_PART;
	}
	char* array = (char*) mmap(NULL, sizeof(uint32_t) * NB_NODES, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	assert(array != (void*)-1);
	int num_part = numa_num_configured_nodes();
	printf("NUMA partitions: %d\n", num_part);

	char* curr_work_t = (char*) mmap(NULL, sizeof(uint32_t) * NUMA_PART, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);

	for(int i = 0; i < num_part; i++) {
		int size = sizeof(uint32_t);
		int page_size = 4096;
		size = size < page_size ? page_size : (size % page_size == 0) ? size : size + size % page_size;
		printf("Size %d\n", size);
		numa_tonode_memory(curr_work_t + i * sizeof(uint32_t), size, i);


		uint32_t end;
		if(i == num_part - 1) end = NB_NODES;
		else end = part_offsets[i+1];
		size = sizeof(uint32_t) * (end - part_offsets[i]);
		size = size < 4096 ? 4096 : (size % 4096 == 0) ? size : size + size % 4096;
		numa_tonode_memory(array + sizeof(uint32_t) * part_offsets[i], size, i);
	}

	parent = (uint32_t*) array;
	curr_work = (uint32_t*) curr_work_t;
	active = (uint32_t*) malloc(NB_NODES * sizeof(uint32_t) );
	active_next = (uint32_t*) malloc(NB_NODES * sizeof(uint32_t) );


	for(uint32_t i = 0 ; i < NUMA_PART; i++) {
		edge_part_degree[i] = (uint32_t*) numa_alloc_onnode(NB_NODES * sizeof(uint32_t), i);
		assert(NULL!= edge_part_degree[i]);
		parallel_for(uint32_t v = 0; v < NB_NODES; v++) 
			edge_part_degree[i][v] = 0;
		edge_part_offsets[i] = (uint64_t*) numa_alloc_onnode(NB_NODES * sizeof(uint64_t), i);
		assert(NULL != edge_part_offsets[i]);
		parallel_for (uint32_t v = 0; v < NB_NODES; v++)
			edge_part_offsets[i][v] = nb_edges;

	}

	printf("Allocated offsets and degree\n");

	parallel_for(int i = 0; i < NB_NODES; i++) {

		for(uint32_t dst_idx = 0; dst_idx < nodes[i].nb_out_edges; dst_idx++) {
			uint32_t dst = edge_array_out[nodes[i].outgoing_edges + dst_idx].dst;
			int node = getNuma_node(dst);
			assert(node < NUMA_PART);	
			edge_part_degree[node][i]++;
		}

	}
	printf("Counted degrees\n");

	for(int i = 0; i < NUMA_PART; i++) { 
		edge_part_offsets[i][0] = 0;
	}
	parallel_for(uint32_t i = 0; i < NUMA_PART; i++) {
		for(uint32_t v = 0; v < NB_NODES; v++) {
			edges_per_node[i] += edge_part_degree[i][v];
			if(v > 0) edge_part_offsets[i][v] = edge_part_offsets[i][v-1] + edge_part_degree[i][v-1]; //edges_per_node[i];
		}
		edge_array_numa[i] = (struct edge*) numa_alloc_onnode(edges_per_node[i] * sizeof(struct edge), i);
		assert(NULL != edge_array_numa[i]);
	}
	//this has to be replace with radix sort	
	parallel_for(uint32_t i = 0; i < NB_NODES; i++) {
		//do the actual edge copy 
		int edge_idx[NUMA_PART];
		edge_idx[0] = edge_idx[1] = 0;
		for(uint32_t idx = 0; idx < nodes[i].nb_out_edges; idx++) {
			int node = getNuma_node(edge_array_out[nodes[i].outgoing_edges + idx].dst);
			struct edge* e  = &edge_array_numa[node][edge_part_offsets[node][i] + edge_idx[node]];
			e->dst = edge_array_out[nodes[i].outgoing_edges + idx].dst;
			edge_idx[node]++;
		}
	}	
	free(edge_array_out);
	init_local();
}

void flush_array(uint32_t* buffer, uint32_t curr) {
	uint32_t idx = __sync_fetch_and_add(&activated, curr);
	assert(curr != 0);
	assert(idx < NB_NODES);
	memcpy(&(active_next[idx]), buffer, curr* sizeof(uint32_t)); 
}


void get_work(uint32_t* start, uint32_t* stop, uint32_t max) {

	int increment = 1024;

	if(work_done + increment > max) {
		increment = (max - work_done ) / ALGO_NB_THREADS;
	}
	if(increment == 0)
		increment = 1;
	*start = __sync_fetch_and_add(&work_done, increment);
	*stop = *start + increment;
	if(*stop > max) *stop = max;

	if(*start >=max ) { *start =  *stop ; return; }

	return;
}
void get_numa_work(uint32_t* start, uint32_t* stop, uint32_t max, int node) {

	int increment = 64;

	if(curr_work[node] + increment > max) {
		increment = (max -curr_work[node] ) /ALGO_NB_THREADS;
	}
	if(increment == 0)
		increment = 1;
	*start = __sync_fetch_and_add(&curr_work[node], increment);
	*stop = *start + increment;
	if(*stop > max) *stop = max;

	if(*start >=max ) { *start =  *stop ; return; }

	return;
}
int has_work(int node, uint32_t max){
	return *(volatile uint32_t*)&curr_work[node] < max;
}

void* process_interleaved(void* c){
	uint32_t stop,start,thread_id, node;
	int tid = (long)c;
	if(tid != 0) {
		cpu_set_t mask;
		CPU_ZERO(&mask);
		CPU_SET(tid, &mask);
		sched_setaffinity(gettid(), sizeof(mask), &mask);

	}
	thread_id = tid;

b: 
	wait_b(&x_sync);
	fence();
	uint32_t buf_size = no_active;

	get_work(&start, &stop, buf_size);
	do {
		for(;start<stop; start++) {
			uint32_t src = active[start];
			for(uint32_t idx = 0; idx < nodes[src].nb_out_edges; idx++) {
				uint32_t dst_id = edge_array_out[nodes[src].outgoing_edges + idx].dst;
				if(parent[dst_id] == NB_NODES && __sync_bool_compare_and_swap(&parent[dst_id], NB_NODES, src))  {

					per_thread_tasks[thread_id].task_array[per_thread_tasks[thread_id].curr++] = dst_id;
					if(per_thread_tasks[thread_id].curr == 1024) {
						per_thread_tasks[thread_id].curr = 0;
						flush_array(per_thread_tasks[thread_id].task_array, 1024);
					}
				}		
			}
		}

		get_work(&start, &stop, buf_size);
	} while(stop!=start);

	if(per_thread_tasks[thread_id].curr > 0)
		flush_array(per_thread_tasks[thread_id].task_array, per_thread_tasks[thread_id].curr);
	per_thread_tasks[thread_id].curr = 0;

	wait_b(&x_sync);

	if(tid != 0 ) goto b;	
}
void* process(void* c) {
	uint32_t stop,start,thread_id, node;
	int tid = (long)c;

	if(tid != 0) {
		cpu_set_t mask;
		CPU_ZERO(&mask);
		CPU_SET(tid, &mask);
		sched_setaffinity(gettid(), sizeof(mask), &mask);

	}
	thread_id = tid;
	node = numa_node_of_cpu(thread_id);

bfs_begin:
	wait_b(&x_sync);


	uint32_t buf_size = no_active;

	get_numa_work(&start, &stop, buf_size, node);

	do {
		for(; start < stop; start++) {
			uint32_t src = active[start]; 
			assert(src < NB_NODES);
			for(uint32_t e_idx = 0; e_idx < edge_part_degree[node][src]; e_idx++) {
				uint32_t dst_id = edge_array_numa[node][edge_part_offsets[node][src] + e_idx].dst;
				if(parent[dst_id] == NB_NODES && __sync_bool_compare_and_swap(&parent[dst_id], NB_NODES, src))  {
					per_thread_tasks[thread_id].task_array[per_thread_tasks[thread_id].curr++] = dst_id;
					if(per_thread_tasks[thread_id].curr == 1024) {
						per_thread_tasks[thread_id].curr = 0;
						flush_array(per_thread_tasks[thread_id].task_array, 1024);
					}
				}
			}
		}
		get_numa_work(&start, &stop, buf_size, node);
	} while(start != stop);


	/*** This was added after the paper to improve BFS running time at the expense of remote accesses
	It might be suboptimal on really big machines or if the partitioning is changed.

	****/ 
	for(int n_steal = 0; n_steal < NUMA_PART;n_steal++) {
		if (n_steal == node) continue;
		get_numa_work(&start, &stop, buf_size, n_steal);

		do {
			for(; start < stop; start++) {
				uint32_t src = active[start];
				assert(src < NB_NODES);
				for(uint32_t e_idx = 0; e_idx < edge_part_degree[n_steal][src]; e_idx++) {
					uint32_t dst_id = edge_array_numa[n_steal][edge_part_offsets[n_steal][src] + e_idx].dst;
					if(parent[dst_id] == NB_NODES && __sync_bool_compare_and_swap(&parent[dst_id], NB_NODES, src))  {
						per_thread_tasks[thread_id].task_array[per_thread_tasks[thread_id].curr++] = dst_id;
						if(per_thread_tasks[thread_id].curr == 1024) {
							per_thread_tasks[thread_id].curr = 0;
							flush_array(per_thread_tasks[thread_id].task_array, 1024);
						}
					}
				}
			}
			get_numa_work(&start, &stop, buf_size, n_steal);
		} while(start != stop);		

	}

	if(per_thread_tasks[thread_id].curr > 0)
		flush_array(per_thread_tasks[thread_id].task_array, per_thread_tasks[thread_id].curr);
	per_thread_tasks[thread_id].curr = 0;
	wait_b(&x_sync);
	if(tid != 0 ) goto bfs_begin;
}
inline  void bfsnuma(struct node* nodes){

	if(numa) {
		for(int i = 0; i < NUMA_PART; i++)
			curr_work[i] = 0;

	}
	__cilkrts_end_cilk();


	pthread_t threads[ALGO_NB_THREADS];

	for(int i = 1; i < ALGO_NB_THREADS; i++) 
		if(numa)
			pthread_create(&threads[i], NULL, process, (void*)i);
		else 
			pthread_create(&threads[i], NULL, process_interleaved, (void*)i);
	cpu_set_t mask;
	CPU_ZERO(&mask);
	CPU_SET(0, &mask);
	sched_setaffinity(gettid(), sizeof(mask), &mask);


	while(no_active ) {

		uint32_t thread_id, node, start,stop;
		uint32_t idx;


		if(numa) {
			process(0);
			for(int i =0; i < NUMA_PART; i++) 
				curr_work[i] = 0;
		}
		else {
			process_interleaved(0);
			work_done = 0;
		}

		printf("Done with iteration %d. Activated %u\n", iterations, activated);	
		iterations++;
		total += activated;
		no_active = activated;
		activated = 0;

		std::swap(active,active_next); 
	}
	for(int i = 1; i < ALGO_NB_THREADS; i++) 
		pthread_cancel(threads[i]);
}

void bfsnuma_destruct(){

	if(numa){
		for(int i = 0; i < NUMA_PART; i++) {
			numa_free(edge_array_numa[i], edges_per_node[i] * sizeof(uint32_t));
			numa_free(edge_part_offsets[i], NB_NODES * sizeof(uint32_t) );
			numa_free(edge_part_degree[i], NB_NODES * sizeof(uint32_t) );

		}
		for(int i = 0; i < ALGO_NB_THREADS; i++){ 
			numa_free(per_thread_tasks[i].task_array, sizeof(uint32_t) * 1024);
		}

	}
	else{
		free(parent);
		for(int i = 0; i < ALGO_NB_THREADS; i++) 
			free(per_thread_tasks[i].task_array);
	}


}
struct algo_func current_algo = {
	.reset = bfsnuma_reset, .main = bfsnuma,   .construct = bfsnuma_construct, .destruct = bfsnuma_destruct
};
