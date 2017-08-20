
/***** Implements NUMA partitioning for PR and corresponding computation
  Has non NUMA aware code and this should be used to benchmark 
  benefit of NUMA-Awareness
 ***/

#include "random.h" 
#include <numa.h>
#include "parallel_ligra.h"
#include "bitmap.h"
#include <cilk/reducer_opadd.h>
#include "barrier.h"
#define ROOT 0
#define NUMA_PART 2 
#define DAMPING_FACTOR 0.85
#define PAGESIZE 4096
float* curr_rank;
float* prev_rank;
uint32_t* degree;
float one_over_n;
struct  per_thread_task {
	uint32_t* task_array;
	uint32_t curr;
};
uint32_t part_offsets[NUMA_PART + 1];
uint32_t* edge_part_degree[NUMA_PART];
uint64_t* edge_part_offsets[NUMA_PART];
struct edge* edge_array_numa[NUMA_PART];
struct per_thread_task per_thread_tasks[ALGO_NB_THREADS];

#define ALPHA 4
void prnuma_reset(struct node*) {

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
struct x_barrier x_sync;
float adding_constant;  
void prnuma_construct(){
	uint64_t part_start, part_stop;

	init_barrier(&x_sync, ALGO_NB_THREADS);


	free(edge_array_out);
	rdtscll(part_start);
	one_over_n = 1.0/(float)NB_NODES;

	if(!numa){
		prev_rank = (float*) malloc(NB_NODES * sizeof(float));
		curr_rank = (float*) malloc(NB_NODES * sizeof(float));
		degree = (uint32_t*) malloc(NB_NODES * sizeof(uint32_t));
		parallel_for(uint32_t i = 0; i < NB_NODES; i++) {
			prev_rank[i] = 0.15; 
			curr_rank[i] = 0.0; 
			degree[i] = nodes[i].nb_out_edges;
		}
		return;
	}
	part_offsets[0] = 0;
	uint32_t remained = nb_edges ; //+ NB_NODES * ALPHA;

	for(int i = 0; i < NUMA_PART ; i++) {
		uint32_t remained_part = NUMA_PART - i;
		uint32_t expected_chunk = remained / remained_part;
		if(remained_part == 1) 
			part_offsets[i+1] = NB_NODES;
		else {
			uint32_t edges_in = 0;
			for( uint32_t v = part_offsets[i]; v < NB_NODES; v++) {
				edges_in += nodes[v].nb_in_edges;// + ALPHA;
				if (edges_in > expected_chunk) {
					part_offsets[i+1] = v;
					break;
				}
			}
			part_offsets[i+1] = (part_offsets[i+1]) / PAGESIZE * PAGESIZE;
		}
		for (uint32_t v = part_offsets[i]; v < part_offsets[i+1]; v++) {
			remained -= nodes[v].nb_in_edges; /// + ALPHA;
		}
	}
	char* array = (char*) mmap(NULL, sizeof(float) * NB_NODES, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	char* array2 = (char*) mmap(NULL, sizeof(float) * NB_NODES, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	char* array3 = (char*) mmap(NULL, sizeof(uint32_t) * NB_NODES, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
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
		size = sizeof(float) * (end - part_offsets[i]);
		size = size < 4096 ? 4096 : (size % 4096 == 0) ? size : size + size % 4096;
		numa_tonode_memory(array + sizeof(float) * part_offsets[i], size, i);
		numa_tonode_memory(array2 + sizeof(float) * part_offsets[i], size, i);
		size = sizeof(uint32_t) * (end - part_offsets[i]);
		size = size < 4096 ? 4096 : (size % 4096 == 0) ? size : size + size % 4096;

		numa_tonode_memory(array3 + sizeof(uint32_t) * part_offsets[i], size, i);

	}

	curr_rank = (float*) array;
	prev_rank = (float*) array2;
	degree = (uint32_t*) array3;
	curr_work = (uint32_t*) curr_work_t;

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
		degree[i] = nodes[i].nb_out_edges;
		for(uint32_t dst_idx = 0; dst_idx < nodes[i].nb_in_edges; dst_idx++) {
			uint32_t dst = edge_array_in[nodes[i].incoming_edges + dst_idx].dst;
			assert(dst < NB_NODES);
			int node = getNuma_node(i);
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
			if(v > 0) edge_part_offsets[i][v] = edge_part_offsets[i][v-1] + edge_part_degree[i][v-1]; 
		}
		edge_array_numa[i] = (struct edge*) numa_alloc_onnode(edges_per_node[i] * sizeof(struct edge), i);
		assert(NULL != edge_array_numa[i]);
	}

	//this has to be replace with radix sort	
	parallel_for(uint32_t i = 0; i < NB_NODES; i++) {
		//do the actual edge copy 
		int edge_idx[NUMA_PART];
		for(int n = 0; n < NUMA_PART; n++) 
			edge_idx[n] = 0;
		for(uint32_t idx = 0; idx < nodes[i].nb_in_edges; idx++) {
			int node = getNuma_node(i); // can be placed by dst, but for PR this has better load balance
			//edge_array_in[nodes[i].incoming_edges + idx].dst);
			assert(node < NUMA_PART);
			struct edge* e  = &edge_array_numa[node][edge_part_offsets[node][i] + edge_idx[node]];
			e->dst = edge_array_in[nodes[i].incoming_edges + idx].dst;
			assert(e->dst < NB_NODES);
			edge_idx[node]++;

		}


	}	


	free(edge_array_in);	
	parallel_for(uint32_t i = 0; i < NB_NODES; i++) {
		prev_rank[i] = 0.15; 
		curr_rank[i] = 0.0; 
	}
	rdtscll(part_stop);
	printf("Partitioning time %.3f s", (float)(part_stop - part_start)/ (double)get_cpu_freq());
}

uint32_t total =0;
int iterations = 0;
uint32_t work_done = 0;
void get_work(uint32_t *start, uint32_t* stop, uint32_t max) {
	int increment =512;

	if(work_done + increment > max)
		increment = (max - work_done) / ALGO_NB_THREADS;
	if(increment == 0)
		increment = 1;

	*start = __sync_fetch_and_add(&work_done, increment);
	*stop = *start + increment;
	if(*stop > max) *stop = max;

	if(*start >=max ) { *start =  *stop ; return; }

}
int stolen[ALGO_NB_THREADS];
void get_numa_work(uint32_t* start, uint32_t* stop, uint32_t max, int node) {

	int increment = 64;

	if(curr_work[node] + increment > max) {
		increment = (max -curr_work[node] ) / ALGO_NB_THREADS;
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

void* process_interleaved(void* c ){


	uint32_t num = NB_NODES/ ALGO_NB_THREADS;
	int tid = (long)c;
	if(tid != 0) {
		cpu_set_t mask;
		CPU_ZERO(&mask);
		CPU_SET(tid, &mask);
		sched_setaffinity(gettid(), sizeof(mask), &mask);

	}
	uint32_t start, stop;

begin:
	wait_b(&x_sync);
	uint32_t idx;
	get_work(&start, &stop, NB_NODES);
	do {
		for(;start < stop; start++) {
			curr_rank[start] = 0;
			for(uint32_t idx = 0; idx < nodes[start].nb_in_edges; idx++) {
				struct edge *e = &edge_array_in[nodes[start].incoming_edges + idx];
				curr_rank[start] += prev_rank[e->dst] / (float)degree[e->dst]; 
			}
			curr_rank[start] = adding_constant +  DAMPING_FACTOR * curr_rank[start];
		}
		get_work(&start, &stop, NB_NODES);
	}while (start != stop);	
	start = tid * num;
	wait_b(&x_sync);
	if(tid != 0) goto begin;
}
void* process_numa(void* c) {

	int tid = (long)c;
	if(tid !=  0) {
		cpu_set_t mask;
		CPU_ZERO(&mask);
		CPU_SET(tid, &mask);
		sched_setaffinity(gettid(), sizeof(mask), &mask);

	}
	uint32_t thread_id, node;
	uint32_t idx;
	uint64_t my_start, my_stop;
	thread_id = tid; 
	node = numa_node_of_cpu(thread_id);
	uint32_t start, stop;

begin:
	wait_b(&x_sync);
	uint32_t buf_size = part_offsets[node+1] - part_offsets[node];

	get_numa_work(&start, &stop, buf_size, node);

	do{
		start += part_offsets[node];
		stop += part_offsets[node];
		for(;start < stop; start++) {
			curr_rank[start] = 0.0;
			for(uint32_t e_idx = 0; e_idx < edge_part_degree[node][start]; e_idx++) {
				uint32_t dst_idx = edge_array_numa[node][edge_part_offsets[node][start] + e_idx].dst;
				curr_rank[start] += prev_rank[dst_idx] / (float)degree[dst_idx];
			}
			/* No need to steal on own node because we changed partitioning compared to BFS to achieve better performance
			   If edges of vertex go to different NUMA nodes this has to be uncommented 
			 */
			/*****
			  for(int n = 0; n < NUMA_PART; n++) {
			  if(n == node) continue;
			  for(uint32_t e_idx = 0; e_idx < edge_part_degree[n][start]; e_idx++) {
			  uint32_t dst_idx = edge_array_numa[n][edge_part_offsets[n][start] + e_idx].dst;
			  curr_rank[start] += prev_rank[dst_idx] / (float)degree[dst_idx];
			  }
			  } ****/
			curr_rank[start] = adding_constant +  DAMPING_FACTOR * curr_rank[start];
		}
		get_numa_work(&start, &stop, buf_size, node);
	} while(start != stop);
	for(uint32_t n_steal = 0; n_steal < NUMA_PART; n_steal++) { //If done steal from cores on other nodes
		if(n_steal == node) continue;
		buf_size = part_offsets[n_steal + 1] - part_offsets[n_steal];
		get_numa_work(&start, &stop, buf_size, node);
		do{
			start += part_offsets[n_steal];
			stop += part_offsets[n_steal];
			for(;start < stop; start++) {
				curr_rank[start] = 0.0;
				for(uint32_t e_idx = 0; e_idx < edge_part_degree[n_steal][start]; e_idx++) {
					uint32_t dst_idx = edge_array_numa[n_steal][edge_part_offsets[n_steal][start] + e_idx].dst;
					curr_rank[start] += prev_rank[dst_idx] / (float)degree[dst_idx];
				}

				for(int n = 0; n < NUMA_PART; n++) {
					if(n == n_steal) continue;
					for(uint32_t e_idx = 0; e_idx < edge_part_degree[n][start]; e_idx++) {
						uint32_t dst_idx = edge_array_numa[n][edge_part_offsets[n][start] + e_idx].dst;
						curr_rank[start] += prev_rank[dst_idx] / (float)degree[dst_idx];
					}
				}
				curr_rank[start] = adding_constant +  DAMPING_FACTOR * curr_rank[start];
			}
			get_numa_work(&start, &stop, buf_size, n_steal);
		} while(start != stop);	

	}
	wait_b(&x_sync);
	if(tid != 0) goto begin;
}
inline  void prnuma(struct node* nodes){
	if(numa)
		parallel_for(int i = 0; i < NUMA_PART; i++)
			curr_work[i] = 0;

	uint64_t t_start, t_stop;
	adding_constant = (1 - DAMPING_FACTOR) * 1 / (float) NB_NODES;

	__cilkrts_end_cilk(); //pthreads allow better thread placemenet 


	pthread_t threads[ALGO_NB_THREADS];
	for(int i = 1; i < ALGO_NB_THREADS; i++){
		if(numa)
			pthread_create(&threads[i],NULL, process_numa, (void*)i);
		else
			pthread_create(&threads[i],NULL, process_interleaved, (void*)i);
	}

	cpu_set_t mask;
	CPU_ZERO(&mask);
	CPU_SET(0, &mask);
	sched_setaffinity(gettid(), sizeof(mask), &mask);

	rdtscll(t_start);

	while(iterations++ < 10 ) {

		int threads_per_node = ALGO_NB_THREADS / NUMA_PART;
		if(numa) {
			process_numa(0);
			for(int i = 0; i < NUMA_PART; i++)
				curr_work[i] = 0;
		}
		else {
			process_interleaved(0);
			work_done = 0;
		}

		printf("Done with iter %d\n", iterations);

		std::swap(curr_rank,prev_rank);
	}

	rdtscll(t_stop);

	for(int i = 1; i < ALGO_NB_THREADS; i++){
		pthread_cancel(threads[i]);
	}
	printf("Algo time %.3f \n" , (double)(t_stop - t_start)/(double)get_cpu_freq());
}

void prnuma_destruct(){
	if(numa) {
		for(int i = 0; i < NUMA_PART; i++) {
			numa_free(edge_array_numa[i], edges_per_node[i] * sizeof(uint32_t));
			numa_free(edge_part_offsets[i], NB_NODES * sizeof(uint32_t) );
			numa_free(edge_part_degree[i], NB_NODES * sizeof(uint32_t) );
		}
		for(int i = 0; i < ALGO_NB_THREADS; i++){ 
			numa_free(per_thread_tasks[i].task_array, sizeof(uint32_t) * 1024);
		}
	}

	else {
		free(degree);
		free(curr_rank);
		free(prev_rank);

	}
}
struct algo_func current_algo = {
	.reset = prnuma_reset, .main = prnuma,   .construct = prnuma_construct, .destruct = prnuma_destruct
};
