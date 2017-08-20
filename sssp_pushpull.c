#include "random.h"
#include "barrier.h"
#include "stdint.h"
#include "parallel_ligra.h"
#define bfs_ROOT_NODE 0
struct node* node_list;
uint32_t* dist;
static uint32_t MAX = UINT32_MAX - 1;
int switched =0;
static void print_stats(void);
static int iterations=0;
static struct thread_stats {
	uint64_t tasks, updates;
} thread_stats[64];

static inline void sssp_algo(int id_, struct thread_buffer *b, struct node *n);
void validate(void);
uint32_t *parents; //for validation purposes
uint32_t *prev_dist;
static struct thread_buffer thread_buffers[ALGO_NB_THREADS];

int getRand(uint32_t i){
	return (i % 10) + 1;

}
int CAS(long * ptr, long oldV, long newV){
	return(__sync_bool_compare_and_swap((long*)ptr, *((long*)&oldV), *((long*)&newV)));
}



inline int writeMin(uint32_t* curr, uint32_t vNew) {
	uint32_t c; int r =0;
	do c = *curr;
	while (c > vNew && !(r = __sync_bool_compare_and_swap(curr, c, vNew))); 
	return r; 

}
static inline void sssp_algo(int id_, struct thread_buffer *b, struct node *n) {
	struct node *dst;
	uint32_t n_id = id(n);
	uint32_t new_dist = prev_dist[n_id]; 	
	uint64_t offset = nodes[n_id].outgoing_edges;
	for(uint32_t i  = 0; i < nodes[n_id].nb_out_edges; i++) {
		uint32_t dst_id = edge_array_out[offset +i].dst;
		new_dist = prev_dist[n_id] +  1; 
		if(dst_id != n_id) {
			if(writeMin(&dist[dst_id], new_dist)){ 
				if(__sync_bool_compare_and_swap(&in_frontier_next[dst_id], 0, 1)) {  
					thread_add_task(b, nodes+dst_id);
				}
			}
		}
	}

}
void validate(void){
	uint32_t nodes_d = 0;
	for(uint32_t i = 0; i < NB_NODES; i++){

		if(prev_dist[i] != MAX) {
			nodes_d++;
			if(prev_dist[parents[i]] > prev_dist[i] ) printf("ERROR Parent of node %d  is %d and has a greater distance\n", i, parents[i]);

		}
	}
	printf("Nodes checked %d \n", nodes_d);
}
static inline void sssp_algo_pull(int id_, struct thread_buffer *b, uint32_t n_id) {
	struct node *dst;
	struct node *n;
	n = &(nodes[n_id]);;
	int add_me = 0;

	uint64_t offset = nodes[n_id].incoming_edges;
	for(uint32_t i = 0; i < nodes[n_id].nb_in_edges; i++) {
		dst = nodes + edge_array_in[offset + i].dst;
		uint32_t dst_id = id(dst);
		if(in_frontier[dst_id] != 1 || dist[dst_id] == MAX) continue;	
		if(dst_id != n_id) {
			uint32_t new_dist =  prev_dist[dst_id] + edge_array_in[offset + i].weight; 
			parents[n_id] = dst_id;
			if(dist[n_id] > new_dist) { 
				dist[n_id] = new_dist;
				add_me = 1;
			}

		}

	}

	if(add_me == 1){  
		thread_add_task(b, n);
		in_frontier_next[n_id] = 1;

	}

}


void sssp_construct(void) {
	uint64_t start,stop;
	rdtscll(start);
	dist = (uint32_t*) malloc(NB_NODES *sizeof(uint32_t));
	prev_dist = (uint32_t*) malloc(NB_NODES *sizeof(uint32_t));

	parallel_for(uint32_t i = 0; i < NB_NODES; i++) {
		dist[i] = MAX ;

		prev_dist[i] = MAX;
	}			
	rdtscll(stop);
	printf ("#Init time for dist array %f\n", (float)(stop - start)/(float)get_cpu_freq());
}

void sssp_destruct(void) {
	uint64_t nodes_discovered=0;

	for(uint32_t i = 0; i < NB_NODES; i++)
		if(dist[i] != MAX) { 
			nodes_discovered++; 


		} 
	printf("Total nodes discovered:: %lu\n", nodes_discovered);
}
static used void iterator(struct node *nodes, algo_fun_t algo) {
	uint64_t start_iter, stop_iter;
	{

		do { // iterationsi

			start_iteration();

			for(int i = 0; i < ALGO_NB_THREADS;i++) {
				init_thread_buffer(&thread_buffers[i]);
			}

			uint32_t items = has_work_to_do();
			rdtscll(start_iter);
			switch(mode) { 

				case PUSH: 
					parallel_for(uint32_t i = 0; i < items; i++) {
						int _id = __cilkrts_get_worker_number();
						struct node* n = &(nodes[current_ids[i]]);
						sssp_algo(_id, &thread_buffers[_id], n);
					} 

					break;

				case PULL:
					parallel_for(uint32_t i = 0; i < NB_NODES; i++) {

						int _id = __cilkrts_get_worker_number();
						struct node* n = &(nodes[i]);
						sssp_algo_pull(_id, &thread_buffers[_id], i);
					} 


					break;
			}

			parallel_for(int i = 0; i < ALGO_NB_THREADS;i++) {
				thread_flush(&thread_buffers[i]);
			}

			rdtscll(stop_iter);
			if(switch_mode) {
				if(has_work_to_do() + wq_old_nb_cleaning() > nb_edges/20.) 
					mode = PULL;
				else
					mode = PUSH;
			}


			parallel_for(uint32_t i = 0; i < NB_NODES; i++) {
				in_frontier[i] = 0;
				prev_dist[i] = dist[i];
			}		

			stop_iteration(); //_only();
			iterations++;

		} while( has_work_to_do());


	}

}

void sssp_reset(struct node *nodes) {


}

/*
 * Default function that launches sssp from node 0
 */
void sssp(struct node *nodes) {
	uint64_t construct_start, construct_stop;
	parents = (uint32_t*) malloc(NB_NODES * sizeof(uint32_t));
	rdtscll(construct_start);
	parallel_for(uint32_t i =0 ; i < NB_NODES; i++) { in_frontier_next[i] = 0; in_frontier[i] = 0; parents[i] = MAX;}
	reset_work();
	node_list = nodes;
	if(!get_task_list())
		init_task_list(NB_NODES);
	dist[BFS_ROOT] = 0;
	prev_dist[BFS_ROOT] = 0;
	parents[BFS_ROOT] = BFS_ROOT;
	add_task(&nodes[BFS_ROOT]);

	in_frontier_next[BFS_ROOT] = 1;

	stop_iteration();

	rdtscll(construct_stop);
	printf ("#Task list time %lu, ( %.3f sec) \n", construct_stop - construct_start, ((float)(construct_stop - construct_start) / (float)(get_cpu_freq())) );
	rdtscll(construct_start);

	iterator(nodes, sssp_algo);

	rdtscll(construct_stop);
	printf ("#Algo time %lu, ( %.3f sec) \n", construct_stop - construct_start, ((float)(construct_stop - construct_start) / (float)(get_cpu_freq())) );

	uint32_t t = 0;
	for(uint32_t i = 0; i < NB_NODES;i++) {
		if(dist[i] != MAX && dist[i] > t) {
			t = dist[i];
		}
	}
	printf("Max dist = %d and dist[1] is %d\n", t, dist[1]);
}


static void print_stats(void) {
}

struct algo_func current_algo = {
	.reset = sssp_reset, .main = sssp,   .construct = sssp_construct, .destruct = sssp_destruct 
};
