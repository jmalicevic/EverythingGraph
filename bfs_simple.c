/*** BFS PUSH, PULL, PUSH/PULL and BFS EDGE CENTRIC  **/
#include "random.h"
#include <cilk/reducer_opadd.h>
#include "parallel_ligra.h"

#define BFS_ROOT 0

struct node* node_list;
uint32_t* dist;



static void print_stats(void);

static struct thread_stats {
	uint64_t tasks, updates;
} thread_stats[64];


static inline void bfs_algo(); 
uint32_t degree_in_frontier = 0;
static short* active;
int switched = 0;
static struct thread_buffer thread_buffers[ALGO_NB_THREADS];
static int iterations = 0;

inline int writeMin(long* curr, long newV) {
	long c; int r =0;
	do c = *curr;
	while (c > newV && !(r = __sync_bool_compare_and_swap(curr,c,newV)) );
	return r; 

}

void bfs_construct(void) {
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

void bfs_destruct(void) {
	uint64_t nodes_discovered=0;

	for(uint64_t i = 0; i < NB_NODES; i++)
		if(dist[i] != 0) nodes_discovered++;

	printf("Total nodes discovered:: %lu\n", nodes_discovered);
}

void bfs_edgecentric(struct node* nodes){

	items_in_frontier = 1;
	uint64_t start_iter, stop_iter;
	while(items_in_frontier != 0) {
		rdtscll(start_iter);
		start_iteration();
		parallel_for(uint32_t i = 0; i < nb_edges; i++) {
			struct edge_t* e = &memblock[i];
			uint32_t src = e->src;
			uint32_t dst = e->dst;

			if(in_frontier[src] && __sync_bool_compare_and_swap(&dist[dst] ,0 , dist[src])) {
				dist[dst] = dist[src] + 1;
				in_frontier_next[dst] = 1;
			}
		}

		rdtscll(stop_iter);
#if CILK
		cilk::reducer_opadd<unsigned long> temp(0);

		parallel_for(uint32_t i = 0; i < NB_NODES; i++) {
			if(in_frontier_next[i] == 1)
				*temp += 1;
			in_frontier[i] = 0;
		}

		items_in_frontier = temp.get_value();
#else
		items_in_frontier = 0;
		parallel_for(uint32_t i = 0; i < NB_NODES; i++) {
			if(in_frontier_next[i] == 1)
				__sync_fetch_and_add(&items_in_frontier,1);
			in_frontier[i] = 0;
		}
#endif
		stop_iteration();

		printf("#iter %d items %u time %f\n", iterations, items_in_frontier, (float)(stop_iter - start_iter)/(float)get_cpu_freq());
		iterations++;
	}



}
static void bfs_pull(uint32_t n_id , struct thread_buffer* b) {

	struct node*n = &nodes[n_id];
	struct node *dst;

	for(uint32_t idx = 0; idx < n->nb_in_edges; idx++) {
		uint32_t dst_id = edge_array_in[n->incoming_edges + idx].dst;
		if(in_frontier[dst_id]) {
			dist[n_id] += 1;
			in_frontier_next[n_id] = 1;
			b->current_buffer_index++;
			b->nb_cleaning+=n->nb_out_edges;
			break;
		}
	}

}
static void bfs_push(uint32_t n_id, struct thread_buffer* b){
	struct node* n = &nodes[n_id];
	struct node* dst;
	for(uint32_t idx = 0; idx < n->nb_out_edges; idx++) {
		uint32_t dst_id = edge_array_out[n->outgoing_edges + idx].dst;
		if(dist[dst_id] == 0) {
			dist[dst_id] += 1;
			if(__sync_bool_compare_and_swap(&in_frontier_next[dst_id], 0, 1)) {
				thread_add_task(b, &nodes[dst_id]);
			}
		}
	}

}

static used void iterator(struct node *nodes) {

	float compute_time, switch_time = 0;
	uint32_t total_items = 1;

	uint64_t iter_start, iter_stop;
	uint64_t prev_nodes=0;
	uint64_t  prev_mode = mode;
	while ( items_in_frontier != 0) {
		rdtscll(iter_start);

		start_iteration();
		uint32_t active_partitions = 0;
		uint32_t total_edges_to_stream = 0;

		uint32_t degree = 0;
		items_in_frontier = 0;
		rdtscll(iter_stop);
		uint32_t items_in_wq = has_work_to_do();
		degree_in_frontier = 0;
		if(mode == PUSH) {

			parallel_for(uint32_t i = 0 ; i < items_in_wq; i++) {
				int id_ = get_threadid(); 
				uint32_t n_id = current_ids[i];
				bfs_push(n_id, &thread_buffers[id_]);	
			}

			parallel_for(uint32_t i = 0; i < ALGO_NB_THREADS; i++) {

				thread_flush(&thread_buffers[i]);
			}	

		}
		else {
			parallel_for(uint32_t i = 0; i < NB_NODES; i++) {
				int id_ = get_threadid();
				if(dist[i] == 0) 
					bfs_pull(i, &thread_buffers[id_]);
			}

		}

		if(mode == PULL)
			for(uint32_t i = 0; i < ALGO_NB_THREADS; i++) {
				items_in_frontier += thread_buffers[i].current_buffer_index;
				degree_in_frontier +=thread_buffers[i]. nb_cleaning;
			}

		stop_iteration_only();
		if(mode == PUSH) {
			items_in_frontier = has_work_to_do();
			degree_in_frontier = wq_old_nb_cleaning();
		}	

		if((switch_mode && mode == PUSH && items_in_frontier + degree_in_frontier > nb_edges/20.) || (mode == PULL && items_in_frontier + degree_in_frontier > nb_edges/20.)) {
			mode = PULL;
			switched = 1;
		}
		else {
			mode = PUSH;
		}
		if(mode == PUSH && switched) {	
			for(int i =0;i < ALGO_NB_THREADS; i++)
				init_thread_buffer(&thread_buffers[i]);

			parallel_for(uint32_t i = 0; i < NB_NODES; i++ ){
				if(in_frontier_next[i]) {
					int id_ = get_threadid();
					thread_add_task(&thread_buffers[id_],&nodes[i]); 
				}
			}
			parallel_for (int i = 0; i < ALGO_NB_THREADS; i++) { 
				thread_flush(&thread_buffers[i]);
			}
			switched = 0;
			stop_iteration_only();

		}

		parallel_for (uint32_t i = 0; i  < NB_NODES; i++)
			in_frontier[i] = 0;
		std::swap(in_frontier, in_frontier_next);
		for(int i =0;i < ALGO_NB_THREADS; i++)
			init_thread_buffer(&thread_buffers[i]);
		total_items += items_in_frontier ;
		rdtscll(iter_stop);
		printf("#Iter %d, items %d , time %f\n", iterations, items_in_frontier, (float)(iter_stop - iter_start)/(float)get_cpu_freq());	
		iterations++;
	}

	printf("Done in %d iterations and found %d nodes\n", iterations, total_items);

}

void bfs_reset(struct node *nodes) {
}

/*
 * Default function that launches a bfs from node 0
 */
void bfs(struct node *nodes) {
	uint64_t construct_start, construct_stop;


	rdtscll(construct_start);

	node_list = nodes;
	dist[BFS_ROOT] = 1;
	in_frontier[BFS_ROOT] = 1;

	items_in_frontier = 1;
	if(!get_task_list())
		init_task_list(NB_NODES);
	reset_task_lists();
	add_task(&nodes[BFS_ROOT]);

	rdtscll(construct_stop);
	printf ("#Task list time %lu, ( %.3f sec) \n", construct_stop - construct_start, ((float)(construct_stop - construct_start) / (float)(get_cpu_freq())) );
	rdtscll(construct_start);
	if(load_mode < 2 || load_mode == 6) {
		printf("Use the BFS grid file, don't know how to run BFS here with this layout\n");
		exit(1);
	}
	if(load_mode > 6) 
		bfs_edgecentric(nodes);
	else 
		iterator(nodes);

	rdtscll(construct_stop);
	printf ("#Algo time %lu, ( %.3f sec) \n", construct_stop - construct_start, ((float)(construct_stop - construct_start) / (float)(get_cpu_freq())) );

	parallel_for(uint32_t i = 0; i < ALGO_NB_THREADS; i++) {
		int id = get_threadid();

		if(id != 0) {
			cpu_set_t mask;
			CPU_ZERO(&mask);
			CPU_SET(id, &mask);
			sched_setaffinity(gettid(), sizeof(mask), &mask);
		}
		init_thread_buffer(&thread_buffers[i]);
	}


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
void *bfs_parallels(void *data){ 

}
struct algo_func current_algo = {
	.reset = bfs_reset, .main = bfs,   .construct = bfs_construct, .destruct = bfs_destruct, 
};
