#include "random.h"
#include "barrier.h"
#define wcc_ROOT_NODE 0

static void print_stats(void);
static struct thread_stats {
	uint64_t tasks, updates;
} thread_stats[64];

static inline void wcc_algo(int id, struct thread_buffer *b, struct node *n);
uint32_t * parents;
long* components;
int phase = 0;
static pthread_t threads[ALGO_NB_THREADS];
struct x_barrier* xsync;

/*
 * Actual wcc algorithm
 */

inline int writeMin(long* curr, long newV) {
	volatile long c; int r =0;
	do c = *curr;
	while (c > newV && !(r = __sync_bool_compare_and_swap(curr,c,newV)) );
	return r; 

}


static inline void wcc_algo(int _id, struct thread_buffer *b, struct node *n) {
	struct node *dst;
	uint32_t n_id = id(n);;

	for(uint32_t idx = 0; idx< n->nb_out_edges; idx++) { 
		uint32_t dst_id = edge_array_out[n->outgoing_edges + idx].dst;
		if(writeMin(&components[dst_id], components[n_id]) ){ 
			if(__sync_bool_compare_and_swap(&in_frontier_next[dst_id], 0, 1)) 
				thread_add_task(b, &nodes[dst_id]);
		}

	}	
}

void *wcc_parallels(void *data) {
	int id = (long)data;
	__id = id + NB_CONCURRENCY;
	if(id != 0) {
		cpu_set_t mask;
		CPU_ZERO(&mask);
		CPU_SET(id, &mask);
		sched_setaffinity(gettid(), sizeof(mask), &mask);
	}

begin:;
      if(!has_work_to_do())
	      usleep(10000);

      __sync_fetch_and_add(&waiting, 1);
      while(*(volatile int*)&parallels_done == 0) NOP10();
	__sync_fetch_and_add(&waiting, -1);

	uint32_t start,end;
	start = NB_NODES / ALGO_NB_THREADS * id;
	end = start + NB_NODES/ALGO_NB_THREADS;
	if ( id == ALGO_NB_THREADS - 1) end = NB_NODES;
	for(;start < end; start++)
		in_frontier_next[start] = 0;


	wait_b(xsync);
	struct thread_buffer b;
	struct node *n;
	init_thread_buffer(&b);

	do {
		struct work w = get_work(id); 
		foreach_task(w, ids, n) {
			wcc_algo(id, &b, n);
		} 
	} while(sub_has_more_work());


	thread_flush(&b);

	__sync_fetch_and_add(&parallels_done, -1);
	while(*(volatile int*)&parallels_done != 0) NOP10();
	if(id != 0)
		goto begin;
	return NULL;
}

void wcc_construct(void) {
	xsync = (struct x_barrier*) malloc(sizeof(struct x_barrier));
	init_barrier(xsync, ALGO_NB_THREADS);
	components = (long*) malloc( NB_NODES * sizeof(long));

	if(load_mode == 2 || load_mode == 3 || load_mode == 4) 
		for(size_t i = 1; i < ALGO_NB_THREADS; i++)
			pthread_create(&threads[i], NULL, wcc_parallels, (void*)i);
}

void wcc_destruct(void) {
	free(components);
	free(xsync);
	if(load_mode  == 2 || load_mode == 3 || load_mode == 4)
		for(size_t i = 1; i < ALGO_NB_THREADS; i++)
			pthread_cancel(threads[i]);
}

static used void iterator(struct node *nodes, algo_fun_t algo) {
	int iterations = 0;
	{
		do { // iterations
			start_iteration();

			while(*(volatile int*)&waiting != ALGO_NB_THREADS - 1); // wait for all threads
			parallels_done = ALGO_NB_THREADS; // go!
			wcc_parallels((void*)0);


			printf("Iter %d work %u\n", iterations, has_work_to_do());
			stop_iteration();
			iterations++;
		} while(has_work_to_do());
	}
}

void wcc_reset(struct node *nodes) {
	print_stats();
	memset(thread_stats, 0, sizeof(thread_stats));
	for(size_t i = 0; i < NB_NODES; i++) {
		components[i] = i;
		in_frontier[i] = 0;
	}
}


void validate() 
{
	uint32_t max_comp = 0;
	uint32_t* diff_comp = (uint32_t*) malloc(NB_NODES * sizeof(uint32_t));
	parallel_for(uint32_t i = 0; i < NB_NODES; i++ ) {
		diff_comp[i] = 0;

	}               
	if(load_mode == 2 || load_mode ==3) {
		for(int i = 0; i < NB_NODES;i++) {
			for(uint32_t idx = 0; idx < nodes[i].nb_out_edges;idx++) {
				struct edge* e = &edge_array_out[nodes[i].outgoing_edges + idx];
				if(components[i] != components[e->dst] ) {
					printf ("ERROR %d in %d and %d in %d \n", i, components[i], e->dst, components[e->dst]);
					exit(1);
				}	
			}

		}

	}
	else {
		for(size_t i =0; i < nb_edges; i++) {
			struct edge_t* e =  &memblock[i];
			if(components[e->src] != components[e->dst] ) { 
				printf ("ERROR %d in %d and %d in %d \n", e->src, components[e->src], e->dst, components[e->dst]);
				exit(1);
			}
		}
	}
	uint32_t total_comp = 0;
	for(uint32_t i =0 ;i < NB_NODES;i++) {
		diff_comp[components[i]]++;
	}               
	for(uint32_t i =0 ;i < NB_NODES;i++) {
		if(diff_comp[i] > max_comp) max_comp = diff_comp[i];
		if(diff_comp[i] !=0) total_comp++;

	}


	printf("Total comp %d , max comp %d \n", total_comp, max_comp);
	free(diff_comp);

}
/*
 * Default function that launches a wcc from all nodes
 */
void wcc(struct node *nodes) {
	uint32_t changes = NB_NODES;
	uint64_t algo_stop, algo_start;
	int iterations = 0;
	rdtscll(algo_start);

	switch(load_mode){
		case 0: //grid
			for(int i = 0; i < NB_NODES; i++)
				in_frontier_next[i] = 0;
			while(changes) {


				printf("Iter %d , active %d\n", iterations++, changes);
				changes = 0;
				for(uint32_t i = 0; i < P; i++) {
					parallel_for(uint32_t j =0; j < P; j++) {
						uint32_t start = row_offsets[i] + offsets[i][j]; 

						uint32_t stop = j == P -1 ? (i == P -1 ? nb_edges:row_offsets[i+1]) : row_offsets[i] + offsets[i][j+1];
						for(;start < stop; start++ ) {
							struct edge_t* e = &memblock[start];
							uint32_t src = e->src;
							uint32_t dst = e->dst;

							if(writeMin(&components[dst], components[src]) ) {
								in_frontier_next[dst] = 1;
							}

							if(writeMin(&components[src], components[dst])) {
								in_frontier_next[src] = 1;
							}


						}


					}

				}
				for(size_t n= 0; n < NB_NODES; n++) {
					if(in_frontier_next[n] == 1) 
						changes +=1 ;
					in_frontier_next[n] = 0;

				}
			}

			break;
		case 8://edge array
			while(changes > 0 ) {
				printf("Iter %d , active %d\n", iterations, changes);
				changes = 0;
				parallel_for(size_t i =0 ; i < nb_edges; i++) {
					struct edge_t* e = &memblock[i];
					uint32_t src = e->src;
					uint32_t dst = e->dst;

					if(writeMin(&components[dst], components[src]) ) {
						in_frontier_next[dst] = 1;	
					}	

					if(writeMin(&components[src], components[dst])) {
						in_frontier_next[src] = 1;
					}
				}	

				for(size_t i = 0; i < NB_NODES; i++) {
					if(in_frontier_next[i] == 1) 
						changes +=1 ;
					in_frontier_next[i] = 0;

				}
				iterations++;
			}
			break;
		case 2:
		case 3:
		case 4: 
			if(!get_task_list())
				init_task_list(NB_NODES);
			reset_task_lists();


			for(size_t i = 0; i < NB_NODES; i++) {
				add_task(&nodes[i]);
				in_frontier_next[i] = 1;
			}

			stop_iteration();
			iterator(nodes, wcc_algo);
			break;
		default:
			printf("Don't know how to run with this layout\n");
			exit(1);
	}

	rdtscll(algo_stop);
	printf("Algo time %f \n", (float)(algo_stop - algo_start)/(float)get_cpu_freq());
	validate();
}




static void print_stats(void) {
	/* Stats */
}

struct algo_func current_algo = {
	.reset = wcc_reset, .main = wcc,  .construct = wcc_construct, .destruct = wcc_destruct, 
};
