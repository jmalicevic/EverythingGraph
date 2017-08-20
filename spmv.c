#include "random.h"
#include "parallel_ligra.h"
static int phase = 0;
static int global_phase = 0;
uint32_t* value_in;
uint32_t* value_out;
static inline void spmv_algo(int id, struct thread_buffer *b, struct node *n);
static inline void spmv_algo_push(int id, struct thread_buffer *b, struct node *n);

pthread_t threads[ALGO_NB_THREADS];
int CAS(uint32_t * ptr, uint32_t oldV, uint32_t newV){
	return(__sync_bool_compare_and_swap((uint32_t*)ptr, *((uint32_t*)&oldV), *((uint32_t*)&newV)));
}
void write_add(uint32_t* _rank, uint32_t val){
	volatile uint32_t newV, oldV;

	do {oldV = *_rank;  newV = oldV+val;}
	while(!CAS(_rank, oldV, newV));
}
static inline used void *parallells(void *data) {
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

	struct thread_buffer b;
	struct node *n;
	init_thread_buffer(&b);
	do {
		struct work w = get_work(id);
		foreach_task(w, ids, n) {
			if(mode == PUSH)
				spmv_algo_push(id, &b, n);
			else 
				spmv_algo(id, &b, n);
		}
	} while(sub_has_more_work());

	thread_flush(&b);

	__sync_fetch_and_add(&parallels_done, -1);
	while(*(volatile int*)&parallels_done != 0) NOP10();
	if(id != 0)
		goto begin;
	return NULL;
}

void spmv_construct(void) {
	if(load_mode !=8)
	for(size_t i = 1; i < ALGO_NB_THREADS; i++)
		pthread_create(&threads[i], NULL, parallells, (void*)i);
}

void spmv_destruct(void) {
	if(load_mode != 8)
	for(size_t i = 1; i < ALGO_NB_THREADS; i++)
		pthread_cancel(threads[i]);
}

static used void iterator(struct node *nodes, algo_fun_t algo) {
	int iterations = 0;
	{
		do { // iterations

			start_iteration();

			while(*(volatile int*)&waiting != ALGO_NB_THREADS - 1) NOP10(); // wait for all threads
			parallels_done = ALGO_NB_THREADS; // go!
			parallells((void*)0);


			stop_iteration();
			iterations++;
			phase = 1 - phase;
			global_phase = iterations;
		} while(has_work_to_do());
	}
}

void spmv_reset(struct node *nodes) {
}
uint32_t edges_seen = 0;
/*
 * Actual spmv algorithm
 */
static inline void spmv_algo(int id_, struct thread_buffer *b, struct node *n) {
	struct node *dst;
	struct edge *e;
	if(global_phase>1) return;
	uint32_t n_id = id(n);
	foreach_incoming_edges(n,dst,e){
		value_out[n_id] += 0.001 * value_in[id(dst)];
	}
}
/*
 * Actual spmv algorithm
 */
static inline void spmv_algo_push(int id_, struct thread_buffer *b, struct node *n) {
	struct node *dst;
	struct edge *e;
	if(global_phase>1) return;
	uint32_t n_id = id(n);
	foreach_outgoing_edges(n,dst,e){
		write_add(&value_out[n_id], 0.001 * value_in[id(dst)]);
	}
}
/*
 * Default function that launches a spmv 
 */
void spmv(struct node *nodes) {
	value_out = (uint32_t*) malloc(NB_NODES * sizeof(uint32_t));
	value_in = (uint32_t*) malloc(NB_NODES * sizeof(uint32_t));
	parallel_for(uint32_t i = 0; i < NB_NODES ; i++){
		value_out[i] = 0.0;
		value_in[i] = i;
	}

	if(load_mode == 8) {
		parallel_for(uint32_t i = 0; i < nb_edges; i++) {
			struct edge_t* e = &memblock[i];
			write_add(&value_out[e->src],  value_in[e->dst]);
		}

	}
	else {
		if(!get_task_list())
			init_task_list(NB_NODES);
		reset_task_lists();

		struct node *n;
		struct node_list l = get_all_nodes();
		foreach_node(l, n)
			add_task(n);
		free(l.starting_nodes);

		iterator(nodes, spmv_algo) ;
	}
}

/* Rerun
 */
void spmv_rerun(struct node *nodes) {
	iterator(nodes, spmv_algo );

}



struct algo_func current_algo = {
	.reset = spmv_reset, .main = spmv,  .construct = spmv_construct, .destruct = spmv_destruct
};
