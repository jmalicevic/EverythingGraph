#include "random.h"

uint32_t *current_ids;
size_t current, next;
static uint32_t *ids[2];
static uint64_t size, indexes[2];
static uint64_t sub_start, sub_stop, biterations;
static uint64_t nb_cleaning_in_wq;
static uint64_t old_nb_cleaning;
algo_fun_t algo;
algo_fun_t2 algo2;
algo_fun_t clean;

#if !NAIVE_SPLIT
#define CHUNK_SIZE 1024
static uint64_t work_done, work_done_all_nodes, work_done_init;
#endif

uint32_t* get_curr(void){
	return ids[current];
}
/* Return sum of degrees of nodes in wq */
uint64_t wq_has_cleaning(void) {
	return nb_cleaning_in_wq;
}
uint64_t wq_old_nb_cleaning(void){
	return old_nb_cleaning;
}
void reset_task_lists(void) {
	memset(indexes, 0, sizeof(indexes));
	//biterations = 0;
	current = 0;
	next = 1;
	current_ids = NULL;
}

void init_task_list(size_t _size) {
	size = _size;
	for(int i = 0; i < 2; i++) {
		ids[i] = (uint32_t*) malloc(size * sizeof(*ids[i]));
	}
	reset_task_lists();
}

/* Say that a task is in the next wq.
 * Check that bit with in_wq(task) */
static inline void set_wq_bit(struct node *dst) {
		in_frontier_next[id(dst)] = 1;
}



int add_to_next(struct node* value, uint32_t items){
	uint32_t* curr_ids = ids[next];	
	curr_ids[indexes[next]] = id(value);
	indexes[next]++;
	if(items <= indexes[next]) return 1;
	return 0;
}
/*
 */
void add_task(struct node* value) {
	if(in_wq(value))
		return;

	uint32_t *curr_ids = ids[next];
	curr_ids[indexes[next]] = id(value);
	indexes[next]++;
	set_wq_bit(value);
}


void switch_lists_only(void){
	current = (current+1)%2;
	next = (next+1)%2;
	indexes[next] = 0;
	nb_cleaning_in_wq = 0;

}
void switch_lists(void) {
	current = (current+1)%2;
	next = (next+1)%2;
	indexes[next] = 0;
	nb_cleaning_in_wq = 0;
	char* tmp = in_frontier;
	in_frontier = in_frontier_next;
	in_frontier_next = tmp;
}

uint32_t *get_task_list(void) {
	return ids[current];
}

size_t get_task_list_size(void) {
	return size;
}

void reset_work(void) {
#if !NAIVE_SPLIT
	work_done = 0;
	work_done_all_nodes = 0;
	work_done_init = 0;
#endif
}
struct work get_work_nodes(int id) {
	uint64_t start_idx, end_idx;

#if NAIVE_SPLIT
	int num = ALGO_NB_THREADS;
	uint64_t nb_tasks, nb_tasks_per_thread;

	nb_tasks = NB_NODES;
	nb_tasks_per_thread = nb_tasks / num;
	start_idx = nb_tasks_per_thread * id;
	end_idx = nb_tasks_per_thread * (id + 1);
	if(id == num - 1)
		end_idx = nb_tasks;
#else
	uint64_t increment = CHUNK_SIZE;
	if(CHUNK_SIZE * ALGO_NB_THREADS > NB_NODES - work_done_all_nodes)
		increment = (NB_NODES-work_done_all_nodes) / ALGO_NB_THREADS;
	if(increment == 0)
		increment = 1;

	start_idx = __sync_fetch_and_add(&work_done_all_nodes,increment);
	end_idx = start_idx + increment;
	if(end_idx > NB_NODES)
		end_idx = NB_NODES;
	if(start_idx >= NB_NODES)
		start_idx = NB_NODES;
#endif

	struct work w = { .start = start_idx, .stop = end_idx };
	return w;
}
int get_task(size_t i) {
	return current_ids[i];
}

struct work get_work(int id) {
	uint64_t start_idx, end_idx;

#if NAIVE_SPLIT
	int num = ALGO_NB_THREADS;
	uint64_t nb_tasks, nb_tasks_per_thread;

	nb_tasks = indexes[current];
	nb_tasks_per_thread = nb_tasks / num;
	start_idx = nb_tasks_per_thread * id;
	end_idx = nb_tasks_per_thread * (id + 1);
	if(id == num - 1)
		end_idx = nb_tasks;
#else
	uint64_t increment = CHUNK_SIZE;
	if(CHUNK_SIZE * ALGO_NB_THREADS > indexes[current])
		increment = indexes[current] / ALGO_NB_THREADS;
	if(increment == 0)
		increment = 1;

	start_idx = __sync_fetch_and_add(&work_done,increment);
	end_idx = start_idx + increment;
	if(end_idx > indexes[current])
		end_idx = indexes[current];
	if(start_idx >= indexes[current])
		start_idx = end_idx;
#endif

	struct work w = { .start = start_idx, .stop = end_idx };
	return w;
}
struct work get_work_init(int id, uint32_t items) {
	uint64_t start_idx, end_idx;

#if NAIVE_SPLIT
	int num = ALGO_NB_THREADS;
	uint64_t nb_tasks, nb_tasks_per_thread;

	nb_tasks = items;
	nb_tasks_per_thread = nb_tasks / num;
	start_idx = nb_tasks_per_thread * id;
	end_idx = nb_tasks_per_thread * (id + 1);
	if(id == num - 1)
		end_idx = nb_tasks;
#else
	uint64_t increment = CHUNK_SIZE;
	if(CHUNK_SIZE * ALGO_NB_THREADS > indexes[current])
		increment = items / ALGO_NB_THREADS;
	if(increment == 0)
		increment = 1;

	start_idx = __sync_fetch_and_add(&work_done_init,increment);
	end_idx = start_idx + increment;
	if(end_idx > items)
		end_idx = items;
	if(start_idx >= items)
		start_idx = end_idx;
#endif

	struct work w = { .start = start_idx, .stop = end_idx };
	return w;
}


int sub_has_more_work(void) {
#if NAIVE_SPLIT
	return 0;
#else
	return *(volatile uint64_t*)&work_done < indexes[current];
#endif
}
#if !NAIVE_SPLIT
int all_nodes_not_done(void){
	return *(volatile uint64_t*)&work_done_all_nodes < NB_NODES;
}
#endif

#if !NAIVE_SPLIT
int init_has_work(uint32_t items){
	return *(volatile uint64_t*)&work_done_init < items;
}
#endif

void start_iteration(void) {
	rdtscll(sub_start);
	current_ids = ids[current];
	reset_work();
}

void stop_iteration_no_switch(void) {
	rdtscll(sub_stop);
	biterations++;
}

//to avoid switching the pointers of in_frontier and in-fronteir_next when moving form PULL to PUSH
void stop_iteration_only(void){ 

rdtscll(sub_stop);
old_nb_cleaning = nb_cleaning_in_wq;
switch_lists_only();
}
void stop_iteration(void) {
	rdtscll(sub_stop);

	old_nb_cleaning = nb_cleaning_in_wq;
	switch_lists();
	biterations++;
}

uint64_t get_nb_iterations(void) {
	return biterations;
}

void skip_iteration(void) {
	switch_lists();
}

uint64_t has_work_to_do(void) {
	return indexes[current];
}
uint64_t next_size(void) {
	return indexes[next];
}
void set_next_0(){
	indexes[current] = 0;
}

void init_thread_buffer(struct thread_buffer * const buffer) {
	buffer->current_buffer_index = 0;
	buffer->nb_cleaning = 0;
}

void thread_add_task(struct thread_buffer * const buffer, struct node* value) {
	buffer->nb_cleaning += value->nb_out_edges;  

	if(buffer->current_buffer_index == PRIV_BUFFER_SIZE) {
		thread_flush(buffer);
	}
	buffer->buffer[buffer->current_buffer_index] = id(value);
	buffer->current_buffer_index++;
	buffer->nb_cleaning += value->nb_out_edges;
   	buffer->nb_cleaning++;
}

void thread_flush(struct thread_buffer * const buffer) {
	if(buffer->current_buffer_index == 0)
		return;

	uint64_t curr_idx = __sync_fetch_and_add(&(indexes[next]), buffer->current_buffer_index);
	if(curr_idx + buffer->current_buffer_index > size) {
		memcpy(&(ids[next][curr_idx]), buffer->buffer, buffer->current_buffer_index * sizeof(*buffer->buffer));
	}
	memcpy(&(ids[next][curr_idx]), buffer->buffer, buffer->current_buffer_index * sizeof(*buffer->buffer));
	buffer->current_buffer_index = 0;

	__sync_fetch_and_add(&nb_cleaning_in_wq, buffer->nb_cleaning);
	buffer->nb_cleaning = 0;
}




