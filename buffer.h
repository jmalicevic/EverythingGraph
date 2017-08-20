#ifndef BUFFER_H
#define BUFFER_H

#define PRIV_BUFFER_SIZE 2048
/*
 * Public interface
 */
void set_next_0();
void init_task_list(size_t size);
void reset_task_lists(void);

// Task adding (directly to the workqueue)
void add_task(struct node* value);

// Efficient task adding (batch - needs flushing)
struct thread_buffer {
   uint32_t buffer[PRIV_BUFFER_SIZE];
   uint64_t current_buffer_index;
   uint64_t nb_cleaning;
};
void init_thread_buffer(struct thread_buffer * const buffer);
void thread_add_task(struct thread_buffer * const buffer, struct node* value);
void thread_flush(struct thread_buffer * const buffer);

// Function to launch the algorithm
typedef void (*algo_fun_t)(int id, struct thread_buffer *b, struct node *n);

static inline used int in_wq(struct node *dst);
uint64_t wq_has_cleaning(void);
uint64_t wq_old_nb_cleaning(void);
int add_to_next(struct node* value, uint32_t items);
void reset_work(void);

/*
 * Private itf
 */
extern uint32_t *current_ids;
extern size_t current, next;
extern struct rwlock work_lock;

void switch_lists(void);
uint32_t *get_task_list(void);
size_t get_task_list_size(void);

struct work {
   uint64_t start;
   uint64_t stop;
};
struct work get_work(int id);
struct work get_work_nodes(int id);
struct work get_work_init(int id, uint32_t items);
void switch_lists_only(void);
void stop_iteration_only(void);
#define foreach_task_all_nodes(w, ids, dst) \
   for(uint32_t ___privi = w.start; ___privi < w.stop; ___privi++) \
	( dst = ___privi ) 
#define foreach_task(w, ids, dst) \
   for(size_t ___privi = w.start; ___privi < w.stop; ___privi++) \
      	if(( dst =  &nodes[current_ids[___privi]]))

int get_task(size_t i);
int sub_has_more_work(void);
int all_nodes_not_done(void);
int init_has_work(uint32_t items);
void start_iteration(void);
void stop_iteration(void);
void stop_iteration_no_switch(void);
void skip_iteration(void);
uint64_t has_work_to_do(void);
uint64_t next_size(void);
uint64_t get_nb_iterations(void);
uint32_t *get_curr (void);
typedef void (*algo_fun_t)(int id, struct thread_buffer *b, struct node *n);
typedef void (*algo_fun_t2) (int id, struct thread_buffer *b, struct node* n, uint32_t n_id);

static inline used int in_wq(struct node *dst) {
	return 0;   
}


#endif
