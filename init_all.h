#ifndef INIT_H
#define INIT_H
#define Q 1000 
extern uint32_t* row_offsets;
extern uint32_t** q_map;
extern char* has_q;

struct input {
	uint32_t src;
	uint32_t dst;
#if  WEIGHTED
  	float  weight;
#endif
};
struct edge_t {
	uint32_t src;
	uint32_t dst;
#if  WEIGHTED
	float weight;
#endif
#if ALS_H
	double error;
#endif
};
struct arg_t{
	struct input *_in;
	uint32_t n;
};
struct p_range_t{
	uint32_t start,end;
};
void init(int full);
void* quickSort(void* arg);
size_t item_size(void);
struct input *get_input(size_t pos);
extern uint32_t P;

extern int sorted_graph;
extern int weighted_graph;
extern int skip_loops;
extern edge_t* memblock;
extern char *filename;
extern size_t nb_edges;
extern uint32_t NB_NODES;
extern struct node *nodes;
extern uint32_t BFS_ROOT;

#if NUMA
// **** NUMA
extern struct node* numa_nodes;

//
#endif
/*extern char *in_frontier;
extern char *in_frontier_next;
extern uint32_t items_in_frontier;
extern uint32_t front_degree;
*/ 

extern struct edge_t*** blocks; //[P][P];
extern uint32_t** offsets;
//extern uint32_t ranges[Q][Q];
uint32_t get_partition_id(uint32_t v_id);
uint32_t get_partition_id_Q(uint32_t v_id, uint32_t v_pid);
struct p_range_t get_p_range(uint32_t p_id);
int sort_compare (struct input in1, struct input in2);
uint32_t get_psize(uint32_t i, uint32_t j);
//std::pair<size_t, size_t> get_partition_range(const size_t n, const size_t no_p, const size_t pid);
#endif
