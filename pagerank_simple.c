#include "random.h"
#include "parallel_ligra.h"
#include <cilk/reducer_opadd.h>
#define DAMPING_FACTOR 0.85
float*prev;
float* rank;
float adding_constant;

uint32_t* degree; //for edge arrays or grids where we dont have the degree as vertex state
static void print_stats(void);

static struct thread_stats {
	uint64_t tasks, updates;
} thread_stats[16];


static inline void pr_algo(); 

static pthread_t threads[ALGO_NB_THREADS];
static struct thread_buffer thread_buffers[ALGO_NB_THREADS];

/*
 * Actual  pr algorithm
 */
int CAS(float * ptr, float oldV, float newV){
	return(__sync_bool_compare_and_swap((long*)ptr, *((long*)&oldV), *((long*)&newV)));
}
void write_add(float* _rank, float  val){
	volatile float newV, oldV;

	do {oldV = *_rank;  newV = oldV+val;}
	while(!CAS(_rank, oldV, newV)); 
}

static uint32_t edges_seen = 0;
//runs PR over all edges stored inside an edgelist so we can still use the loading mode 2 or 3 here
static inline void pr_algo_edgecentric() { 
	//For this to be corect, the degree needs to be counted , it is counted towards the pre-processing time
	parallel_for(uint32_t i  = 0; i < nb_edges;i++) {

		struct edge_t* e = &memblock[i];
		uint32_t src = e->src;
		uint32_t dst = e->dst;
		write_add(&rank[dst], prev[src] / nodes[src].nb_out_edges); 
	}	
}
static inline void pr_algo_push(){
	parallel_for(uint32_t i = 0; i < NB_NODES; i++) {
		struct node* n = &nodes[i];
		for(uint32_t j = 0; j < n->nb_out_edges; j++) {
			uint32_t dst_id = edge_array_out[n->outgoing_edges + j].dst;
			write_add(&rank[dst_id], prev[i]/n->nb_out_edges); // Update SUM of destination
		}
	}
}

static inline void pr_algo_pull() {
	parallel_for(uint32_t i = 0; i < NB_NODES; i++) {
		rank[i] = 0.0;
		struct node* n = &nodes[i];
		for(uint32_t j = 0 ; j < n->nb_in_edges; j++) {
			uint32_t dst_id = edge_array_in[n->incoming_edges + j].dst;
			rank[i] += (float)prev[dst_id] / (float)(degree[dst_id]);  
		}
		rank[i] = adding_constant + DAMPING_FACTOR * rank[i];

	}

}
void pr_construct(void) {
	uint64_t start,stop;
	rdtscll(start);
	adding_constant = (1 - DAMPING_FACTOR) * 1/(float)NB_NODES;
	degree = (uint32_t*) malloc(NB_NODES * sizeof(float));
	rank = (float*) malloc(NB_NODES *sizeof(float));

	prev = (float*) malloc(NB_NODES * sizeof(float));
	float one_over_n = 1.0/(float)NB_NODES;
	parallel_for(uint32_t i = 0; i < NB_NODES; i++){
		prev[i] = 0.15;// one_over_n;
		rank[i] = 0.0;
		degree[i] = nodes[i].nb_out_edges;
	}
	rdtscll(stop);
	printf ("#Init time for state %f\n", (double)(stop - start)/(double)get_cpu_freq());
}

void pr_destruct(void) {

	free(degree);
	free(prev);
	free(rank);

}
static used void iterator(struct node *nodes) {
	int iterations = 0;

	float total_error = 0;
	uint64_t iter_start, iter_stop;
	while(iterations < 10) {
		rdtscll(iter_start);  
		if(mode == PUSH) {
			pr_algo_push();
			parallel_for(uint32_t i = 0; i< NB_NODES; i++) {
				rank[i] = adding_constant + DAMPING_FACTOR * rank[i];
				prev[i] = rank[i];
				rank[i] = 0.0;
			}
		}
		else
			pr_algo_pull();

		iterations++;
		rdtscll(iter_stop)
			printf(" Iter %d time %f\n", iterations, (double)(iter_stop - iter_start) / (double)get_cpu_freq());
	}

}
/*
 * Reset the graph (clean distance, father and workqueue presence
 */
void pr_reset(struct node *nodes) {
}

/*
 * Default function that launches pagerank 
 */
void prgrid(struct node *nodes) {
	uint64_t construct_start, construct_stop;

	rdtscll(construct_start);

	rdtscll(construct_stop);
	printf ("#Time to set the active array %lu, ( %.3f sec) \n", construct_stop - construct_start, ((float)(construct_stop - construct_start) / (float)(get_cpu_freq())) );
	rdtscll(construct_start);

	if(load_mode > 6) 
		pr_algo_edgecentric();
	else	
		iterator(nodes);

	rdtscll(construct_stop);
	printf ("#Algo time %lu, ( %.3f sec) \n", construct_stop - construct_start, ((float)(construct_stop - construct_start) / (float)(get_cpu_freq())) );
	for(int i =0 ;i <10;i++) {
		printf("Rank[%d] = %.5f\n", i, prev[i]);
	}
}



static void print_stats(void) {
	/* Stats */
	uint64_t tasks = 0, updates = 0;
	for(size_t i = 0; i < ALGO_NB_THREADS; i++) {
		tasks += thread_stats[i].tasks;
		updates += thread_stats[i].updates;
	}
	printf("\t[PR - TOTAL] %lu tasks done %lu updates pushed\n", tasks, updates);
}


struct algo_func current_algo = {
	.reset = pr_reset, .main = prgrid,   .construct = pr_construct, .destruct = pr_destruct,
};
