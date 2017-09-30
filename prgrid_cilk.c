#include "random.h"
#include "parallel_ligra.h"
#include <cilk/reducer_opadd.h>
#define DAMPING_FACTOR 0.85
float*prev;
float* rank;
float adding_constant;
int converged = 0;
uint32_t no_in = 0;

uint32_t* degree;
static void print_stats(void);

static struct thread_stats {
	uint64_t tasks, updates;
} thread_stats[16];


static inline void prgrid_algo(); 
static short* active;

static pthread_t threads[ALGO_NB_THREADS];
static struct thread_buffer thread_buffers[ALGO_NB_THREADS];

/*
 * Actual  pr algorithm
 */
template <class ET>
inline bool cas(ET *ptr, ET oldv, ET newv) {
if (sizeof(ET) == 8) {
return __sync_bool_compare_and_swap((long*)ptr, *((long*)&oldv), *((long*)&newv));
} else if (sizeof(ET) == 4) {
return __sync_bool_compare_and_swap((int*)ptr, *((int*)&oldv), *((int*)&newv));
} else {
assert(false);
}
}
template <class ET>
inline void write_add(ET *a, ET b) {
volatile ET newV, oldV;
do {oldV = *a; newV = oldV + b;}
while (!cas(a, oldV, newV));
}

static uint32_t edges_seen = 0;

//offsets here stores the size of the cell to be used with load mode 1 where the grid is not contigious in memory
static inline void prgrid_algo_nosort(){
	for(uint32_t i = 0; i < P; i++) {
		parallel_for(uint32_t j = 0; j < P; j++) 
			if(offsets[i][j] != 0) {
				for(uint32_t start = 0; start < offsets[i][j]; start++) {
					struct edge_t* e = &blocks[i][j][start];
					rank[e->dst] += prev[e->src]/nodes[e->src].nb_out_edges;
				}
			}
	}

}

//offsets here store the actuall offset of a cell, not its size and the grid is actually still in memblock
static inline void prgrid_algo() { 
	for(uint32_t i = 0; i < P; i++) {
		parallel_for(int j = 0; j < P; j++) {
			uint32_t start = row_offsets[i] + offsets[i][j];
			uint32_t stop =  (j == P - 1 ? (i == P - 1? nb_edges : row_offsets[i+1] ) : row_offsets[i] +  offsets[i][j+1] ); 
			for( ; start < stop; start++) 			{
				struct edge_t* e = &memblock[start]; 
				uint32_t src = e->src;
				uint32_t dst = e->dst;
				rank[dst] += (prev[src]/(float)nodes[src].nb_out_edges);
			}

		}
	}
}

//runs PR over all edges stored inside an edgelist so we can still use the loading mode 2 or 3 here
static inline void prgrid_algo_col() { 

	parallel_for(uint32_t i  = 0; i < nb_edges;i++) {
		struct edge_t* e = &memblock[i];
		uint32_t src = e->src;
		uint32_t dst = e->dst;
		write_add(&rank[dst], prev[src] / nodes[src].nb_out_edges);

	}	
}
void prgrid_construct(void) {
	uint64_t start,stop;
	rdtscll(start);
	adding_constant = (1 - DAMPING_FACTOR) * 1/(float)NB_NODES;
	rank = (float*) malloc(NB_NODES *sizeof(float));

	prev = (float*) malloc(NB_NODES * sizeof(float));
	float one_over_n = 1.0/(float)NB_NODES;
	parallel_for(uint32_t i = 0; i < NB_NODES; i++){
		prev[i] = 0.15;//one_over_n;
		rank[i] = 0.0;
	}
	rdtscll(stop);
	printf ("#Init time for state %f\n", (double)(stop - start)/(double)get_cpu_freq());
}

void prgrid_destruct(void) {
	free(prev);
	free(rank);
}


static used void iterator(struct node *nodes) {
	int iterations = 0;

	uint64_t iter_start, iter_stop;
	rdtscll(iter_start);

	while(iterations < 10) {

		//	if(load_mode == 0 || load_mode == 5)
		//		prgrid_algo(); //prgrid_algo_col(); //prgrid_algo();
		//	else if(load_mode == 2 || load_mode == 3 || load_mode >= 6)
		//		 prgrid_algo_col();
		if(load_mode == 1)
			prgrid_algo_nosort();	
		else if(load_mode == 0 || load_mode == 6)
			prgrid_algo();
		else prgrid_algo_col(); // Will run PR over memblock (and edge array) that can be sorted either way depending init_all.c
					// as long as it has not been freed 


		parallel_for(uint32_t i = 0; i< NB_NODES; i++) {
			rank[i] = adding_constant + DAMPING_FACTOR * rank[i];
			prev[i] = rank[i];
			rank[i] = 0.0;
		}
		iterations++;
	}
	rdtscll(iter_stop);
	printf("Iter %d time %f \n", iterations, (double)(iter_stop - iter_start) / (double)get_cpu_freq());

}

void prgrid_reset(struct node *nodes) {
}

/*
 * Default function that launches a bfs from node 0
 */
void prgrid(struct node *nodes) {
	uint64_t construct_start, construct_stop;

	rdtscll(construct_start);
	active = (short*) malloc( P * sizeof(short));
	memset(active, 1, P * sizeof(short));


	rdtscll(construct_stop);
	printf ("#Time to set the active array %lu, ( %.3f sec) \n", construct_stop - construct_start, ((double)(construct_stop - construct_start) / (double)(get_cpu_freq())) );
	rdtscll(construct_start);
	iterator(nodes);
	//	print_stats();

	rdtscll(construct_stop);
	printf ("#Algo time %lu, ( %.3f sec) \n", construct_stop - construct_start, ((double)(construct_stop - construct_start) / (double)(get_cpu_freq())) );

}


/*
 * Rerun
 */
void prgrid_rerun(struct node *nodes) {
	iterator(nodes);
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
	.reset = prgrid_reset, .main = prgrid,   .construct = prgrid_construct, .destruct = prgrid_destruct,
};
