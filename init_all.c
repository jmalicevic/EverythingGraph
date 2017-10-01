#include "random.h"
#include <math.h>

//#if GRID
#include"radixSort_ligra.h"
int sorted_graph = 0, weighted_graph = 0, skip_loops = 1;
char *filename;
typedef pair<uint32_t,uint32_t> intPair;
char *memblock1; // raw pointer to the file content
struct edge_t *memblock;
uint64_t nb_edges;
uint32_t NB_NODES;
struct node *nodes;
#if NUMA
struct node* numa_nodes;
#endif
uint32_t P;
uint32_t* row_offsets;
char *has_q;
uint32_t** q_map;
struct edge_t*** blocks;
uint32_t** offsets;
uint32_t BFS_ROOT;
uint64_t init_start,init_stop;
uint64_t parallel_adds;
uint64_t count_on_load = 1; //TODO parametrize
//uint64_t partition_size, split_point, mod_val;

inline std::pair<size_t, size_t> get_partition_range(const size_t vertices, const size_t partitions, const size_t partition_id) {
	const size_t split_partition = vertices % partitions;
	const size_t partition_size = vertices / partitions + 1;
	if (partition_id < split_partition) {
		const size_t begin = partition_id * partition_size;
		const size_t end = (partition_id + 1) * partition_size;
		return std::make_pair(begin, end);
	}
	const size_t split_point = split_partition * partition_size;
	const size_t begin = split_point + (partition_id - split_partition) * (partition_size - 1);
	const size_t end = split_point + (partition_id - split_partition + 1) * (partition_size - 1);
	return std::make_pair(begin, end);
}

size_t item_size(void) {
	if(weighted_graph) {
		return sizeof(struct input);
	} else {
		return 2*sizeof(uint32_t);
	}
}
uint32_t get_psize(uint32_t i,uint32_t j ) {
	return (j == P - 1 ? (i == P - 1 ? nb_edges : row_offsets[i + 1]) - (row_offsets[i] + offsets[i][j]) : offsets[i][j+ 1] - offsets[i][j]);
}

uint32_t get_partition_id(uint32_t v_id){
	if(mod_val == 0 ) { //NB_NODES % P == 0)  
		//	size_t partition_size = NB_NODES / P;
		return v_id / partition_size;
	}
	return (v_id < split_point) ? v_id / partition_size : (v_id - split_point) / (partition_size -1) + mod_val; // (NB_NODES % P);
}
uint32_t get_partition_id_Q(uint32_t v_id, uint32_t v_pid){

	std::pair<size_t,size_t> range = get_partition_range(NB_NODES, P, v_pid);
	uintE n = range.second - range.first;
	if(n % Q == 0 ) { 
		uintE p_size = n / Q;
		return (v_id-range.first) / p_size ; 
	}   
	v_id -= range.first;
	size_t partition_size = n / Q + 1;  
	size_t split_point =  n % Q * partition_size;
	return (v_id < split_point) ? ((v_id) / partition_size ) : ((v_id - split_point) / (partition_size -1) + (n % Q)) ;   
}


struct input *get_input(size_t pos) {
	return (struct input *) &(memblock1[pos * item_size()]);
}
uint64_t load_start,load_stop;

void preload_graph() {


	parallel_for(uint64_t i  = 0; i < nb_edges; i++) {
		struct input* in = get_input(i);
		struct edge_t* e = &memblock[i];
		e->src = in->src;
		e->dst = in->dst;


		if(count_on_load) { // Degree counting  overlapped with load, should be used when loading graph from slow storage
			__sync_fetch_and_add(&nodes[e->src].nb_out_edges,1);
			if (!isSymmetric) __sync_fetch_and_add(&nodes[e->dst].nb_in_edges,1);
		}
#if defined ALS_H
		e->error = 0.0;
#endif

#if WEIGHTED
		if(weighted_graph) e->weight = in->weight;
#endif
		if(createUndir) {
			struct edge_t* e_rev = &memblock[i + nb_edges];
			e_rev->src = in->dst;
			e_rev->dst = in->src;
		}

	}

}
void init_grid_sort_src(int full, int sort_by_src) {
	offsets = (uint32_t**) malloc (P * sizeof(uint32_t*));

	for(size_t i =0; i< P; i++) {
		offsets[i] = (uint32_t*) malloc(P* sizeof(uint32_t));
	}


	for(size_t i = 0; i < P; i++) {
		for(size_t j = 0; j < P; j++) {
			offsets[i][j] = 0;
		}  
	}


	rdtscll(init_start);
	if(not_processed) {


		parallel_for(uint32_t i = 0; i < P; i++) { row_offsets[i] = nb_edges;}

		//Fully sort by source first
		rdtscll(load_start);

		intSort::iSort(memblock, nb_edges, NB_NODES + 1, getEdgeSrc<struct edge_t>()); //PartitionIdSrc<struct edge_t>());//EdgeSrc<struct edge_t>());
		rdtscll(load_stop);
		printf("#Sort by src time by src %f \n" , (float)(load_stop-load_start) / (float)get_cpu_freq());
		row_offsets[0] = 0; 

		rdtscll(load_start);
		parallel_for(uint32_t i = 1 ; i < nb_edges; i++) {
			struct edge_t* e = &memblock[i];
			uint32_t psrc = get_partition_id(e->dst);
			uint32_t pprev = get_partition_id(memblock[i-1].dst);
			if(psrc != pprev) {
				row_offsets[psrc] = i;
			} 
		}


		rdtscll(load_stop);

		printf("# Offset big  time by src %f \n" , (float)(load_stop-load_start) / (float)get_cpu_freq());			
		rdtscll(load_start);
		sequence::scanIBack(row_offsets,row_offsets,(int)P,minF<uintT>(),(uintT)nb_edges);   
		rdtscll(load_stop);
		printf("#Sequence big time by src %f \n" , (float)(load_stop-load_start) / (float)get_cpu_freq());

		uint32_t* node_offsets = (uint32_t* ) malloc(NB_NODES * sizeof(uint32_t));

		rdtscll(load_start);
		parallel_for (uint32_t i = 0; i < NB_NODES; i ++) node_offsets[i] = nb_edges;
		node_offsets[memblock[0].src] = 0;
		if(memblock[0].src != 0) {
			uint32_t i = 0;
			while (i < memblock[0].src) 
				node_offsets[i++] = 0;
		}
		parallel_for(uint32_t i = 1; i < nb_edges; i++ ) {
			uint32_t src = memblock[i].dst;
			uint32_t prev = memblock[i-1].dst;
			if(src != prev) {
				node_offsets[src] = i;

			}
		}
		sequence::scanIBack(node_offsets,node_offsets,(int)NB_NODES,minF<uintT>(),(uintT)nb_edges);   


		parallel_for(uint32_t i = 0; i < NB_NODES; i++) {
			nodes[i].nb_out_edges = i == NB_NODES -1 ? nb_edges : node_offsets[i+1] - node_offsets[i];
		}
		rdtscll(load_stop);
		printf("#Out degree count %f \n" , (float)(load_stop-load_start) / (float)get_cpu_freq());  
		free(node_offsets);
		rdtscll(load_start);


		parallel_for(uint32_t i = 0; i < P; i++) {
			uint32_t o = row_offsets[i];
			uint32_t l = ((i == P - 1) ? nb_edges - row_offsets[i] : row_offsets[i+1] - row_offsets[i]) ;  
			if (l != 0) 
				intSort::iSort(memblock + row_offsets[i],l, NB_NODES + 1, getEdgeSrc<struct edge_t>()); // getPartitionIdDst<struct edge_t>());
		}
		rdtscll(load_stop);
		printf("#Sort by dst time - by src %f \n" , (float)(load_stop-load_start) / (float)get_cpu_freq());
		parallel_for(uint32_t i = 0; i < P; i++) 
			for(uint32_t j = 0; j < P; j++) 
				offsets[i][j] = nb_edges+1  ;

		rdtscll(load_start);	
		for(uint32_t i = 0; i < P; i++) {
			uint32_t start, idx;
			offsets[0][i] =0 ;
			uint32_t stop = ( i == P - 1?  nb_edges  : row_offsets[i + 1] ) - row_offsets[i];
			parallel_for(uint32_t start = 1;  start < stop; start++) {
				struct edge_t* e = &memblock[row_offsets[i] + start];
				uint32_t pdst = get_partition_id(e->src);
				uint32_t pprev = get_partition_id(memblock[row_offsets[i] + start - 1].src);
				if(pdst != pprev) {
					offsets[pdst][i] = start;
				}

			} 
		}


		rdtscll(load_stop);
		printf("#Offset dst time - by src %f \n" , (float)(load_stop-load_start) / (float)get_cpu_freq());
		rdtscll(load_start);
		for(int j = P - 1; j >=0 ; j--) {
			for(int i = P - 1; i>= 0 ; i--) {
				if(offsets[i][j] == nb_edges  + 1){
					offsets[i][j] = ( i == P-1 ? (j == P-1 ? nb_edges : row_offsets[j+1]) : offsets[i+1][j]); 
				}
				//	sequence::scanIBack(offsets[i], offsets[i], (int)P, minF<uintT>(), (uintT)(i == P - 1 ? nb_edges - row_offsets[i] : row_offsets[i+1] - row_offsets[i]));
			}
		}
		rdtscll(load_stop);
		printf("#Seq dst time - by src %f \n" , (float)(load_stop-load_start) / (float)get_cpu_freq());
		rdtscll(init_stop);
		printf("#Total grid create time %lu ( %f s)\n", init_stop - init_start, ((float)(init_stop - init_start))/(float)get_cpu_freq());
		/*	if(sort_by_src){
			rdtscll(load_start);
			parallel_for(int i = 0; i < P; i++) {
			for(int j = 0; j < P; j ++) {
			uint32_t start = row_offsets[i] + offsets[i][j];
			uint32_t stop = (j == P-1 ? (i == P -1 ? nb_edges : row_offsets[i+1]) : row_offsets[i] + offsets[i][j+1] ) ;
			intSort::iSort(memblock + start,stop - start, NB_NODES + 1, getEdgeSrc<struct edge_t>());
			}
			}
			rdtscll(load_stop);
			printf("#Sort cells by src time %lu ( %f s)\n", load_stop - load_start, ((float)(load_stop - load_start))/(float)get_cpu_freq());


			} */
	}

}


void init_grid_sort(int full, int sort_by_src) {
	offsets = (uint32_t**) malloc (P * sizeof(uint32_t*));

	for(size_t i =0; i< P; i++) {
		offsets[i] = (uint32_t*) malloc(P* sizeof(uint32_t));
	}


	for(size_t i = 0; i < P; i++) {
		for(size_t j = 0; j < P; j++) {
			offsets[i][j] = 0;
		}  
	}
	edge_array_out = (struct edge*) calloc(nb_edges ,sizeof (struct edge));


	rdtscll(init_start);
	if(true) {

		parallel_for(uint32_t i = 0; i < P; i++) { row_offsets[i] = nb_edges;}

		//Fully sort by source first
		rdtscll(load_start);
		intSort::iSort(memblock, nb_edges, NB_NODES + 1, getEdgeSrc<struct edge_t>());//getPartitionSrc<struct edge_t>()); //change back to memblock , als modification
		rdtscll(load_stop);
		printf("#Sort by src time - by src %f \n" , (float)(load_stop-load_start) / (float)get_cpu_freq());
		row_offsets[0] = 0; 

		rdtscll(load_start);
		parallel_for(uint32_t i = 1 ; i < nb_edges; i++) {
			struct edge_t* e = &memblock[i];
			uint32_t psrc = get_partition_id(e->src);
			uint32_t pprev = get_partition_id(memblock[i-1].src);
			if(psrc != pprev) {
				row_offsets[psrc] = i;
			} 
		}


		rdtscll(load_stop);

		printf("# Offset big  time - by src %f \n" , (float)(load_stop-load_start) / (float)get_cpu_freq());			
		rdtscll(load_start);
		sequence::scanIBack(row_offsets,row_offsets,(int)P,minF<uintT>(),(uintT)nb_edges);   
		rdtscll(load_stop);
		printf("#Sequence big time - by src %f \n" , (float)(load_stop-load_start) / (float)get_cpu_freq());

		rdtscll(load_start);


		parallel_for(uint32_t i = 0; i < P; i++) {
			uint32_t o = row_offsets[i];
			//	if( o == nb_edges) continue;
			uint32_t l = ((i == P - 1) ? nb_edges - row_offsets[i] : row_offsets[i+1] - row_offsets[i]) ;  
			if (l != 0) 
				if(!sort_by_src) intSort::iSort(memblock + row_offsets[i],l, P + 1, getPartitionIdDst<struct edge_t>());
				else intSort::iSort(memblock + row_offsets[i],l, P + 1, getPartitionIdDst<struct edge_t>()); // getPartitionIdDst<struct edge_t>());
		}
		rdtscll(load_stop);
		printf("#Sort by dst time - by src %f \n" , (float)(load_stop-load_start) / (float)get_cpu_freq());
		parallel_for(uint32_t i = 0; i < P; i++) 
			for(uint32_t j = 0; j < P; j++) 
				offsets[i][j] = nb_edges  ;

		rdtscll(load_start);	
		for(uint32_t i = 0; i < P; i++) {
			uint32_t start, idx;
			offsets[i][0] = 0;
			uint32_t stop = ( i == P - 1?  nb_edges  : row_offsets[i + 1] ) - row_offsets[i];
			parallel_for(uint32_t start =  1;  start < stop; start++) {
				struct edge_t* e = &memblock[row_offsets[i] + start];
				uint32_t pdst = get_partition_id(e->dst);
				uint32_t pprev = get_partition_id(memblock[row_offsets[i] + start - 1].dst);
				if(pdst != pprev) {
					offsets[i][pdst] = start;
				}

			} 
		}


		rdtscll(load_stop);
		printf("#Offset dst time - by src %f \n" , (float)(load_stop-load_start) / (float)get_cpu_freq());
		rdtscll(load_start);
		parallel_for(uint32_t i = 0; i < P; i++) {
			sequence::scanIBack(offsets[i], offsets[i], (int)P, minF<uintT>(), (uintT)(i == P - 1 ? nb_edges - row_offsets[i] : row_offsets[i+1] - row_offsets[i]));
		}

		rdtscll(load_stop);
		printf("#Seq dst time - by src %f \n" , (float)(load_stop-load_start) / (float)get_cpu_freq());
		rdtscll(init_stop);
		printf("#Total grid create time %lu ( %f s)\n", init_stop - init_start, ((float)(init_stop - init_start))/(float)get_cpu_freq());
		if(sort_by_src){
			rdtscll(load_start);
			parallel_for(int i = 0; i < P; i++) {
				for(int j = 0; j < P; j ++) {
					uint32_t start = row_offsets[i] + offsets[i][j];
					uint32_t stop = (j == P-1 ? (i == P -1 ? nb_edges : row_offsets[i+1]) : row_offsets[i] + offsets[i][j+1] ) ;
					if(stop - start > 0)	
						intSort::iSort(memblock + start,stop - start, NB_NODES + 1, getEdgeSrc<struct edge_t>());
				}
			}
			rdtscll(load_stop);
			printf("#Sort cells by src time %lu ( %f s)\n", load_stop - load_start, ((float)(load_stop - load_start))/(float)get_cpu_freq());


		}
	}

}

void init_grid_nosort(int  full) {

	pthread_spinlock_t** b_locks;	
	uint32_t**  sizes;
	sizes = (uint32_t **) malloc(P * sizeof(uint32_t*));
	b_locks = (pthread_spinlock_t**) malloc(P * sizeof(pthread_spinlock_t*));
	offsets = (uint32_t**) malloc (P * sizeof(uint32_t*));
	blocks = (struct edge_t***) malloc(P * sizeof(struct edge_t**));



	for(size_t i =0; i< P; i++) {
		sizes[i] = (uint32_t* ) malloc(P*sizeof(uint32_t));
		offsets[i] = (uint32_t*) malloc(P* sizeof(uint32_t));
		b_locks[i] = (pthread_spinlock_t*) malloc( P * sizeof(pthread_spinlock_t));
		blocks[i] = (struct edge_t**) malloc(P * sizeof(struct edge_t*));
	}


	for(size_t i = 0; i < P; i++) {
		for(size_t j = 0; j < P; j++) {
			pthread_spin_init(&b_locks[i][j], PTHREAD_PROCESS_PRIVATE);
			offsets[i][j] = 0;
			uint32_t size = 4;
			sizes[i][j] = 16;
			blocks[i][j] = (struct edge_t* ) malloc(sizes[i][j]  * sizeof(struct edge_t));
		}  
	}


	rdtscll(init_start);
	if(not_processed) {

		parallel_for(size_t i = 0 ; i < NB_NODES; i++) nodes[i].nb_out_edges = 0;

		parallel_for(size_t index = 0; index < parallel_adds; index ++) {
			struct edge_t *i = &memblock[index]; //get_input(index);
			uint64_t src = i->src, dst = i->dst;

			uint32_t row = get_partition_id(src);
			uint32_t col = get_partition_id(dst);

			size_t tmpsize;
			pthread_spin_lock(&b_locks[row][col]);

			tmpsize = offsets[row][col]++;
			if(offsets[row][col] ==  sizes[row][col]) {
				sizes[row][col] *= 2;
				blocks[row][col] = (struct edge_t*) realloc(blocks[row][col], sizes[row][col] * sizeof(struct edge_t)); 
			}
			struct edge_t * e = &blocks[row][col][tmpsize];
			e->src = src;
			e->dst = dst;
			pthread_spin_unlock(&b_locks[row][col]);

		}
		rdtscll(init_stop);
		printf("#Total grid create time %lu ( %f s)\n", init_stop - init_start, ((float)(init_stop - init_start))/(float)get_cpu_freq());
	}


}
void count_degree() {
	uint64_t start_count,stop_count;
	rdtscll(start_count);

	parallel_for(size_t i  = 0; i < nb_edges; i++) {

		struct edge_t* e = &memblock[i];
		__sync_fetch_and_add(&nodes[e->src].nb_out_edges, 1);

		if(!isSymmetric) {
			__sync_fetch_and_add(&nodes[e->dst].nb_in_edges, 1);
		}


	}
	rdtscll(stop_count);
	printf("# Degree count time  %lu ( %f s)\n", stop_count - start_count, ((float)(stop_count - start_count))/(float)get_cpu_freq());

}
void init_adj_count(){

	uint64_t start_adj, stop_adj;
	if(!count_on_load) {

		count_degree();
	}

	edge_array_out = (struct edge*) malloc(nb_edges * sizeof(struct edge));
	if(!isSymmetric) edge_array_in = (struct edge*) malloc(nb_edges * sizeof(struct edge));

	rdtscll(start_adj);


	uint64_t current_offsets = 0; uint64_t curr_offsets_in = 0;

	// Compute offsets for the beginning of adjacency lists

	for(uint32_t i = 0; i < NB_NODES; i++) {
		nodes[i].outgoing_edges = current_offsets;
		current_offsets += nodes[i].nb_out_edges;
		nodes[i].nb_out_edges = 0;

		if(!isSymmetric)  {
			nodes[i].incoming_edges = curr_offsets_in;
			curr_offsets_in += nodes[i].nb_in_edges;
			nodes[i].nb_in_edges = 0;
		}
	}

	// Place edges in corresponding adjacency lists
	parallel_for(size_t i = 0; i < nb_edges;i++) {
		struct edge_t * e = &memblock[i];

		uint32_t src = e->src;
		uint32_t dst = e->dst;

		pthread_spin_lock(&nodes[src].lock);
		edge_array_out[nodes[src].outgoing_edges + nodes[src].nb_out_edges++].dst = dst;
		pthread_spin_unlock(&nodes[src].lock);

		if(!isSymmetric) {
			pthread_spin_lock(&nodes[dst].lock);
			edge_array_in[nodes[dst].incoming_edges + nodes[dst].nb_in_edges++].dst = src;
			pthread_spin_unlock(&nodes[dst].lock);
		}
	}

	rdtscll(stop_adj);
	printf("# Ad out edges create time  %lu ( %f s)\n", stop_adj - start_adj, ((float)(stop_adj - start_adj))/(float)get_cpu_freq());
	free(memblock); 

}
void init_adj_dynamic() {
	uint64_t start_adj, stop_adj;
	struct realloc_stats {
		uint32_t realloc_out;
		uint32_t realloc_in;
	};
	struct realloc_stats* thread_stats = (struct realloc_stats*) calloc(ALGO_NB_THREADS, sizeof(struct realloc_stats));	
	uint64_t edges_read = 0;
	uint64_t** lists = (uint64_t**) malloc ( NB_NODES * sizeof(uint64_t*));
	uint64_t* sizes = (uint64_t *) malloc(NB_NODES * sizeof(uint64_t));
	uint64_t** lists_in;
	uint64_t* sizes_in;
	if(!isSymmetric) {
		lists_in = (uint64_t**) malloc ( NB_NODES * sizeof(uint64_t*));
		sizes_in = (uint64_t *) malloc(NB_NODES * sizeof(uint64_t));

	}
	for(uint32_t i = 0; i < NB_NODES; i++) {
		sizes[i] = 4;
		lists[i] = (uint64_t*) malloc(4 * sizeof(uint64_t));
		nodes[i].nb_out_edges = 0;
		if(!isSymmetric) {
			sizes_in[i] = 4;
			lists_in[i] = (uint64_t*) malloc(4 * sizeof(uint64_t));
			nodes[i].nb_in_edges = 0;

		}
	}
	rdtscll(start_adj);

	parallel_for(size_t i = 0; i < nb_edges; i++) {
		//			struct input* in = get_input(i);
		struct edge_t* e = &memblock[i];
		uint32_t src = e->src;
		uint32_t dst = e->dst;

		pthread_spin_lock(&nodes[src].lock);
		int _id = get_threadid(); 
		if(nodes[src].nb_out_edges >= sizes[src]) {
			sizes[src] *= 2;
			lists[src] = (uint64_t*) realloc(lists[src], sizes[src] * sizeof(uint64_t));
			//			thread_stats[_id].realloc_out++;
		}
		lists[src][nodes[src].nb_out_edges++] = dst;
		pthread_spin_unlock(&nodes[src].lock);		 

		if(!isSymmetric) {
			pthread_spin_lock(&nodes[dst].lock);
			if(nodes[dst].nb_in_edges >= sizes_in[dst]) {
				sizes_in[dst] *= 2;
				lists_in[dst] = (uint64_t*) realloc(lists_in[dst], sizes_in[dst] * sizeof(uint64_t));
				//				thread_stats[_id].realloc_in++;

			}
			lists_in[dst][nodes[dst].nb_in_edges++] = src;
			pthread_spin_unlock(&nodes[dst].lock);
		}

	}
	rdtscll(stop_adj);
	//	for(uint32_t i =0; i < ALGO_NB_THREADS; i++) 
	//		printf("T%d::reallocs::out::%u::in::%u\n",i, thread_stats[i].realloc_out,thread_stats[i].realloc_in);
	printf("\n");
	printf("# Dynamic adj create  %lu ( %f s)\n", stop_adj - start_adj, ((float)(stop_adj - start_adj))/(float)get_cpu_freq());
	parallel_for(uint32_t i = 0; i < NB_NODES; i++) {
		free(lists[i]);
		if(!isSymmetric)	free(lists_in[i]);
	}
	free(lists);
	free(thread_stats);
	if(!isSymmetric) { free(lists_in); free(sizes_in); }
	free(sizes);


	free(memblock);
}
void init_edgelist(int full, int full_sort, int sort_edge_array) {

	uint64_t start_adj, stop_adj;

	if(!count_on_load) count_degree();


	rdtscll(start_adj);

	if(!sort_edge_array) return;
	intSort::iSort(memblock, nb_edges, NB_NODES + 1, getEdgeSrc<struct edge_t>());

	uint64_t* edge_offsets = (uint64_t*) malloc(NB_NODES * sizeof(uint64_t)); 
	parallel_for(uint32_t i = 0; i < NB_NODES; i++) edge_offsets[i] = nb_edges;

	//edge_array_out[0].dst = memblocs[0].dst;

	edge_offsets[memblock[0].src] = 0;
	if(memblock[0].src != 0) {
		uint32_t i = 0;
		while (i != memblock[0].src) 
			edge_offsets[i++] = 0;

	}
	parallel_for(size_t i = 1 ; i < nb_edges; i++) {
		struct edge_t* e = &memblock[i];
		e->dst = memblock[i].dst;
		if(memblock[i].src != memblock[i-1].src) 
			edge_offsets[memblock[i].src] = i;
	}
	sequence::scanIBack(edge_offsets, edge_offsets,(int)NB_NODES,minF<uintT>(),(uint64_t)(nb_edges));    

	parallel_for(uint32_t i = 0; i < NB_NODES; i++) {

		nodes[i].nb_out_edges = (i == NB_NODES-1 ? nb_edges : edge_offsets[i+1]) - edge_offsets[i]; 
		nodes[i].outgoing_edges = edge_offsets[i];
	}
	rdtscll(stop_adj);
	printf("# Ad out edges create time  %lu ( %f s)\n", stop_adj - start_adj, ((float)(stop_adj - start_adj))/(float)get_cpu_freq());
	//		intSort::iSort(memblock, nb_edges, NB_NODES + 1, getEdgeDst<struct edge_t>());

	intSort::iSort(memblock, nb_edges, NB_NODES + 1, getEdgeSrc<struct edge_t>());

	if(full_sort) {
		rdtscll(start_adj);
		parallel_for(uint32_t i = 0; i < NB_NODES; i++)
			if(nodes[i].nb_out_edges!= 0)
				intSort::iSort(memblock + nodes[i].outgoing_edges, nodes[i].nb_out_edges, NB_NODES + 1, getEdgeDst<struct edge_t>());
		rdtscll(stop_adj);
		printf("# Sort out lists by dst time  %lu ( %f s)\n", stop_adj - start_adj, ((float)(stop_adj - start_adj))/(float)get_cpu_freq());

	}
	free(edge_offsets);	
	//		free(memblock);

}
void init_adj_sort(int full, int full_sort) {

	uint32_t s,d;
	uint64_t start_adj, stop_adj;

	rdtscll(start_adj);

	if(createUndir) nb_edges *= 2;
	edge_array_out = (struct edge*) malloc(nb_edges * sizeof(struct edge));
	//		rdtscll(start_adj);

	intSort::iSort(memblock, nb_edges, NB_NODES + 1, getEdgeSrc<struct edge_t>());
	rdtscll(stop_adj);
	printf("# Sort time  %lu ( %f s)\n", stop_adj - start_adj, ((float)(stop_adj - start_adj))/(float)get_cpu_freq());
	rdtscll(start_adj);
#if CREATE_WEIGHT
	weights = (uint32_t*) malloc(nb_edges * sizeof(uint32_t));
	weights_in = (uint32_t *) malloc(nb_edges * sizeof(uint32_t));
#endif
	edge_array_out[0].dst = memblock[0].dst;
	uint64_t* edge_offsets = (uint64_t*) malloc(NB_NODES * sizeof(uint64_t)); 


	parallel_for(uint32_t i = 0; i < NB_NODES; i++) edge_offsets[i] = nb_edges;


	edge_array_out[0].dst = memblock[0].dst;
#if  WEIGHTED
	if(weighted_graph) edge_array_out[0].weight = memblock[0].weight;
	//	edge_offsets[0] = 0;
#endif

#if CREATE_WEIGHT

	s = memblock[0].src % 10;
	d =memblock[0].dst % 10;

	edge_array_out[0].weight = s+d;                 
#endif
	if(memblock[0].src != 0) {
		for(uint32_t i = 0; i < memblock[0].src; i++) {
			edge_offsets[i] = 0;
		}
	}	
	edge_offsets[memblock[0].src] = 0;

	parallel_for(uint64_t i = 1 ; i < nb_edges; i++) {
		struct edge* e = &edge_array_out[i];
#if WEIGHTED
		if(weighted_graph) e->weight = memblock[i].weight;
#endif
		e->dst = memblock[i].dst;
#if CREATE_WEIGHT
		uint32_t s1 = memblock[i].src % 10;
		uint32_t d1 =memblock[ i].dst % 10;

		edge_array_out[i].weight = s1+d1;
#endif
		if(memblock[i].src != memblock[i-1].src) 
			edge_offsets[memblock[i].src] = i;

	}
	sequence::scanIBack(edge_offsets, edge_offsets,(long)NB_NODES ,minF<uintT>(),(uint64_t)(nb_edges));    

	parallel_for(uint32_t i = 0; i < NB_NODES; i++) {

		nodes[i].nb_out_edges = (i == NB_NODES-1 ? nb_edges : edge_offsets[i+1]) - edge_offsets[i]; 
		nodes[i].outgoing_edges = edge_offsets[i];
	}
	rdtscll(stop_adj);
	//intSort::iSort(memblock, nb_edges, NB_NODES + 1, getEdgeDst<struct edge_t>());

	printf("# Ad out edges create time  %lu ( %f s)\n", stop_adj - start_adj, ((float)(stop_adj - start_adj))/(float)get_cpu_freq());

	if(full_sort) {
		rdtscll(start_adj);
		parallel_for(uint32_t i = 0; i < NB_NODES; i++) {
			if(nodes[i].nb_out_edges != 0)
				intSort::iSort(&(edge_array_out[nodes[i].outgoing_edges]), nodes[i].nb_out_edges, NB_NODES + 1, getEdgeDst<struct edge>());
		}
		rdtscll(stop_adj);
		printf("# Sort out lists by dst time  %lu ( %f s)\n", stop_adj - start_adj, ((float)(stop_adj - start_adj))/(float)get_cpu_freq());

	}

	if(!isSymmetric){

		rdtscll(stop_adj);
		printf("# Load time  %lu ( %f s)\n", stop_adj - start_adj, ((float)(stop_adj - start_adj))/(float)get_cpu_freq());

		edge_array_in = (struct edge*) malloc(nb_edges * sizeof(struct edge));
		rdtscll(start_adj);

		intSort::iSort(memblock, nb_edges, NB_NODES + 1, getEdgeDst<struct edge_t>());

		edge_array_in[0].dst = memblock[0].dst;
		uint64_t* edge_offsets_in = (uint64_t*) malloc(NB_NODES * sizeof(uint64_t)); 
		parallel_for(uint32_t i = 0; i < NB_NODES; i++) edge_offsets_in[i] = nb_edges;

		edge_array_in[0].dst = memblock[0].src;
#if WEIGHTED
		if(weighted_graph) edge_array_in[0].weight = memblock[0].weight;
#endif
		edge_offsets_in[memblock[0].dst] = 0;

#if CREATE_WEIGHT
		s = memblock[0].src % 10; 
		d =memblock[0].dst % 10;

		edge_array_in[0].weight = s+d ; //(NB_NODES > s ? NB_NODES - s : s  - NB_NODES) % 10;
#endif
		if(memblock[0].dst != 0) {
			for(uint32_t i = 0; i < memblock[0].dst; i++)
				edge_offsets_in[i] = 0;

		}
		parallel_for(uint64_t i = 1 ; i < nb_edges; i++) {
			struct edge* e = &edge_array_in[i];
			e->dst = memblock[i].src;
#if WEIGHTED
			if(weighted_graph) e->weight = memblock[i].weight;
#endif
#if CREATE_WEIGHT
			uint32_t s1 = memblock[i].src % 10; 
			uint32_t d1 =memblock[ i].dst % 10; 

			edge_array_in[i].weight = s1+d1;		

#endif
			if(memblock[i].dst != memblock[i-1].dst) 
				edge_offsets_in[memblock[i].dst] = i;
		}
		sequence::scanIBack(edge_offsets_in, edge_offsets_in,(int)NB_NODES ,minF<uintT>(),(uint64_t)(nb_edges));    

		parallel_for(uint32_t i = 0; i < NB_NODES; i++) {

			nodes[i].nb_in_edges = (i == NB_NODES-1 ? nb_edges : edge_offsets_in[i+1]) - edge_offsets_in[i]; 
			nodes[i].incoming_edges = edge_offsets_in[i];
		}
		rdtscll(stop_adj);
		//intSort::iSort(memblock, nb_edges, NB_NODES + 1, getEdgeDst<struct edge_t>());


		printf("# Ad in edges create time  %lu ( %f s)\n", stop_adj - start_adj, ((float)(stop_adj - start_adj))/(float)get_cpu_freq());

		if(full_sort) {
			rdtscll(start_adj);
			parallel_for(uint32_t i = 0; i < NB_NODES; i++) {
				if(nodes[i].nb_in_edges != 0)
					intSort::iSort(&(edge_array_in[nodes[i].incoming_edges]), nodes[i].nb_in_edges, NB_NODES + 1, getEdgeDst<struct edge>());
			}
			rdtscll(stop_adj);
			printf("# Sort out lists by dst time  %lu ( %f s)\n", stop_adj - start_adj, ((float)(stop_adj - start_adj))/(float)get_cpu_freq());

		}
	}

	free(edge_offsets);	
	free(memblock);


}
void init(int full) {
	parallel_for(int i = 0;  i  < ALGO_NB_THREADS;i++) {

              int id = get_threadid();
                 cpu_set_t mask;
                CPU_ZERO(&mask);
                CPU_SET(id, &mask);
                sched_setaffinity(gettid(), sizeof(mask), &mask);

        }
	nodes = (struct node*) malloc(NB_NODES*sizeof(*nodes));
	locks = (pthread_spinlock_t* ) malloc(NB_NODES * sizeof(pthread_spinlock_t));

	parallel_for(size_t i = 0; i < NB_NODES; i++) {
		memset(&nodes[i], 0, sizeof(nodes[i]));
		pthread_spin_init(&locks[i], PTHREAD_PROCESS_PRIVATE);
		pthread_spin_init(&nodes[i].lock, PTHREAD_PROCESS_PRIVATE);
	}

	printf("Load start");

	struct stat sb;
	int fd = open(filename, O_RDONLY);
	if(!fd == -1)
		die("Cannot open %s\n", filename);
	fstat(fd, &sb);
	printf("#Size: %lu\n", (uint64_t)sb.st_size);
	assert(fd);
	memblock1 = (char*) mmap(NULL, sb.st_size, PROT_WRITE, MAP_PRIVATE, fd, 0);
	if(not_processed) {
		nb_edges = sb.st_size/item_size();
	}
	else {
		nb_edges = (sb.st_size - 2*sizeof(uint64_t) * NB_NODES)/2/item_size();
	}
	if(createUndir)
		memblock = (struct edge_t*) malloc( 2* nb_edges * sizeof(struct edge_t));
	else
		memblock = (struct edge_t*) malloc( nb_edges * sizeof(struct edge_t));
	row_offsets = (uint32_t*) calloc(P, sizeof(uint32_t));

	printf("#NB edges = %lu\n", nb_edges);

	parallel_adds = nb_edges;
	if(full)
		parallel_adds = nb_edges;

	mod_val = NB_NODES % P;
	if(mod_val != 0)
		partition_size = NB_NODES / P + 1;
	else partition_size = NB_NODES / P;
	split_point = mod_val * partition_size;

	//Preload graph into memory
	uint64_t load_start,load_stop;
	rdtscll(load_start);

	preload_graph();
	rdtscll(load_stop);
	printf("#Load time %lu ( %f s )\n", load_stop - load_start, ((float)(load_stop - load_start))/(float)get_cpu_freq());
	if(load_mode < 2 || load_mode >= 6) 
		count_degree(); //For all these modes it can be optimized so that we use the fast scan operation	


	switch(load_mode) {

		case 0: //grid sort 
			init_grid_sort(0, 0);
	
			break;
		case 1: //grid nosort
			init_grid_nosort(full);
			break;
		case 2: //adj created with sort
			init_adj_sort(0, 0);
			break;
		case 3: //adj created with sort and fully sorted within 
			init_adj_sort(full, 0);
			break;

		case 4: //adj list created on load; arg = 0 -> dynamic; arg = 1 -> count sort 
			init_adj_count(); //full);
			break;
		case 5: init_adj_dynamic();
		 break;
		case 6: //grid using sorting but cells are sorted by src (not dst)
			init_grid_sort(full, 1);
			break;
		case 7: //edgelista is adjacency list (sorted by source)
			init_edgelist(full, 0, 1);
			break;
		case 8: //edge list is sorted by source and then by destination within edges having the same source
			init_edgelist(full, 1, 1);
			break;
		case 9: init_edgelist(full, 0, 0);


	}
	//	exit(1);

	munmap(memblock1, sb.st_size);

}


