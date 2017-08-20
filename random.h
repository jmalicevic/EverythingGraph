#ifndef RANDOM_H
#define RANDOM_H


#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <numa.h>

#include "parallel_ligra.h"
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/time.h>
#include <stdint.h>
#include <locale.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <pthread.h>
#include <omp.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <sys/sysinfo.h>
#include <sys/syscall.h>
#include <sched.h>
#include <signal.h>

#define NUMA 1
#define NAIVE_SPLIT 0   // Do smarter load repartition if 0
#define WEIGHTED 0 // set flag to 1  
extern int silent;
extern int isSymmetric;
extern int createUndir;
extern int not_processed;
extern int load_mode;
extern int numa;

// Graph representation
struct node;
struct edge;
struct thread_buffer;

#include ALGO_H
#define GRID 0 
#if WEIGHTED 
#define weight_t typeof(((struct edge*)NULL)->weight)
#else
#define weight_t float

#endif
#define ALGO_NB_THREADS 16   // number of threads for the algorithm
#define NB_CONCURRENCY 16 // number of concurrent threads adding/removing edges

#if defined(SSSPADJ_H) 
#define CREATE_WEIGHT 1 //THe graphs used were unweighted so we used this to create them on the fly. Remove if graph is weighted
extern uint32_t* weights;
extern uint32_t* weights_in;
#endif

struct algo_func {
   void (*reset)(struct node *);
   void (*main)(struct node *);
     void (*construct)(void);
   void (*destruct)(void);
};

struct edge_buffer {
   uint32_t nb_edges;
   uint32_t nb_max_edges;
   struct edge *edges;
};

struct edge_buffer_list {
   uint32_t nb_edges;
   uint32_t nb_buffers;
   struct edge_buffer *buffers;
};
extern struct edge* edge_array_out;
extern struct edge* edge_array_in;


extern pthread_spinlock_t* locks;
struct node {
 uint32_t nb_out_edges;
  uint64_t outgoing_edges;
  uint32_t nb_in_edges;
  uint64_t incoming_edges;
  pthread_spinlock_t lock;
};


typedef enum { ALGO, BEST } algo_t;
extern algo_t algo_phase;
extern volatile int parallels_done;
extern volatile int waiting;
extern __thread int __id;

extern struct algo_func current_algo;
#include "init_all.h"
#include "utils.h"
#include "buffer.h"

static inline used uint32_t id(struct node *n) {
   return (uint32_t)(long)(n - &nodes[0]);
}

extern uint32_t items_in_frontier;
extern uint32_t front_degree;
extern char* in_frontier;
extern char* in_frontier_next;
typedef enum {
		PUSH,
		PULL
} algo_mode;
extern algo_mode mode;
extern int switch_mode; //should we switch between push and pull
#endif
