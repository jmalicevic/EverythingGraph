#ifndef UTILS_H
#define UTILS_H

#include "random.h"

/*
 * Addition and deletion functions
 */
struct node_list {
   uint64_t nb_starting_nodes;
   uint64_t *starting_nodes;
};
struct node_list get_all_nodes(void);

void print_utils_stats(void);

/*
 * Misc
 */
uint64_t get_cpu_freq(void);
void shuffle(size_t *array, size_t n);
int sort_uint(const void *_a, const void *_b);

uint64_t fast_atollu(const char * str, size_t *new_index);
#define NB_MAX_EDGES -1UL 

/*
 * Macros to iterate over nodes or edges
 */
#define GET_MACRO(_1,_2,_3,NAME,...) NAME

#define foreach_node(l, dst) \
   for(size_t ___i = 0; ___i < l.nb_starting_nodes; ___i++) \
      if((dst = &nodes[l.starting_nodes[___i]]))


uint32_t lock_edge_buffer(struct edge_buffer *buffer);
int unlock_edge_buffer(struct edge_buffer *buffer, uint32_t nb_max);

#define foreach_outgoing_edges2(b, _dst) \
   for( size_t ___j = 0,  ___size = n->nb_out_edges;  ( ___size <= NB_MAX_EDGES && ___j < ___size) /*|| unlock_edge_buffer( &n->outgoing_edges, ___max)*/; ___j++) \
      for( size_t ___k = 0; ( ___size <= NB_MAX_EDGES && ___k < 1) ; ___k++) \
         if(  ___size <= NB_MAX_EDGES && ( /*n->outgoing_edges.edges[___j]*/edge_array_out[n->outgoing_edges+___j].dst < NB_NODES) && ( _dst = &nodes[edge_array_out[n->outgoing_edges+___j].dst/*n->outgoing_edges.edges[___j].dst*/])   )

#define foreach_outgoing_edges3(b, _dst, e) \
   for( size_t ___j = 0, /*___max = lock_edge_buffer( &n->outgoing_edges),*/ ___size = n->nb_out_edges;/*outgoing_edges.nb_edges;*/ ( ___size <= NB_MAX_EDGES && ___j < ___size) /* || unlock_edge_buffer( &n->outgoing_edges, ___max) */; ___j++) \
      for (size_t ___k = 0; (___size <= NB_MAX_EDGES && ___k < 1); ___k++) \
         if( ___size <= NB_MAX_EDGES && (edge_array_out[n->outgoing_edges+___j].dst < NB_NODES) && (_dst = &nodes[edge_array_out[n->outgoing_edges + ___j].dst]) && (e = &edge_array_out[n->outgoing_edges+___j])  )

#define foreach_outgoing_edges(...) \
   GET_MACRO(__VA_ARGS__, foreach_outgoing_edges3, foreach_outgoing_edges2)(__VA_ARGS__)

#define foreach_incoming_edges2(b, _dst) \
   for( size_t ___j = 0, /* ___max = lock_edge_buffer(&n->incoming_edges), */ ___size = n->nb_in_edges; (___size <= NB_MAX_EDGES && ___j < ___size) /*  || unlock_edge_buffer(&n->incoming_edges, ___max) */; ___j++) \
      for(size_t ___k = 0; (___size <= NB_MAX_EDGES && ___k < 1)  ; ___k++) \
         if(  ___size <= NB_MAX_EDGES && (edge_array_in[n->incoming_edges+___j].dst < NB_NODES) && (_dst = &nodes[edge_array_in[n->incoming_edges+___j].dst]) )

#define foreach_incoming_edges3(b, _dst, e) \
   for(size_t ___j = 0, /* ___max = lock_edge_buffer(&n->incoming_edges), */ ___size = n->nb_in_edges;  (___size <= NB_MAX_EDGES && ___j < ___size)  /*|| unlock_edge_buffer(&n->incoming_edges, ___max)*/; ___j++) \
      for(size_t ___k = 0; (___size <= NB_MAX_EDGES && ___k < 1) ; ___k++) \
         if( ___size <= NB_MAX_EDGES && (edge_array_in[n->incoming_edges+___j].dst < NB_NODES) && (_dst = &nodes[edge_array_in[n->incoming_edges+___j].dst]) && (e = &edge_array_in[n->incoming_edges + ___j])  )

#define foreach_incoming_edges(...) \
   GET_MACRO(__VA_ARGS__, foreach_incoming_edges3, foreach_incoming_edges2)(__VA_ARGS__)

#define foreach_edges(n, _dst, code) \
   foreach_outgoing_edges(n, _dst) { code; } \
   foreach_incoming_edges(n, _dst) { code; }


/*
 * For most algorithms, we only update a node if its value meets a condition.
 * Check that optimistically and then lock if the condition is met and recheck.
 */
#define LOCKED_IF(dst, cond) \
   if(cond) { \
      pthread_spin_lock(&(dst)->lock); \
      if(!(cond)) { pthread_spin_unlock(&(dst)->lock); } else

#define LOCKED_IF_END }


/*
 * Misc
 */
#ifdef __x86_64__
#define rdtscll(val) {                                           \
       unsigned int __a,__d;                                        \
       asm volatile("rdtsc" : "=a" (__a), "=d" (__d));              \
       (val) = ((unsigned long)__a) | (((unsigned long)__d)<<32);   \
}
#else
#define rdtscll(val) __asm__ __volatile__("rdtsc" : "=A" (val))
#endif

#define used __attribute__((unused))

#define die(msg, args...) \
do {                         \
            fprintf(stderr,"(%s,%d) " msg "\n", __FUNCTION__ , __LINE__, ##args); \
            assert(0);                 \
         } while(0)


static inline pid_t gettid(void) {
   return syscall(__NR_gettid);
}
#define NOP10() asm("nop;nop;nop;nop;nop;nop;nop;nop;nop;nop")

#define fence() asm volatile ("" : : : "memory")
#endif
