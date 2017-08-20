#include "random.h"

uint64_t get_cpu_freq(void) {
   FILE *fd;
   uint64_t freq = 0;
   float freqf = 0;
   char *line = NULL;
   size_t len = 0;

   fd = fopen("/proc/cpuinfo", "r");
   if (!fd) {
      fprintf(stderr, "failed to get cpu frequency\n");
      perror(NULL);
      return freq;
   }

   while (getline(&line, &len, fd) != EOF) {
      if (sscanf(line, "cpu MHz\t: %f", &freqf) == 1) {
         freqf = freqf * 1000000UL;
         freq = (uint64_t) freqf;
         break;
      }
   }

   fclose(fd);
   return freq;
}

uint64_t fast_atollu(const char * str, size_t *new_index) {
   uint64_t val = 0;
   const char *beginning = str;
   while(1) {
      char c = *str;
      if(c > '9' || c < '0')
         break;
      val = val*10LLU + (c - '0');
      str++;
   }
   *new_index += (str - beginning);
   return val;
}

int sort_uint(const void *_a, const void *_b) {
   const uint64_t a = *(uint64_t*)_a, b = *(uint64_t*)_b;
   if(a > b)
      return 1;
   if(b > a)
      return -1;
   return 0;
}

int sort_edges_dst(const void *a, const void *b) {
   const struct edge * const e1 = (const struct edge*)a, * const e2 = (const struct edge*) b;
   if(e1->dst > e2->dst)
      return 1;
   if(e1->dst < e2->dst)
      return -1;
   return 0;
}

int sort_edges_src_dst(const void *a, const void *b) {
   const struct edge * const e1 =(const struct edge*) a, * const e2 =(const struct edge*) b;
#if HAS_EDGE_SRC
   if(e1->src > e2->src)
      return 1;
   if(e1->src < e2->src)
      return -1;
#endif
   if(e1->dst > e2->dst)
      return 1;
   if(e1->dst < e2->dst)
      return -1;
   return 0;
}

uint32_t lock_edge_buffer(struct edge_buffer *buffer) {
begin:
   fence();
   uint32_t nb_max = *(volatile uint32_t*)&buffer->nb_max_edges;
   fence();
   return nb_max;
}

int unlock_edge_buffer(struct edge_buffer *buffer, uint32_t nb_max) {
   fence();
   buffer->nb_max_edges = nb_max; //unlocks the node
   return 0;
}

struct node_list get_all_nodes(void) {
      size_t i;
      uint64_t nb_starting_nodes;
      uint64_t *starting_nodes = (uint64_t*) malloc(NB_NODES * sizeof(*starting_nodes));
#pragma omp parallel for
      for(i = 0; i < NB_NODES; i++)
         starting_nodes[i] = i;
      nb_starting_nodes = NB_NODES;
      struct node_list n = { .nb_starting_nodes = nb_starting_nodes, .starting_nodes = starting_nodes };
      return n;
}

int sort_edges(const void *_a, const void *_b) {
   const struct edge *a = (const struct edge*)_a, *b = (const struct edge*) _b;
   if(a->dst > b->dst)
      return 1;
   else if(a->dst < b->dst)
      return -1;
   else
      return 0;
}

void sort_buffer(struct edge_buffer *b) {
   qsort(b->edges, b->nb_edges, sizeof(*b->edges), sort_edges);
}

void print_utils_stats(void) {
}
