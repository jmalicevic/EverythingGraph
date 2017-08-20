#ifndef WCC_H
#define WCC_H

void wcc_construct(void);
void wcc_destruct(void);
void wcc_reset(struct node *nodes);
void wcc(struct node *nodes);
void *wcc_parallels(void* id);

struct edge {
   uint32_t dst;
#if WEIGHTED
   uint32_t weight;
#endif
};

struct metadata {
   uint64_t component;
   struct node *father;
};

#endif
