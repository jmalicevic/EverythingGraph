#ifndef SPMV_H
#define SPMV_H

void spmv_construct(void);
void spmv_destruct(void);
void spmv_reset(struct node *nodes);
void spmv(struct node *nodes);

extern uint32_t* value_in;
extern uint32_t* value_out;
struct edge {
   uint32_t dst;
#if WEIGHTED
   float weight;
#endif
};

#endif
