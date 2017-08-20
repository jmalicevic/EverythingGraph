//
// Created by JaSmiNa on 04/10/16.
//

#ifndef SSSPADJ_H
#define SSSPADJ_H

void sssp_construct(void);
void sssp_destruct(void);
void sssp_reset(struct node *nodes);
void sssp(struct node *nodes);
struct edge {
    uint32_t dst;
    float weight;
};

extern uint32_t* dist;




#endif 
