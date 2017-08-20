//
// Created by JaSmiNa on 04/10/16.
//
#pragma once
#ifndef PRGRID_H
#define PRGRID_H

void prgrid_construct(void);
void prgrid_destruct(void);
void prgrid_reset(struct node *nodes);
void prgrid(struct node *nodes);
void prgrid_rerun(struct node *nodes);

struct edge {
    uint32_t dst;
#if WEIGHTED
    float weight;
#endif
};

extern struct node* node_list;
extern float* prev;
extern float* rank;



#endif //GRAPH_UP_BFS_LIGRA_H
