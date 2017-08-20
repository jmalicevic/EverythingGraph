//
// Created by JaSmiNa on 04/10/16.
//
#pragma once
#ifndef PRSIMPLE_H
#define PRSIMPLE_H

void pr_construct(void);
void pr_destruct(void);
void pr_reset(struct node *nodes);
void pr(struct node *nodes);
void pr_rerun(struct node *nodes);

struct edge {
    uint32_t dst;
#if WEIGHTED
    float weight;
#endif
};

extern struct node* node_list;
extern float* prev;
extern float* rank;



#endif //
