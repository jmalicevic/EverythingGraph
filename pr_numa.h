#pragma once
#ifndef PRNUMA_H
#define PRNUMA_H

void prnuma_construct(void);
void prnuma_destruct(void);
void prnuma(struct node *nodes);
void prnuma_reset(struct node *nodes);

struct edge{
	uint32_t dst;
#if WEIGHTED
	float weight;
#endif
};

extern float* prev_rank;
extern float* curr_rank;
extern uint32_t* degree;
#endif
