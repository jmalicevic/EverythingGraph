#pragma once
#ifndef BFSNUMA_H
#define BFSNUMA_H

void bfsnuma_construct(void);
void bfsnuma_destruct(void);
void bfsnuma(struct node *nodes);
void bfsnuma_reset(struct node *nodes);

struct edge{
	uint32_t dst;
#if WEIGHTED
	float weight;
#endif
};

extern uint32_t* parent;
 
#endif
