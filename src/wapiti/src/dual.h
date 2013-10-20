#ifndef DUAL_H
#define DUAL_H

#include <stdint.h>


struct dual_s
{
  uint32_t N;  //sentence length
  uint32_t Y;  //label count
  double* unary_u;   //penalization vector: each element is a
                     //penalization and depends on position i and tag
                     //t (represented by an int)
  double* binary_u; //binary penalties indexed by position and tag1, tag2
};
typedef struct dual_s dual_t;


dual_t* dual_init(uint32_t, uint32_t);
void dual_set_unary_penalty(dual_t *, uint32_t, uint32_t, double);
void dual_add_unary_penalty(dual_t *, uint32_t, uint32_t, double);
double dual_get_unary_penalty(dual_t *, uint32_t, uint32_t);

void dual_set_binary_penalty(dual_t *, uint32_t, uint32_t, uint32_t, double);
void dual_add_binary_penalty(dual_t *, uint32_t, uint32_t, uint32_t, double);
double dual_get_binary_penalty(dual_t *, uint32_t, uint32_t, uint32_t);

#endif
