#ifndef DUAL_H
#define DUAL_H

#include <stdint.h>


struct dual_s{
   uint32_t N;  //sentence length
   uint32_t Y;  //label count
   double* u;   //penalization vector: each element is a penalization and depends on position i and tag t (represented by an int) 
 //  int* z; //boolean vector: each value indicates if the tag t is the tag provided by analyser at position i
};
typedef struct dual_s dual_t;


dual_t* dual_init(uint32_t, uint32_t);
void dual_set_penalty(dual_t *, uint32_t, uint32_t, double);
void dual_add_penalty(dual_t *, uint32_t, uint32_t, double);
//void dual_set_in_z(dual_t *, uint32_t, uint32_t);
double dual_get_penalty(dual_t *, uint32_t, uint32_t);

#endif
