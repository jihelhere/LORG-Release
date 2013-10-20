#include "dual.h"
#include "tools.h"

#include <stdlib.h>

void dual_init_vectors(dual_t *d){
  uint32_t i,t,t2;

   for(i = 0 ; i < d-> N ; ++i)
   {
     for(t = 0 ; t < d->Y ; ++t)
     {
       d->unary_u[i*(d->N)+t] = 0.0;
       for(t2 = 0 ; t2 < d->Y ; ++t2)
       {
         d->binary_u[i*(d->N)*t*(d->Y) + t2] = 0.0;
       }
     }
   }
}



dual_t* dual_init(uint32_t slength, uint32_t nlabels)
{
  dual_t * d = (dual_t *)malloc(sizeof(dual_t));
  d->N = slength;
  d->Y = nlabels;
  d->unary_u = (double *)malloc(sizeof(double)*slength*nlabels);
  if(d->unary_u == NULL)
  {
    fatal("memory error");
  }
  d->binary_u = (double*)malloc(sizeof(double)*slength*nlabels*nlabels);
  if(d->binary_u == NULL)
  {
    free(d->unary_u);
    fatal("memory error");
  }

  dual_init_vectors(d);
  return d;
}


double dual_get_unary_penalty(dual_t *d, uint32_t i,uint32_t t)
{
  uint32_t index = i * (d->N) + t;
   return (d->unary_u[index]);
}

void dual_set_unary_penalty(dual_t *d, uint32_t i,uint32_t t,double val)
{
     d->unary_u[i*(d->N)+t] = val;
}

void dual_add_unary_penalty(dual_t *d, uint32_t i,uint32_t t,double val)
{
     d->unary_u[i*(d->N)+t] += val;
}


void dual_free(dual_t *d)
{
  if(d != NULL)
  {
    if(d->unary_u != NULL)
    {
      free(d->unary_u);
    }
    if(d->binary_u != NULL)
    {
      free(d->binary_u);
    }
    free(d);
  }
}


double dual_get_binary_penalty(dual_t *d, uint32_t i,uint32_t t1, uint32_t t2)
{
  uint32_t index = i * (d->N) * t1 * (d->Y) + t2;
  return (d->binary_u[index]);
}

void dual_set_binary_penalty(dual_t *d, uint32_t i,uint32_t t1, uint32_t t2, double val)
{
  d->binary_u[i*(d->N)*t1*(d->Y)+t2] = val;
}

void dual_add_binary_penalty(dual_t *d, uint32_t i,uint32_t t1, uint32_t t2, double val)
{
  d->binary_u[i*(d->N)*t1*(d->Y)+t2] += val;
}
