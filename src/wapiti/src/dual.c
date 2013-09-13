#include "dual.h"
#include "tools.h"

#include <stdlib.h>

void dual_init_vectors(dual_t *d){
   uint32_t i,t;
   for(i = 0 ; i < d-> N ; i++){
      for(t = 0 ; t < d->Y ; t++){
         d->u[i*(d->Y)+t] = 0.0;
         //d->z[i*(d->Y)+t] = 0;
      }
   }
}



dual_t* dual_init(uint32_t slength, uint32_t nlabels){
    dual_t * d = (dual_t *)malloc(sizeof(dual_t));
    d->N = slength;
    d->Y = nlabels;
    d->u = (double *)malloc(sizeof(double)*slength*nlabels);
    if(d->u == NULL){
      fatal("memory error");
    }
 /*   d->z = (int *)malloc(sizeof(int)*slength*nlabels);
    if(d->z == NULL){
      fatal("memory error");
    }*/
    dual_init_vectors(d);
    return d;
}


double dual_get_penalty(dual_t *d, uint32_t i,uint32_t t){
   uint32_t index = i * (d->Y) + t;
   //return (d->u[index])*(d->z[index]);
   return (d->u[index]);
}

void dual_set_penalty(dual_t *d, uint32_t i,uint32_t t,double val){
     d->u[i*(d->Y)+t] = val;
}

void dual_add_penalty(dual_t *d, uint32_t i,uint32_t t,double val){
     d->u[i*(d->Y)+t] += val;
}

/*
void dual_set_in_z(dual_t *d, uint32_t i,uint32_t t){
     d->z[i*(d->Y)+t] = 1;
}*/


void dual_free(dual_t *d){
   if(d != NULL){
      if(d->u != NULL){
        free(d->u);
      }
 /*    if(d->z != NULL){
        free(d->z);
      }*/
      free(d);
   }
}
