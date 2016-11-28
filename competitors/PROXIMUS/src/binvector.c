/********************************************************************************
* Authors     : Mehmet Koyuturk & Ananth Grama                                  *
* Institution : Purdue University, Department of Computer Sciences              *
* Last update : 04/30/2004                                                      *
* Copyright   : See README file for copyright information.                      *
********************************************************************************/

#include "../include/binvector.h"

BinVector *
initBinVector(unsigned int n){
  BinVector *x = (BinVector *) tmalloc (sizeof(BinVector));
  x->n         = n;
  x->nz        = 0;
  x->list      = (unsigned int *) tmalloc ((inWord(n)+1)*sizeof(unsigned int));
  return x;
}

void
freeBinVector(BinVector *x){
  tfree (x->list, (inWord(x->n)+1)*sizeof(unsigned int));
  tfree (x, sizeof(BinVector)); 
}

unsigned int 
element(BinVector *x, unsigned int i){
  return((x->list[inWord(i)]&setBit(bitOff(i)))>0);
}

unsigned int 
setBit(unsigned int i){
  unsigned int z = ONE;
  return (z << i);
}

unsigned int 
noOfOnes(unsigned int z){
  if (z < 2)
    return z;
  else
    return ((z%2)+noOfOnes(z>>1));
}

void
displayBinVector(BinVector *x, char *name){
  unsigned int i;
  mprintf("%s[%d, %d] ", name, x->n, x->nz); 
  for (i=0; i<x->n; i++){
    if (element(x,i))
      mprintf("%d ",i);
  }
  mprintf("\n");
}
