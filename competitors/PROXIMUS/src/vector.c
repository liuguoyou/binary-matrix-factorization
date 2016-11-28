/********************************************************************************
* Author      : Mehmet Koyuturk & Ananth Grama                                  *
* Institution : Purdue University, Department of Computer Sciences              *
* Last update : 03/30/2004                                                      *
* Copyright   : See README file for copyright information.                      *
********************************************************************************/

#include "../include/vector.h"

Vector *
initVector(unsigned int n){
  Vector *s = (Vector *) tmalloc (sizeof(Vector));
  s->n      = n;
  s->max    = 0;
  s->list   = (unsigned int *) tmalloc (n*sizeof(unsigned int));
  return s;
}

void
freeVector(Vector *s){
  tfree(s->list, s->n*sizeof(unsigned int));
  tfree(s, sizeof(Vector));
}

unsigned int *
countingSort(Vector *s){
  unsigned int *index, *count;
  unsigned int i, v;
  unsigned int max = s->max+1;  // bound on elements in the array 

  index = (unsigned int *) tmalloc (s->n*sizeof(unsigned int));
  count = (unsigned int *) tmalloc (max*sizeof(unsigned int));
  memset(count, 0, max*sizeof(unsigned int));
  
  // Count occurence of each possible value in the array
  for (i=0; i<s->n; i++)
    count[s->list[i]]++;
  
  // Compute prefix sum
  for (i=1; i< max; i++)
    count[i] += count[i-1];
  
  // Locate each value in the sorted array
  for(i=0; i<s->n; i++){
    v = s->list[i];
    count[v]--;
    index[count[v]] = i;
  }
  // index now contains position of each entry in the sorted array

  tfree(count, max*sizeof(unsigned int));
  
  return index;
}


double
continuousSolveXS(BinVector *x, Vector *s){
  unsigned int *index;
  int i,n,k;
  Bool stop;
  unsigned int nw, sum;
  double obj, nobj;

  if ( x->n != s->n)
    errexit("solveXS: Vector sizes don't match x:%d s%d\n", x->n, s->n);

  n = x->n;
  nw = inWord(n)+1;

  index = countingSort(s); // First sort s

  memset(x->list, 0, nw*sizeof(unsigned int));

  sum   = 0;
  obj   = 0.0;
  x->nz = 0;

  for(i=n-1, stop=FALSE; i>=0 && !stop; i--){
    k = index[i];                     // Visiting each entry of x in sorted or
    sum += s->list[k];
    nobj = (1.0*sum*sum)/((n-i));     // Check if improves the objective function
    if ( nobj > obj ){                // if so, set that entry of x
      x->list[inWord(k)] += setBit(bitOff(k));
      x->nz++;
      obj = nobj;
    }
    else                              // else, stop
      stop = TRUE;
  }

  tfree(index, s->n*sizeof(unsigned int));

  return(obj);
}


int
discreteSolveXS(BinVector *x, Vector *s, unsigned int y){
  int i,k;
  Bool stop;
  unsigned int nw;
  int sum, max;
  int argmax;

  if ( x->n != s->n)
    errexit("solveXS: Vector sizes don't match x:%d s%d\n", x->n, s->n);

  nw = inWord(s->n)+1;

  memset(x->list, 0, nw*sizeof(unsigned int));
  sum= 0;
  max = 0;
  x->nz = 0;

  for(i=0; i<x->n ; i++){
    if ((s->list[i]+s->list[i]) >= y){
      x->list[inWord(i)] += setBit(bitOff(i));
      sum += s->list[i]+s->list[i]-y;           
      x->nz++;
    }
    if (s->list[i] > max){
      max = s->list[i];
      argmax = i;
    }
  }

  if ( !(x->nz) ){
    x->list[inWord(argmax)] += setBit(bitOff(argmax));
    x->nz = 1;
    sum = max;
  }

  return(sum);
}

void
displayVector(Vector *s){
  unsigned int i;
  for (i=0; i<s->n; i++){
    if (s->list[i] > 0)
      mprintf("%d:%d ",i,s->list[i]);
  }
  mprintf("\n");
}
