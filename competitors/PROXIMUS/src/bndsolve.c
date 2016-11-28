/********************************************************************************
* Authors     : Mehmet Koyuturk & Ananth Grama                                  *
* Institution : Purdue University, Department of Computer Sciences              *
* Last update : 03/31/2004                                                      *
* Copyright   : See README file for copyright information.                      *
********************************************************************************/

#include "../include/bndsolve.h"

void
discreteRank1App(Matrix *A, BinVector *x, BinVector *y, Byte init){
  Vector *s, *z;                    // The vectors Ay and A^Tx respectively
  unsigned int i,j;
  int obj, nobj;
  Row *nrow;

  s = initVector(A->m);
  z = initVector(A->n);

  // Initialization of pattern vector
  if (init != GGG && init!=NEIGHBOR)
    initialize(A, y, init);
  else{ // Initialization targets presence vector for these methods
    initialize(A, x, init);
    lMVMult(A, x, z);
    discreteSolveXS(y, z, x->nz);
  }

  obj     = -1;
  nobj    = 0;
  for (j=0; j<MAXIT && nobj>obj; j++){
    obj = nobj;
    rMVMult(A, y, s);
    discreteSolveXS(x, s, y->nz);             // Solve for x
    lMVMult(A, x, z);
    nobj    = discreteSolveXS(y, z, x->nz);  // Solve for y
  }

  // Record presence into matrix
  A->ones  = 0;
  for (nrow=A->rlist, i=0; nrow != ((Row *)(NULLROW)); nrow=nrow->next, i++){
    nrow->sel = element(x, i);
    A->ones  += nrow->sel;
  }
  A->zeros = A->m-A->ones;

  freeVector(z);
  freeVector(s);
}

void
continuousRank1App(Matrix *A, BinVector *x, BinVector *y, Byte init){
  Vector *s, *z;                     // The vectors Ay and A^Tx respectively
  unsigned int i,j;
  double obj, nobj;
  Row *nrow;

  s = initVector(A->m);
  z = initVector(A->n);
  
  // Initialization of pattern vector
  if (init != GGG && init!=NEIGHBOR)
    initialize(A, y, init);
  else{ // Initialization targets presence vector for these methods
    initialize(A, x, init);
    lMVMult(A, x, z);
    continuousSolveXS(y, z);
  }

  
  // Alternating Iterative Heuristic
  obj = -1.0;
  nobj = 0.0;
  for (j=0; j<MAXIT && nobj>obj; j++){
    obj  = nobj;
    rMVMult(A, y, s);
    continuousSolveXS(x, s);            // Solve for x
    lMVMult(A, x, z);
    nobj    = continuousSolveXS(y, z);  // Solve for y
  }

  // Record presence into matrix
  A->ones  = 0;
  for (nrow=A->rlist, i=0; nrow != ((Row *)(NULLROW)); nrow=nrow->next, i++){
    nrow->sel = element(x, i);
    A->ones  += nrow->sel;
  }
  A->zeros = A->m-A->ones;

  freeVector(z);
  freeVector(s);
}


int
hamming(Row *nrow, BinVector *y){
  unsigned int i;
  int h;
  Col *ncol;
  unsigned int c = 0;

  // Compute x^Ty
  ncol = nrow->clist;
  for (i=0; i< nrow->nz; i++)
    c += element(y, *(ncol+i));
  // Compute Hamming distance
  h = y->nz + nrow->nz - 2*c; 

  if (h<0) // This mustn't happen if everything's going OK!
    errexit("Negative hamming distance %d.\n",h);

  return h;
}

int
hammingRadius(Matrix *A, BinVector *y){
  Row *nrow;
  int h;
  int hr = 0;
  // Find maximum Hamming distance to pattern vector among all rows selected
  for (nrow = A->rlist; nrow != (Row *) NULLROW; nrow = nrow->next)
    if ( nrow->sel && (h=hamming(nrow,y)) > hr )
      hr = h;
  return hr;
}

int
partition(Matrix *A, BinVector *y, int epsilon, unsigned int mcs){
  unsigned int i,j;
  Row *nrow;
  Col *ncol;
  unsigned int *index;
  Vector *h;
  Bool stop;
  int hr, hd;
  unsigned int e;
  
  h = initVector(A->m);

  // Compute Hamming distances for each row
  h->max  = 0;
  for (i=0, nrow = A->rlist; nrow != (Row *) NULLROW; i++, nrow = nrow->next){
    nrow->sel = FALSE;
    h->list[i] = hamming(nrow,y);
    if (h->list[i] > h->max) h->max = h->list[i];
  }
  
  // Sort rows based on Hamming distance to pattern
  index = countingSort(h);

  // Progressing in sorted order, select rows closest to pattern vector
  A->ones = 0;
  stop = FALSE;
  hr = 0;

  for (i=0; i<A->m && !stop; i++){
    j    = index[i];
    nrow = A->index[j];
    hd   = h->list[j];
    if ( hd <= epsilon ){
      hr = hd;
      nrow->sel = TRUE;
      A->ones++;
    }
    else
      stop = TRUE;
  }

  if (!(A->ones)){
    memset(y->list, 0,(inWord(y->n)+1)*sizeof(unsigned int));
    nrow = A->index[index[0]];
    for (i=0; i<nrow->nz; i++){
      j = nrow->clist[i];
      y->list[inWord(j)] += setBit(bitOff(j));
    }
    y->nz = nrow->nz;
    for (nrow = A->rlist; nrow != (Row *) NULLROW; nrow = nrow->next){
      hd = hamming(nrow, y);
      if ( hd <= epsilon ){
	hr = hd;
	nrow->sel = TRUE;
	A->ones++;
      }
    }
  }

  tfree(index, h->n*sizeof(unsigned int));
  A->zeros = A->m - A->ones;
  
  freeVector(h);
  
  return hr;
}

void
recordLeaf(Matrix *A, BinVector *x, BinVector *y, Matrix *X, Matrix *Y, 
	   int hr, double *avg_hamming){
  unsigned int i,j;
  Row *crow, *nrow;
  Col *cptr;
  int found;
  
  if (X->m != Y->m)
    errexit("Presence and pattern matrices have unequal number of rows: X(%d), Y(%d)\n",
	    X->m,Y->m);
    
   // Check if found pattern is already in the list
  for (i=0, found=-1; i<Y->m && found<0; i++){
    int h;
    h = hamming(&(Y->rlist[i]), y);
    if ( !h )
      found = i;
  }
 
  // If so, record the presence vector in appropriate place
  if ( found >= 0 ){
    unsigned int movesize = 0;
    for (i=found+1; i<=X->m; i++){
      nrow = (X->rlist)+i;
      nrow->clist += x->nz;
      movesize    += nrow->nz;
    }
    crow = (X->rlist)+found;
    cptr = (crow->clist)+(crow->nz);
    memmove(cptr+(x->nz), cptr, movesize*sizeof(Col));
    for (i=0; i<x->n; i++){
      if (element(x,i)){
	(*cptr) = (A->index[i])->id;
	cptr++;
      }
    }
    crow->nz += x->nz;
    X->nz += x->nz;
    (*avg_hamming) += 1.0*hr*x->nz;
  }
  // Else, record presence and ... 
  else{
    crow = (X->rlist)+X->m;
    cptr = crow->clist;
    crow->nz = 0; 
    for (i=0; i<x->n; i++){
      if ( element(x,i) ){
	(*cptr) = (A->index[i])->id;
	cptr++;
	crow->nz++;
      }
    }
    nrow = crow+1;
    crow->next  = nrow;
    nrow->clist = cptr;
    nrow->nz    = 0;
    nrow->next  = (Row *) NULLROW;
    X->m++;
    X->nz += crow->nz;
    
    // ... pattern vectors
    crow = (Y->rlist)+Y->m;
    cptr = crow->clist;
    crow->nz = 0; 
    for (i=0; i<y->n; i++){
      if ( element(y,i) ){
	(*cptr) = i;
	cptr++;
	crow->nz++;
      }
    }
    nrow = crow+1;
    crow->next  = nrow;
    nrow->clist = cptr;
    nrow->nz    = 0;
    nrow->next  = (Row *) NULLROW;
    Y->m++;
    Y->nz += crow->nz;
    
    (*avg_hamming) += 1.0*hr*x->nz;
  }    
}

void
bndSolve(Matrix *A, Matrix *X, Matrix *Y, Byte init, int epsilon, unsigned int mcs, 
	 double *avg_hamming, Byte algorithm){
  Matrix *A0, *A1;
  Row *nrow;
  int hr;
  BinVector *x; 
  BinVector *y; 
  unsigned int i;

  if (A->m == 0)
    return;

  y = initBinVector(A->n);
  x = initBinVector(A->m);

  if (algorithm == DISCRETE)
    discreteRank1App(A, x, y, init);
  else
    continuousRank1App(A, x, y, init);

  if (!A->ones)
    errexit("Zero presence vector. A.m=%d A.n=%d A.nz=%d\n", A->m, A->n, A->nz);

  hr = hammingRadius(A, y);
  
  if ( !A->zeros ){ // and the hamming radius is small enough record pattern
    if (hr <= epsilon || A->m <= mcs ){
      recordLeaf(A, x, y, X, Y, hr, avg_hamming);
      freeBinVector(x);
      freeBinVector(y);
      return;
    }
    else{ // else partition based on hamming distances
      hr = partition(A, y, epsilon, mcs);
      if ( !A->zeros ){
	recordLeaf(A, x, y, X, Y, hr, avg_hamming);	
	freeBinVector(x);
	freeBinVector(y);
	return;
      }
      else{
	freeBinVector(x);
	A1 = (Matrix *) smalloc (sizeof(Matrix));
	A0 = (Matrix *) smalloc (sizeof(Matrix));
	splitMatrix(A, A0, A1);
	x = initBinVector(A1->m);
	x->nz = A1->m;
	memset(x->list, 0, (inWord(x->n)+1)*sizeof(unsigned int));
	for (i=0; i<x->n; i++){
	  x->list[inWord(i)] += setBit(bitOff(i));
	  A1->index[i]->sel = TRUE;
	}
	recordLeaf(A1, x, y, X, Y, hr, avg_hamming);	
	freeBinVector(x);
	freeBinVector(y);
	bndSolve(A0, X, Y, init, epsilon, mcs, avg_hamming, algorithm);
      }
    }
  }
  else if (hr <= epsilon){
    freeBinVector(x);
    A1 = (Matrix *) smalloc (sizeof(Matrix));
    A0 = (Matrix *) smalloc (sizeof(Matrix));
    splitMatrix(A, A0, A1);
    x = initBinVector(A1->m);
    x->nz = A1->m;
    memset(x->list, 0, (inWord(x->n)+1)*sizeof(unsigned int));
    for (i=0; i<x->n; i++){
      x->list[inWord(i)] += setBit(bitOff(i));
      A1->index[i]->sel = TRUE;
    }
    recordLeaf(A1, x, y, X, Y, hr, avg_hamming);	
    freeBinVector(x);
    freeBinVector(y);
    bndSolve(A0, X, Y, init, epsilon, mcs, avg_hamming, algorithm);
  }
  else{
    freeBinVector(x);
    freeBinVector(y);
    
    // split matrix based on presence of discovered pattern
    A1 = (Matrix *) smalloc (sizeof(Matrix));
    A0 = (Matrix *) smalloc (sizeof(Matrix));
    splitMatrix(A, A0, A1);
    
    // recursively try again for each submatrix
    bndSolve(A1, X, Y, init, epsilon, mcs, avg_hamming, algorithm);
    bndSolve(A0, X, Y, init, epsilon, mcs, avg_hamming, algorithm);
  }
}

unsigned int
bnd(Matrix *A, Matrix *X, Matrix *Y, Byte init, int epsilon, unsigned int mcs, 
    double *avghamming, Byte algorithm){
  
   // Initialize presence and pattern matrices
  X->m     = 0;
  X->n     = A->m;
  X->nz    = 0;
  X->rlist = (Row *) smalloc((1+A->m)*sizeof(Row));
  X->rlist[0].clist = (Col *) smalloc(A->m*sizeof(Col));
  
  Y->m     = 0;
  Y->n     = A->n;
  Y->nz    = 0;
  Y->rlist = (Row *) smalloc((1+A->m)*sizeof(Row));
  Y->rlist[0].clist = (Col *) smalloc((2*A->nz)*sizeof(Col));
  
  bndSolve(A, X, Y, init, epsilon, mcs, avghamming, algorithm);
  
  if (X->m != Y->m)
    errexit("Presence and pattern matrices have unequal number of rows: X(%d), Y(%d)\n",
	    X->m,Y->m);
  
  X->rlist[X->m-1].next = (Row *) NULLROW;
  Y->rlist[Y->m-1].next = (Row *) NULLROW;

  return X->m;
}
