/********************************************************************************
* Authors     : Mehmet Koyuturk & Ananth Grama                                  *
* Institution : Purdue University, Department of Computer Sciences              *
* Last update : 03/31/2004                                                      *
* Copyright   : See README file for copyright information.                      *
********************************************************************************/

#include "../include/initialize.h"

void
initialize(Matrix *A, BinVector *y, Byte method){

  srand( (unsigned int) time(NULL));
  memset(y->list,0,(inWord(y->n)+1)*sizeof(unsigned int));
  y->nz = 0;

  switch (method){
  case ALLONES    : initAllOnes(y);      break; 
  case CENTER     : initCenter(A, y);    break;
  case INITMAX    : initMaximum(A,y);    break;
  case PARTITION  : initPartition(A,y);  break;
  case GGG        : initGGG(A,y);        break;
  case NEIGHBOR   : initNeighbor(A,y);   break;
  case RANDOMROW  : initRandomRow(A,y);        break;
  case RANDOM     : initRandom(A,y);     break;
  default         : errexit("Invalid initialization method %d.\n",method);
  }

  if (!(y->nz))
    errexit("0-vector initialization.\n");
}

void
initAllOnes(BinVector *y){
  unsigned int i;
  for (i=0; i<=inWord(y->n); i++)
    y->list[i] = ONES;
  y->nz = y->n;
}

void
initCenter(Matrix *A, BinVector *y){
  unsigned int i, j;
  Row *nrow;
  double avg;
  unsigned int *count;

  count = (unsigned int *) tmalloc(y->n * sizeof(unsigned int));
  memset(count, 0, y->n*sizeof(unsigned int));

  for(nrow=A->rlist; nrow != ((Row *)NULLROW); nrow=nrow->next){
    for (j=0; j< nrow->nz; j++)
      count[nrow->clist[j]]++;
  }

  avg = (1.0*A->nz)/y->n;
  y->nz = 0;
  for (i=0; i<y->n; i++){
    if (count[i] > avg){
      y->list[inWord(i)] += setBit(bitOff(i));
      y->nz++;
    }
  }

  tfree(count, y->n*sizeof(unsigned int));
}

void
initMaximum(Matrix *A, BinVector *y){
  unsigned int i, j;
  Row *nrow;
  double avg;
  unsigned int *count;
  unsigned int max;

  count = (unsigned int *) tmalloc(y->n * sizeof(unsigned int));
  memset(count, 0, y->n*sizeof(unsigned int));

  for(nrow=A->rlist; nrow != ((Row *)NULLROW); nrow=nrow->next){
    for (j=0; j< nrow->nz; j++)
      count[nrow->clist[j]]++;
  }
  
  max = 0;
  for(i=0; i<y->n;i++){
    if (count[i]>count[max])
      max = i;
  }

  y->list[inWord(max)] += setBit(bitOff(max));
  y->nz = 1;

  tfree(count, y->n*sizeof(unsigned int));
}

void
initPartition(Matrix *A, BinVector *y){
  unsigned int i, j, c;
  Row *nrow;
  double avg;
  unsigned int *count;
  unsigned int sep,sum;
  int dif;
  Bool flag;
  int t = (A->m/2)+1;

  count = (unsigned int *) tmalloc(y->n*sizeof(unsigned int));
  memset(count, 0, y->n*sizeof(int));

  for(nrow=A->rlist; nrow != ((Row *)NULLROW); nrow=nrow->next){
    for (j=0; j< nrow->nz; j++)
      count[nrow->clist[j]]++;
  }

  dif = A->m+1;
  for(i=0; i<y->n;i++){
    int d = abs(count[i]-t);
    if ( d < dif){
      dif = d;
      sep = i;
    }
  }

  memset(count, 0, y->n*sizeof(int));
  for(nrow=A->rlist; nrow != ((Row *)NULLROW); nrow=nrow->next){
    flag = FALSE;
    for (j=0; j< nrow->nz; j++){      
      if ( nrow->clist[j] == sep) 
	flag = TRUE;
    }
    if (flag){
      for (j=0; j< nrow->nz; j++)      
	count[nrow->clist[j]]++;
    }
  }

  sum = 0;
  c = 0;
  for(i=0; i<y->n;i++){
    sum += count[i];
    if (count[i] > 0)
      c++;
  }

  avg = (1.0*sum)/c;
  for(i=0; i<y->n;i++){
    if (count[i] >= avg){  
      y->list[inWord(i)] += setBit(bitOff(i));
      y->nz++;
    }
  }
  
  tfree(count, y->n*sizeof(unsigned int));
}

void
initGGG(Matrix *A, BinVector *y){
  unsigned int i, j, r;
  int set;
  Row *nrow;
  unsigned int *setscore[2], rowscore[2], totscore[2], setcount[2];
  double avg;
  double rating[2];

  setscore[0] = (unsigned int *)tmalloc(A->n*sizeof(unsigned int));
  setscore[1] = (unsigned int *)tmalloc(A->n*sizeof(unsigned int));
  memset(setscore[0], 0, A->n*sizeof(unsigned int));
  memset(setscore[1], 0, A->n*sizeof(unsigned int));
  totscore[0] = totscore[1] = 0;
  setcount[0] = setcount[1] = 0;

  for (nrow = A->rlist, r=0; nrow != (Row *) NULLROW; nrow = nrow->next,r++){
    rowscore[0] = rowscore[1] = 0;
    for (i=0; i<nrow->nz; i++){
      j = *(nrow->clist+i);
      rowscore[0] += setscore[0][j];
      rowscore[1] += setscore[1][j];
    }
    if ( r >= 10){
      rating[0] = (1.0*rowscore[0]/totscore[0]);
      rating[1] = (1.0*rowscore[1]/totscore[1]);
    }
    if (r<10)
      set = r % 2;
    else if ( rating[1] > rating[0] )
      set = 1;
    else if ( rating[1] < rating[0] )
      set = 0;
    else
      set = rand() % 2;
    for (i=0; i<nrow->nz; i++){
      j = *(nrow->clist+i);
      (setscore[set][j])++;
      totscore[set]++;
    }
    setcount[set]++;
    nrow->sel = set;
  }
  if (setcount[0] > setcount[1] ) set = 0;
  else set = 1;
  if (setcount[set] == 0) set = 1-set; 
  avg = (1.0*totscore[set])/y->n;

  for (nrow = A->rlist, r=0; nrow != (Row *) NULLROW; nrow = nrow->next,r++){
    if ( nrow->sel == set ){
      y->list[inWord(r)] += setBit(bitOff(r));
      y->nz++;
    }
    nrow->sel = 0;
  }

  tfree(setscore[1], A->n*sizeof(unsigned int));
  tfree(setscore[0], A->n*sizeof(unsigned int));
}

void
initNeighbor(Matrix *A, BinVector *y){
  unsigned int i,j,r;
  Row *nrow, *srow;
  unsigned int *mark;
  int size = 0;

  i = rand() % A->m;
  srow = A->index[i];
  mark = (unsigned int *)tmalloc(A->n*sizeof(unsigned int));
  memset(mark, 0, A->n*sizeof(unsigned int));
  for (i=0; i<srow->nz; i++){
      j = *(srow->clist+i);
      mark[j] = 1;
  }

  for ( r=0, nrow = A->rlist; nrow != (Row *) NULLROW; r++, nrow = nrow->next ){
    nrow->sel = FALSE;  
    for (i=0; i<nrow->nz; i++){
	j = *(nrow->clist+i);
	if (mark[j]) nrow->sel = TRUE; 
    }
  
    if ( nrow->sel ){
      y->list[inWord(r)] += setBit(bitOff(r));
      y->nz++;
    }
    nrow->sel = FALSE;
  }
  tfree(mark, A->n*sizeof(unsigned int));
}

void
initRandomRow(Matrix *A, BinVector *y){
  unsigned int i,j,r;
  Row *srow;

  i = rand() % A->m;
  srow = A->index[i];
  for (i=0; i<srow->nz; i++){
      j = *(srow->clist+i);
      y->list[inWord(j)] += setBit(bitOff(j));
      y->nz++;
  }
}

unsigned int 
l2(unsigned int n){
  unsigned int m;
  unsigned int k = 0;
  for (m=n; m>0; m=m>>1)
    k++;
  return k;
}
  
void
initRandom(Matrix *A, BinVector *y){
  unsigned int i, j, c, r, k;
  Row *nrow;
  unsigned int *count;

  count = (unsigned int *) tmalloc(A->n*sizeof(int));
  memset(count, 0, A->n*sizeof(unsigned int));
 
  for(nrow=A->rlist; nrow != ((Row *)NULLROW); nrow=nrow->next){
    for (j=0; j<nrow->nz; j++)
      count[nrow->clist[j]]++;
  }
 
  c = 0;
  for (i=0; i<y->n;i++)
    if (count[i])
      count[c++] = i;
 
  r = A->nz/A->m;
  if (r>c) r=c;
  for (i=0;i<r;i++){
    j = rand() % c;
    k = count[j];
    y->list[inWord(k)] += setBit(bitOff(k));
    y->nz++;
    count[j] = count[--c];
  }
  tfree(count, A->n*sizeof(int));
}
