#include "../include/system.h"
#include "../include/binvector.h"
#include "../include/matrix.h"

Matrix *matrixMultiply(Matrix *X, Matrix *Y){
  Row *brow, *xrow, *yrow;
  Col *cptr;
  unsigned int *rnz;
  unsigned int i, j, nz;
  Matrix *B;
  
  rnz = (unsigned int *) smalloc(X->n * sizeof(unsigned int));
  memset(rnz, 0, X->n*sizeof(unsigned int) );
  xrow = X->rlist;
  yrow = Y->rlist;
  nz = 0;
  for (i = 0; i<X->m; i++){
    for (j = 0; j<xrow->nz; j++){
      if (rnz[xrow->clist[j]])
	errexit("Row %d duplicated at vector %d\n", xrow->clist[j], i); 
      rnz[xrow->clist[j]] = yrow->nz;
      nz += yrow->nz;
      //printf("%d(%d) => %d\n",i,j,nz); fflush(stdout);
    }
    xrow = xrow->next;
    yrow = yrow->next;
  }
  //printf("Total nonzeros in B=%d\n",nz); fflush(stdout);
  B = initMatrix(X->n, Y->n, nz);
  //printf("Matrix B initialized %d\n",B->rlist[0].clist); fflush(stdout);
  B->rlist[0].nz = rnz[0];
  for (i=1; i<B->m; i++){
    brow = &(B->rlist[i]);
    B->rlist[i-1].next = brow;
    brow->clist = (B->rlist[i-1]).clist + rnz[i-1];
    brow->nz = rnz[i];
    //printf("%d => %d [%d]\n",i,B->rlist[i].clist,B->rlist[i].nz);
  }
  brow->next = (Row *) NULLROW;
  //printf("Matrix rows done..\n"); fflush(stdout);
  xrow = X->rlist;
  yrow = Y->rlist;
  for (i = 0; i<X->m; i++){
    unsigned int size= yrow->nz * sizeof(unsigned int);
    for (j = 0; j<xrow->nz; j++){
      //printf("%d(%d)=>%d [%d->%d] to %d\n", i,j,xrow->clist[j],yrow->nz,size,(B->rlist[xrow->clist[j]]).clist); fflush(stdout);
      memcpy((B->rlist[xrow->clist[j]]).clist, yrow->clist, size);
    }
    xrow = xrow->next;
    yrow = yrow->next;
  }
  return B;
}

unsigned int
rowMatch(Row *a, Row *b){
  unsigned int i,j;
  unsigned int match = 0;
  
  for (i=0; i<a->nz; i++){
    unsigned int key = a->clist[i];
    Bool stop = FALSE;
    for (j=0; j<b->nz && !stop; j++){
      if (b->clist[j] == key){
	match++;
	stop = TRUE;
      }
    }
  }
  
  return match;
}

void 
compare(Matrix *A, Matrix *B, double *precision, double *recall, unsigned int *error, double *hamming, double *density){
  Row *arow, *brow;
  unsigned int i;
  unsigned int match = 0;
  int h = 0;
  double d = 0;
  
  arow = A->rlist;
  brow = B->rlist;
  for (i=0; i<A->m; i++){
    unsigned int mt = rowMatch(arow, brow);
    match += mt;
    h     += arow->nz+brow->nz-2*mt;
    d     += (1.0*mt)/brow->nz;
    arow = arow->next;
    brow = brow->next;
  }
  (*precision) = (100.0*match)/B->nz;
  (*recall)    = (100.0*match)/A->nz;
  (*error)     = A->nz+B->nz-2*match;
  (*hamming)   = (1.0*h)/A->m;
  (*density)   = d/A->m;
}


int
main(int argc, char *argv[]){
  Matrix *A, *X, *Y, *B;
  char fname[100];
  double precision, recall, hamming, density;
  unsigned int error;

  if ( argc <2 )
    errexit("Run as\n%s <fname>\n",argv[0]);
  
  initMemoryPool(pow2(POOLSIZE),pow2(TMPSIZE));

  strcpy(fname, argv[1]);
  A = readMatrix(fname);
  strcat(fname, ".X.out");
  X = readMatrix(fname);
  strcpy(fname, argv[1]);
  strcat(fname, ".Y.out");
  Y = readMatrix(fname);
  
  B = matrixMultiply(X, Y);
  
  compare(A, B, &precision, &recall, &error, &hamming, &density);

  printf("Original      Matrix(%d X %d): %d nonzeros\n",A->m, A->n, A->nz); 
  printf("Approximation Matrix(rank %d): %d nonzeros\n",X->m, B->nz);
  printf("Average Hamming Distance = %6.2lf\n",hamming);
  printf("Average Density          = %6.2lf\n",density);
  printf("Error                    = %6d\n", error);
  printf("Precision                = %6.2lf\n", precision);
  printf("Recall                   = %6.2lf\n", recall);
  
  return(1);
}
