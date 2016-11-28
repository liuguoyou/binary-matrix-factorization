/********************************************************************************
* Author      : Mehmet Koyuturk                                                 *
* Supervisor  : Ananth Grama                                                    *
* Institution : Purdue University, Department of Computer Sciences              *
* Last update : 04/31/2004                                                      *
* Copyright   : See README file for copyright information.                      *
********************************************************************************/

#include "../include/matrix.h"

Matrix *
initMatrix(unsigned int m, unsigned int n, unsigned int nz){
  Matrix *A = (Matrix *) smalloc(sizeof(Matrix));
  A->m      = m;   
  A->n      = n;
  A->nz     = nz;
  A->rlist  = (Row *) smalloc(m*sizeof(Row));
  A->index  = (Row **) smalloc(m*sizeof(Row *));
  A->rlist->clist = (Col *) smalloc(nz*sizeof(Col));
  return A;
}

Matrix *
readMatrix(char *fname){
  Matrix *A;
  FILE *fp;
  unsigned int m, n, nz;
  unsigned int i, k, cl;
  char *line, *oldstr, *newstr;
  Row *next, *nrow;
  Col *cptr;
  
  fp = fopen(fname, "r");

  if ( !((int) fp) )
    errexit("Can't open file %s for reading.\n",fname); 

  line = (char *) tmalloc(MAXLINE);

  do{
    fgets(line, MAXLINE, fp);
  } while (line[0] == '%');

  k = sscanf(line,"%d %d %d", &m, &n, &nz);
  if (k < 3)
    errexit("Can't read matrix parameters from %s. Top line has only %d(<3) entries.\n", fname,k);
  else if ( k > 3)
    errexit("Too many matrix parameters in file %s: %d(>3).\n",fname,k);

  A = initMatrix(m, n, nz); 
  
  nrow = A->rlist;
  cptr = nrow->clist;

  for (i=0; i<m; i++){                // read non-zeros of row-by-ros (so line-by-line) 
    A->index[i] = nrow;
    nrow->id = i;
    nrow->nz = 0;
    do{
      fgets(line, MAXLINE, fp);
    } while (line[0] == '%');
    oldstr = line;

    if (strlen(line) >= MAXLINE-5)
      errexit("Too many non-zeros in row %d. Maximum number of characters allowed in a line is %d.\n",
	      i,MAXLINE-5);
    nrow->clist = cptr;

    for (;;){
      cl = (int)strtol(oldstr, &newstr, 10);
      if (oldstr==newstr)
        break;
      oldstr = newstr;

      if ( ( cl < 0 ) || ( cl >= n ) )
        errexit("Invalid column id %d in list of row %d\n", cl, i);

      *cptr = cl;
      cptr++;
      nrow->nz++;
    }
    nrow->sel = 0;
    if (nrow->nz == 0) 
      errexit("Empty row %d on line %d of file %s.\n",i,i+1,fname);

    if ( (nrow-A->rlist) < (A->m-1) ){
      next = nrow;
      next++;
      nrow->next = next;
      nrow = next;
    }
  }
  nrow->next = (Row *) NULLROW;

  if ((cptr-(A->rlist->clist)) != A->nz)
    errexit("Number of nonzeros read: %d, Number of nonzeros declared: %d\n",
	    cptr-A->rlist->clist,A->nz);

  tfree(line, MAXLINE);
  return A;
}

void 
writeMatrix(Matrix *A, char *fname){
  unsigned int i,c;
  Row *nrow;
  FILE *ofp = fopen(fname, "w");
  if (ofp == NULL){
    mprintf("***WARNING: Can't open file %s for writing.\n",fname);
    return;
  }
  fprintf(ofp, "%d %d %d\n", A->m, A->n, A->nz);
  c = 0;
  for(nrow = A->rlist; nrow != ((Row *)NULLROW); nrow=nrow->next){
    for (i=0; i<nrow->nz; i++){
      fprintf(ofp,"%d ", nrow->clist[i]);
    }
    fprintf(ofp, "\n");
  }
  fclose(ofp);
}

void 
displayMatrix(Matrix *A, char *fname){
  mprintf("*******************************************************\n");
  mprintf("Input matrix (%s):\n", fname); 
  mprintf("-------------------------------------------------------\n");
  mprintf(" Number of rows     = %10d\n", A->m);
  mprintf(" Number of columns  = %10d\n", A->n);
  mprintf(" Number of nonzeros = %10d\n", A->nz);
  mprintf(" Density            = %10.4lf\n",(1.0*A->nz)/(A->m*A->n)); 
  mprintf("*******************************************************\n");
}

void
rMVMult(Matrix *A, BinVector *x, Vector *s){
  unsigned int i, j, jptr, jeptr;
  unsigned int wid, wof;
  unsigned int Aw, xw, sw;
  Col c;
  int cw;
  Row *nrow;

  if ( s->n != A->m || x->n != A->n)
    errexit("Right matrix-vector multiplication: sizes don't match!\nA=%dx%d s=%dx1 x=%dx1\n",
	    A->m,A->n,s->n,x->n);

  memset(s->list, 0, s->n*sizeof(unsigned int)); // Initialize s to zero vector
  s->max = 0;
  nrow = A->rlist;

  for (i=0; i<A->m; i++){                        // For each row of A
    for (j=0; j<nrow->nz; j++){                 
      c = nrow->clist[j];                        // for each nonzero of row
      s->list[i] += element(x,c);              // increment corresponding entry of s 
    }
    nrow=nrow->next;
    if (s->list[i] > s->max) s->max = s->list[i];
  }
   
}


void
lMVMult(Matrix *A, BinVector *x, Vector *s){
  unsigned int i,j;
  Row *nrow;

  if ( s->n != A->n || x->n != A->m)
    errexit("Left matrix-vector multiplication: sizes don't match!\nA=%dx%d s=%dx1 x=%dx1\n",
	    A->m,A->n,s->n,x->n);


  memset(s->list, 0, s->n*sizeof(unsigned int)); // initialize s to zero vector
  
  nrow = A->rlist;
  for (i=0; i<A->m; i++){                        // for each row of A
    if ( element(x,i) ){                         // if corresponding entry of x is nonzero
      for (j=0; j<nrow->nz; j++)                 // for each nonzero of row
        s->list[nrow->clist[j]]++;               // increment corresponding entry of s 
    }
    nrow = nrow->next;
  }

  s->max = 0;                                    // Maximum of s will be used in counting sort
  for (i=0; i<s->n; i++)
    if (s->list[i] > s->max) s->max = s->list[i];
}

void 
splitMatrix(Matrix *A, Matrix *A0, Matrix *A1){
  int i;
  Row *nrow;
  Row **prev0, **prev1;

  if (!A->ones || !A->zeros){
    errexit("Matrix with no rows[%d %d].\n", A->zeros, A->ones);
  }

  // Initialize child matrices
  A0->m  = A1->m  = 0;
  A0->n  = A1->n  = A->n;
  A0->nz = A1->nz = 0;
   
  prev0 = &(A0->rlist);
  prev1 = &(A1->rlist);
  
  // for each row of x, identify corresponding child matrix, put that row in that matrix
  for (nrow=A->rlist; nrow != ((Row *)(NULLROW)); nrow=nrow->next){
    if (nrow->sel){
      (*prev1)     = nrow;
      A1->m++;
      A1->nz      += nrow->nz;
      nrow->sel    = FALSE;
      prev1        = &(nrow->next);
    }
    else{
      (*prev0)     = nrow;
      A0->m++;
      A0->nz      += nrow->nz;
      nrow->sel    = FALSE;
      prev0        = &(nrow->next);
    }
  }
 
  (*prev0) = (Row *) NULLROW;
  (*prev1) = (Row *) NULLROW;

  // Create index of rows for child matrices
  A0->index = A->index;
  A1->index = A->index+A0->m;
  for (i=0, nrow=A0->rlist; nrow != ((Row *)(NULLROW)); nrow=nrow->next, i++)
    A0->index[i] = nrow;
  for (i=0, nrow=A1->rlist; nrow != ((Row *)(NULLROW)); nrow=nrow->next, i++)
    A1->index[i] = nrow;
}
