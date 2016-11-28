/********************************************************************************
* Authors     : Mehmet Koyuturk & Ananth Grama                                 *
* Institution : Purdue University, Department of Computer Sciences              *
* Last update : 03/31/2004                                                      *
* Copyright   : See README file for copyright information.                      *
********************************************************************************/

#include "../include/bnd.h"

int
main(int argc, char *argv[]){
  Matrix *A;                         // Original matrix
  Matrix *X, *Y;                     // Presence & Pattern matrices 
  unsigned int nleaves;              // # of identified patterns (leaves of recursion tree)
  ArgList *arguments;                // List of input parameters
  double timestart, timeend;         // Record time
  double avghamming;                 // Average hamming radius of resulting decumposition

  arguments = processArguments(argc, argv);
  initMemoryPool(pow2(arguments->poolsize),pow2(arguments->tmpsize));
  A = readMatrix(arguments->fname);
  displayMatrix(A, arguments->fname);

  avghamming = 0.0;
  X = (Matrix *) smalloc (sizeof(Matrix));
  Y = (Matrix *) smalloc (sizeof(Matrix));
  timestart=((double)(clock()))*1.0e-06;
  nleaves = bnd(A, X, Y, arguments->init, arguments->epsilon,arguments->mcs, 
		&avghamming, arguments->algorithm);
  timeend=((double)(clock()))*1.0e-06;
  
  avghamming = avghamming/A->m; 
  
  mprintf("Decomposition with (threshold=%d, init=",arguments->epsilon);
  switch (arguments->init){
  case ALLONES    : mprintf("Allones)\n"); break;
  case CENTER     : mprintf("Center)\n"); break;
  case INITMAX    : mprintf("Maximum)\n"); break;
  case PARTITION  : mprintf("Partition)\n"); break;
  case GGG        : mprintf("GGG)\n"); break;
  case NEIGHBOR   : mprintf("Neighbor)\n"); break;
  case RANDOMROW  : mprintf("Random-row)\n"); break;
  case RANDOM     : mprintf("Random)\n"); break;
  default : mprintf("None?)\n");
  }
  mprintf("-------------------------------------------------------\n");
  mprintf(" Number of patterns                = %d\n",nleaves);
  mprintf(" Average hamming radius            = %2.2lf\n",avghamming);
  mprintf(" Run time                          = %2.2lf secs.\n",timeend-timestart); 
  mprintf("*******************************************************\n");

  if (arguments->write)
    writeOutput(arguments->fname, X, Y);

  return 1;
}
