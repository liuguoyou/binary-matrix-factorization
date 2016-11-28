/********************************************************************************
* Authors     : Mehmet Koyuturk & Ananth Grama                                  *
* Institution : Purdue University, Department of Computer Sciences              *
* Last update : 03/31/2004                                                      *
* Copyright   : See README file for copyright information.                      *
********************************************************************************/
#ifndef _IO_H
#define _IO_H

#include "system.h"
#include "initialize.h"
#include "bnd.h"

typedef struct{
  char fname[100];   // input file name containing input matrix
  Byte algorithm;    // algorithm for rank-one approximation
  Byte init;         // initialization strategy
  int epsilon;       // bound on hamming radius
  Bool write;        // flag for writing output
  unsigned int mcs;  // minimum cluster size
  int poolsize;      // log of static memory pool size
  int tmpsize;       // log of temporary memory pool size
} ArgList;           // argument list

ArgList *processArguments(int argc, char *argv[]);
void writeOutput(char *fname, Matrix *X, Matrix *Y);

#endif
