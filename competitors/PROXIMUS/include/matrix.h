/********************************************************************************
* Authors     : Mehmet Koyuturk & Ananth Grama                                  *
* Institution : Purdue University, Department of Computer Sciences              *
* Last update : 04/31/2004                                                      *
* Copyright   : See README file for copyright information.                      *
********************************************************************************/

#ifndef _MATRIX_H
#define _MATRIX_H

#include "system.h"
#include "row.h"
#include "binvector.h"
#include "vector.h"

typedef struct{
  unsigned int m;                     // number of rows
  unsigned int n;                     // number of columns
  unsigned int nz;                    // number of nonzeros
  Row *rlist;                         // pointer to first row
  int ones;                           // # of rows containing discovered pattern
  int zeros;                          // not containing
  Row **index;                        // pointer to index of rows
} Matrix;       
// Initializes a matrix
Matrix *initMatrix(unsigned int m, unsigned int n, unsigned int nz);
// Reads matrix from file fname
Matrix *readMatrix(char *fname);
// Writes matrix A into file fname
void writeMatrix(Matrix *A, char *fname);
// Displays parameters of matrix A
void displayMatrix(Matrix *A, char *fname);
// Computes s=Ax
void rMVMult(Matrix *A, BinVector *x, Vector *s);
// Computes s=A^Tx
void lMVMult(Matrix *A, BinVector *x, Vector *s);
// Splits matrix A into two submatrices A0 and A1 based on rows
void splitMatrix(Matrix *A, Matrix *A0, Matrix *A1);

#endif
