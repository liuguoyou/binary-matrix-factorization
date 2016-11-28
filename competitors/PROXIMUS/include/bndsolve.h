/********************************************************************************
* Authors     : Mehmet Koyuturk & Ananth Grama                                  *
* Institution : Purdue University, Department of Computer Sciences              *
* Last update : 03/31/2004                                                      *
* Copyright   : See README file for copyright information.                      *
********************************************************************************/

#ifndef _BNDSOLVE_H
#define _BNDSOLVE_H

#include <math.h>
#include "system.h"
#include "binvector.h"
#include "vector.h"
#include "matrix.h"
#include "initialize.h"

#define MAXIT      100

// Finds rank-one approximation xy^T for binary matrix A using initialization strategy init
void discreteRank1App(Matrix *A, BinVector *x, BinVector *y, Byte init);
void continuousRank1App(Matrix *A, BinVector *x, BinVector *y, Byte init);
// Computes normalized Hamming distance between two binary vectors (each row is a binary vector)
int hamming(Row *nrow, BinVector *y);
// Computes normalized Hamming radius of selected rows of A centered around binary vector y
int hammingRadius(Matrix *A, BinVector *y);
// Partitions (rows of) matrix A based on normalized Hamming distances of rows to pattern vector y
int partition(Matrix *A, BinVector *y, int epsilon, unsigned int mcs);
// Records discovered pattern with its corresponding presence vector
void recordLeaf(Matrix *A, BinVector *x, BinVector *y, Matrix *X, Matrix *Y, int hr, double *avg_hamming);
// Recursively finds rank-one approximations to binary matrix A
void bndSolve(Matrix *A, Matrix *X, Matrix *Y, Byte init, int epsilon, unsigned int mcs, double *avg_hamming, Byte algorithm);
// Computes non-orthogonal approximation XY^T for binary matrix A
unsigned int bnd(Matrix *A, Matrix *X, Matrix *Y, Byte init, int epsilon, unsigned int mcs, double *avg_hamming, Byte algorithm);

#endif
