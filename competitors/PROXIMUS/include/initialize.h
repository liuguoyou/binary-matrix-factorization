/********************************************************************************
* Author      : Mehmet Koyuturk                                                 *
* Supervisor  : Ananth Grama                                                    *
* Institution : Purdue University, Department of Computer Sciences              *
* Last update : 05/02/2003                                                      *
* Copyright   : See README file for copyright information.                      *
********************************************************************************/

#ifndef _INITIALIZE_H
#define _INITIALIZE_H

#include <math.h>
#include <time.h>
#include "system.h"
#include "binvector.h"
#include "matrix.h"

#define ALLONES     1                  // initialization schemes
#define CENTER      2
#define INITMAX     3
#define PARTITION   4
#define GGG         5
#define NEIGHBOR    6
#define RANDOMROW   7
#define RANDOM      8

// Resets vector, calls appropriate initialization routine
void initialize(Matrix *A, BinVector *y, Byte method);
// Initializes the pattern vector to vector of all ones
void initAllOnes(BinVector *y);
// Initializes the pattern vector to the center of points(rows of A)
// in the n-dimensional hypercube (n is number of columns of A)
void initCenter(Matrix *A, BinVector *y);
// Initializes the pattern vector to a vector with a one
// in the entry that corresponds to the column of A that has
// maximum number of nonzeros
void initMaximum(Matrix *A, BinVector *y);
// Partitions the rows of A into two parts along the dimension
// that results in the most balanced partition, initializes 
// the pattern vector  to the center of the smaller part
void initPartition(Matrix *A, BinVector *y);
// Implements an adapted version of Greedy Graph Growing
// to get a balanced partition the rows of A, initializes the pattern
// vector to the center of one part
void initGGG(Matrix *A, BinVector *y);
// Initializes the pattern vector to the center of the row set
// that includes rows that share a nonzero with a randomly selected
// row
void initNeighbor(Matrix *A, BinVector *y);
// Initializes the pattern vector to a randomly selected row of A
void initRandomRow(Matrix *A, BinVector *y);
// Initializes the pattern vector by setting randomly
// selected entries
void initRandom(Matrix *A, BinVector *y);

//log2
unsigned int l2(unsigned int n);

#endif
