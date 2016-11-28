/********************************************************************************
* Authors     : Mehmet Koyuturk & Ananth Grama                                  *
* Institution : Purdue University, Department of Computer Sciences              *
* Last update : 04/06/2004                                                      *
* Copyright   : See README file for copyright information.                      *
********************************************************************************/

#ifndef _VECTOR_H
#define _VECTOR_H

#include "system.h"
#include "binvector.h"

#define DISCRETE                 1
#define CONTINUOUS               2

typedef struct{
  unsigned int n;                     // vector size
  unsigned int *list;                 // array of elements
  unsigned int max;                   // upper bound on elements, used for counting sort
} Vector;                             // integer vector

// Initializes vector
Vector *initVector(unsigned int n);
// Deallocates space for vector
void freeVector(Vector *s);
// Returns an index that contains a sorted ordering of elements of s
unsigned int *countingSort(Vector *s);
// Computes x to maximize 2*x^Ts-||x||||y|| heuristicly
int discreteSolveXS(BinVector *x, Vector *s, unsigned int y);
// Computes x to maximize x^Ts/||x|| heuristicly
double continuousSolveXS(BinVector *x, Vector *s);
// Displays vector for debugging
void displayVector(Vector *s);

#endif
