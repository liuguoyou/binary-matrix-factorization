/*****************************************************************************
* Authors     : Mehmet Koyuturk & Ananth Grama                               *
* Institution : Purdue University, Department of Computer Sciences           *
* Last update : 05/02/2003                                                   *
* Copyright   : See README file for copyright information and terms of use.  *
******************************************************************************/

#ifndef _BINVECTOR_H
#define _BINVECTOR_H

#include "system.h"

#define OFSMASK    0x0000001f          // 31
#define ZERO       0x00000000          // 0
#define ONE        0x00000001          // 1
#define ONES       0xffffffff          // !0

#define bitOff(x) (x&OFSMASK)         // returns the bit offset of given index in a binary vector
#define inWord(x) (x>>5)              // returns the word id of a given index in a binary vector
#define xor(x, y) ( (x & ~y) | (~x & y) ) 
#define indexOf(i,j,m) (i*m+j)        // maps two dimensional index to one dimension

typedef struct{
  unsigned int n;                     // vector size
  unsigned int nz;                    // number of nonzeros in the vector
  unsigned int *list;                 // bit array of elements
} BinVector;                          // binary vector

// Initializes binary vector
BinVector *initBinVector(unsigned int n);
// Deallocates space for binary vector
void freeBinVector(BinVector *x);
// Returns i'th element of binary vector x
unsigned int element(BinVector *x, unsigned int i);
// Returns an unsigned integer with only i'th bit equal to one
unsigned int setBit(unsigned int i);
// Counts the number of ones in a word
unsigned int noOfOnes(unsigned int z);
// Displays binary vector for debugging
void displayBinVector(BinVector *x, char *name);

#endif
