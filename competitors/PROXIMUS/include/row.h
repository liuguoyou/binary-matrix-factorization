/********************************************************************************
* Author      : Mehmet Koyuturk                                                 *
* Supervisor  : Ananth Grama                                                    *
* Institution : Purdue University, Department of Computer Sciences              *
* Last update : 05/02/2003                                                      *
* Copyright   : See README file for copyright information.                      *
********************************************************************************/
#ifndef _ROW_H
#define _ROW_H

#include "system.h"

#define NULLROW    -1

typedef unsigned int Col;

struct Row{
  Col *clist;                         // pointer to list of nozeros
  struct Row *next;                   // pointer to next row
  unsigned int nz;                    // number of nonzeros in this row
  unsigned char sel;                  // is the pattern present in this row?
  unsigned int id;                    // id of row in the original matrix
};
typedef struct Row Row;               // row of a matrix


#endif
