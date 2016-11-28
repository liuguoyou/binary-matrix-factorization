/********************************************************************************
* Authors     : Mehmet Koyuturk, Ananth Grama                                   *
* Institution : Purdue University, Department of Computer Sciences              *
* Last update : 03/30/2004                                                      *
* Copyright   : See README file for copyright information and terms of use.     *
********************************************************************************/

#ifndef _SYSTEM_H
#define _SYSTEM_H

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#define MAXLINE 1024*1024  // maximum length of line that can be read from a file

// boolean constants
#define TRUE       1
#define FALSE      0

#define NONE      -1

// Memory alllocation constants
#define POOLSIZE  28
#define TMPSIZE   24

typedef unsigned char Bool;
typedef unsigned char Byte;

//exported functions 
void errexit(const char *format, ...);
void mprintf(const char *format, ...);
void initMemoryPool(int ps, int ts);
void *smalloc(int size);
void *tmalloc(int size);
void tfree(void *ptr, int size);
void sfree(void *ptr, int size);
void printArray(int *list, int n, char *name);

#endif

