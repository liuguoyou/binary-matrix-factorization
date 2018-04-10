#ifndef UTIL_H
#define UTIL_H
#include "binmat.h"

typedef std::pair<idx_t,idx_t> aux_t;


void inc_verbosity();

char get_verbosity();

/**
 * Saves an image where each column of D is transformed into a square and shown
 * as a tile.
 * @param D matrix to be rendered as a mosaic of square tiles
 * @param fname file name of image to write to
 */
void render_mosaic(const  binary_matrix& D, const char* fname);

/** 
 * utility for sorting according to counts 
 */
void counting_sort(aux_t* s, idx_t n);

/** 
 * write a binary matrix A as a PBM image of file name fname 
 * @param A matrix to be saved as image
 * @param fname file name of image to write to
 */
int write_pbm(binary_matrix& A, const char* fname);

#endif

