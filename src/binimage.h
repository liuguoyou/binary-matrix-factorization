#ifndef BINIMAGE_H
#define BINIMAGE_H
#include <stdlib.h>
#include "binmat.h"
#include "intmat.h"

#define CCLIP(x,a,b) ( (x) > (a) ? ( (x) < (b) ? (x) : (b) ) : (a) )

typedef enum {
  EXTRACT_EXACT, /// only extract patches which contain true pixels, possibly leaving bordering pixels out
  EXTRACT_FULL, ///  extract patches so that whole image is covered, extrapolating border pixels as needed
} extract_t;

size_t compute_grid_size(const size_t size, const size_t width, const size_t stride, extract_t extract_type);


void extract_patches(const binary_matrix& I,
		     const size_t width,
		     const size_t stride,
		     const extract_t e,
		     binary_matrix& P);


// I must be preallocated 
void stitch_patches(const binary_matrix& P,
		    const size_t M,
		    const size_t N,
		    const size_t stride,
		    const extract_t e,
		    binary_matrix& I,
		    integer_matrix& A,
		    integer_matrix& C);
#endif
