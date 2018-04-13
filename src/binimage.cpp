#include <cassert>
#include <cmath>
#include "binimage.h"
#include "intmat.h"

size_t compute_grid_size(const size_t size, const size_t width, const size_t stride, extract_t extract_type) {
  // return M, the number of patches to be extracted along a dimension
  // EXTRACT_EXACT: last index (li) must be < m so
  // li = (M-1)*stride + width -1 <= m -1 => M = floor [ (m - width + stride) / stride ]
  // EXTRACT_FULL: first index of last patch (fi) must be < m so
  // fi = (M-1)*stride  <= m - 1 => stride*M <= m + stride - 1 => M = floor [(m -1 + stride) / stride ]
  return extract_type == EXTRACT_EXACT ? (size+stride-width)/stride : (size + stride - 1)/stride;
}


void extract_patches(const binary_matrix& I,
		     const size_t width,
		     const size_t stride,
		     const extract_t e,
		     binary_matrix& X) {
  const size_t M = I.get_rows();
  const size_t N = I.get_cols();
  const size_t mg = compute_grid_size(M,width,stride,e);
  const size_t ng = compute_grid_size(N,width,stride,e);
  const size_t n = mg*ng;
  const size_t m = width*width;
  assert(X.get_rows() == n);
  assert(X.get_cols() == m);
  binary_matrix P(width,width),V(1,m);
  size_t li = 0;
  X.clear(); // debug
  for (size_t i = 0, ig = 0; ig < mg; ig++, i += stride) { 
    for (size_t j = 0, jg = 0; jg < ng; jg++, j += stride) {
      I.copy_submatrix_to(i,i+width,j,j+width,P);
      P.copy_vectorized_to(V);
      X.set_row(li,V);
      li++;
    } // for j
  } // for i
  P.destroy();
  V.destroy();
}

void stitch_patches(const binary_matrix& P,
		    const size_t M,
		    const size_t N,
		    const size_t stride,
		    const extract_t e,
		    binary_matrix& I,
		    integer_matrix& A,
		    integer_matrix& C) {
  assert(I.get_rows() == M);
  assert(I.get_cols() == N);
  assert(A.get_rows() == M);
  assert(A.get_cols() == N); 
  assert(C.get_rows() == M);
  assert(C.get_cols() == N); 
  const size_t width = (size_t) sqrt(double(P.get_cols()));
  const size_t mg = compute_grid_size(M,width,stride,e);
  const size_t ng = compute_grid_size(N,width,stride,e);
  const size_t n = mg*ng;
  assert(P.get_rows() == n);
  C.clear(); // clear image
  A.clear(); // clear image
  size_t k = 0;
  for (size_t ig = 0, i = 0; ig < mg; ++ig, i+= stride ) {
    for (size_t jg = 0, j = 0; jg < ng; ++jg, j += stride) {
      size_t l = 0;
      for (size_t ii = 0; ii < width; ++ii) {
	for (size_t jj = 0; jj < width; ++jj) {
	  C.inc(i+ii,j+jj);
	  if (P.get(k,l))
	    A.inc(i+ii,j+jj);
	  l++;
	} // for ii
      }// for jj
      k++;
    } // for j
  } // for i
  //
  // normalization
  //
  for (size_t i = 0; i < M; ++i ) {
    for (size_t j = 0; j < N; ++j ) {
      I.set(i, j, 2*A.get(i,j) >= C.get(i,j) );
    } // for j
  } // for i
}
