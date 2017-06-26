#ifndef ENCODE_SAMPLES_H
#define ENCODE_SAMPLES_H

#include "binmat.h"
#include <omp.h>


/** 
 * Given a current error E, dictionary D and coefficients A, update coefficients and error, 
 * so that the total weight of the error E is reduced. The algorithm updates one row of A at a time.
 *
 * E = AD is an n x m matrix
 * D is a p x m matrix, where each row is an atom of dimension m
 * A is a n x p matrix, where each row contains the coefficients for representing the corresponding row of X=AD+E 
*/
idx_t encode_samples_basic(binary_matrix& E, const binary_matrix& D, binary_matrix& A);

idx_t encode_samples_omp(binary_matrix& E, const binary_matrix& D, binary_matrix& A);

/** 
 * Given a current error E, dictionary D and coefficients A, update coefficients and error, 
 * so that the total weight of the error E is reduced. The algorithm updates one row of A at a time.
 *
 * E = AD is an n x m matrix
 * D is a p x m matrix, where each row is an atom of dimension m
 * A is a n x p matrix, where each row contains the coefficients for representing the corresponding row of X=AD+E 
*/
idx_t encode_samples_fast(binary_matrix& E, const binary_matrix& D, binary_matrix& A);

#endif
