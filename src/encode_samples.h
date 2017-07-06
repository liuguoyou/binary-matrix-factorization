#ifndef ENCODE_SAMPLES_H
#define ENCODE_SAMPLES_H
#include "binmat.h"

/** 
 * Given a current error E, dictionary D and coefficients A, update coefficients and error, 
 * so that the total weight of the error E is reduced. The algorithm updates one row of A at a time.
 *
 * E = AD is an n x m matrix
 * D is a p x m matrix, where each row is an atom of dimension m
 * A is a n x p matrix, where each row contains the coefficients for representing the corresponding row of X=AD+E 
 * max_a_weight is the maximum allowed number of ones in any row of A
 * max_e_weight is the maximum allowed number of ones in any row of E
*/
idx_t encode_samples_basic(binary_matrix& E,
			   const binary_matrix& H,
			   const binary_matrix& D,
			   binary_matrix& A,
			   const idx_t max_a_weight,
			   const idx_t max_e_weight);

/** 
 * Given a current error E, dictionary D and coefficients A, update coefficients and error, 
 * so that the total weight of the error E is reduced. The algorithm updates one row of A at a time.
 *
 * E = AD is an n x m matrix
 * D is a p x m matrix, where each row is an atom of dimension m
 * A is a n x p matrix, where each row contains the coefficients for representing the corresponding row of X=AD+E 
 * max_a_weight is the maximum allowed number of ones in any row of A
 * max_e_weight is the maximum allowed number of ones in any row of E
*/
idx_t encode_samples_omp(binary_matrix& E,
			 const binary_matrix& H,
			 const binary_matrix& D,
			 binary_matrix& A,
			 const idx_t max_a_weight,
			 const idx_t max_e_weight);

#endif
