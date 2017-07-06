#ifndef UPDATE_DICTIONARY_H
#define UPDATE_DICTIONARY_H

#include "binmat.h"

/** 
 * Given a current error E = X - AD, dictionary D and coefficients A, update dictionary and error, 
 * so  that the total weight of the error E is reduced.
 * @param E = AD is an n x m matrix
 * @param H is an n x m  mask where 0 indicates a datum in E is missing in the corresponding place
 * @param D is a p x m matrix, where each row is an atom of dimension m
 * @param A is a n x p matrix, where each row contains the coefficients for representing the corresponding row of X=AD+E
 */
idx_t update_dictionary_steepest(binary_matrix& E,
				 const binary_matrix& H,
				 binary_matrix& D,
				 binary_matrix& A);


/** 
 * Given a current error E = X - AD, dictionary D and coefficients A, update dictionary and error, 
 * so  that the total weight of the error E is reduced.
 * @param E = AD is an n x m matrix
 * @param H is an n x m  mask where 0 indicates a datum in E is missing in the corresponding place
 * @param D is a p x m matrix, where each row is an atom of dimension m
 * @param A is a n x p matrix, where each row contains the coefficients for representing the corresponding row of X=AD+E
 */
idx_t update_dictionary_steepest_omp(binary_matrix& E,
				 const binary_matrix& H,
				     binary_matrix& D,
				     binary_matrix& A);

/**
 * dictionary update using PROXIMUS step, which is in some way similar to the K-SVD
 * rank one update step.
ror, 
 * @param E = AD is an n x m matrix
 * @param H is an n x m  mask where 0 indicates a datum in E is missing in the corresponding place
 * @param D is a p x m matrix, where each row is an atom of dimension m
 * @param A is a n x p matrix, where each row contains the coefficients for representing the corresponding row of X=AD+E
 */
idx_t update_dictionary_proximus(binary_matrix& E,
				 const binary_matrix& H,
				 binary_matrix& D,
				 binary_matrix& A);

/**
 * dictionary update using PROXIMUS step, which is in some way similar to the K-SVD
 * rank one update step.
ror, 
 * @param E = AD is an n x m matrix
 * @param H is an n x m  mask where 0 indicates a datum in E is missing in the corresponding place
 * @param D is a p x m matrix, where each row is an atom of dimension m
 * @param A is a n x p matrix, where each row contains the coefficients for representing the corresponding row of X=AD+E
 */
idx_t update_dictionary_proximus_omp(binary_matrix& E,
				     const binary_matrix& H,
				     binary_matrix& D,
				     binary_matrix& A);

#endif
