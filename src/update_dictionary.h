#ifndef UPDATE_DICTIONARY_H
#define UPDATE_DICTIONARY_H

#include "binmat.h"

/** 
 * Given a current error E = X - AD, dictionary D and coefficients A, update dictionary and error, 
 * so  that the total weight of the error E is reduced.
 * E = AD is an n x m matrix
 * D is a p x m matrix, where each row is an atom of dimension m
 * A is a n x p matrix, where each row contains the coefficients for representing the corresponding row of X=AD+E
 */
idx_t update_dictionary_steepest(binary_matrix& E, binary_matrix& D, binary_matrix& A);
idx_t update_dictionary_steepest_omp(binary_matrix& E, binary_matrix& D, binary_matrix& A);

/**
 * dictionary update using PROXIMUS step
 */
idx_t update_dictionary_proximus(binary_matrix& E, binary_matrix& D, binary_matrix& A);

/**
 * dictionary update using PROXIMUS step; multiprocessor enabled
 */
idx_t update_dictionary_proximus_omp(binary_matrix& E, binary_matrix& D, binary_matrix& A);

#endif
