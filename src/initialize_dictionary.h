#ifndef INITIALIZE_DICTIONARY_H
#define INITIALIZE_DICTIONARY_H

#include "binmat.h"



/**
 * Each atom $d_k$ is initialized using a randomly drawn sample from $X$,
 * $c_k$.  Then $d_k$ is computed as the Hamming average of all samples $x$
 * so that $c_k\transp x > 0$, that is, so that they have at least one
 * dimension where both $c_k$ and $x$ are $1$. This is a
 * lame method for very fat matrices, as there will be so many of such
 * samples, and if noise is present. Some more intelligent criterion must be
 * devised for the dictionary learning case.
 */
void initialize_dictionary_neighbor(const binary_matrix& E, 
				    const binary_matrix& H,
				    binary_matrix& D, 
				    binary_matrix& A);

/**
 * Works only when $K \leq M$. This implementation ranks the dimensions of
 * the data matrix $X$ (that is, the number of columns of $X$) in descending
 * weight order. Then, each atom $k$ is initialized as the Hamming average of
 * all samples in $X$ which have a value of $1$ in the $k$-th ranked
 * dimension.
 */
void initialize_dictionary_partition(const binary_matrix& E, 
				     const binary_matrix& H,
				     binary_matrix& D, 
				     binary_matrix& A);


/** Sign of the DCT  **/
void initialize_dictionary_sdct(const binary_matrix& E, 
				const binary_matrix& H,
				binary_matrix& D, 
				binary_matrix& A);


void initialize_dictionary_sdct2d(const binary_matrix& E, 
				  const binary_matrix& H,
				  binary_matrix& D, 
				  binary_matrix& A);

/** 
 * Hamming overcomplete basis: columns are the binary representations of '1','2','etc. We can get up to 2^n elements **/
void initialize_dictionary_hamming(const binary_matrix& E, 
				  const binary_matrix& H,
				  binary_matrix& D, 
				   binary_matrix& A);


/**
 * Each atom is initialized with samples from a Bernoulli(1/2) IID process
 */
void initialize_dictionary_random(const binary_matrix& E,
				  const binary_matrix& H,
				  binary_matrix& D,
				  binary_matrix& A);

/**
 * set expected proportion of ones in each random atom
 */
void set_random_dictionary_density(double a);

//#define initialize_dictionary initialize_dictionary_partition //  GOOD
//#define initialize_dictionary initialize_dictionary_neighbor  // BEST SO FAR
//#define initialize_dictionary initialize_dictionary_graph_grow VERY SLOW AND DOES NOT WORK WELL, DONT KNOW WHY
//#define initialize_dictionary initialize_dictionary_random_centroids WORST


#endif
