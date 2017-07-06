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
 * Given a number $K$ of atoms, each data sample is randomly assigned to one
 * atom. The $K$ atoms are then the average (element-wise majority vote) of
 * all the samples to which they were assigned.
 * I don't think this kind of initialization makes ANY sense.
 */
void initialize_dictionary_random_centroids(const binary_matrix& E, 
					    const binary_matrix& H,
					    binary_matrix& D, 
					    binary_matrix& A);

/**
 * Given a number $K$ of atoms, each data sample is randomly assigned to one
 * atom. The $K$ atoms are then the element-wise parity (XOR) of all the
 * samples to which they were assigned. I don't think this kind of
 * initialization makes ANY sense either.
 */
void initialize_dictionary_random_centroids_xor(const binary_matrix& E, 
						const binary_matrix& H,
						binary_matrix& D, 
						binary_matrix& A);
/**
 * This is actually a variant of the original graph grow algorithm. Given
 * $K$ as the number of atoms to initialize, the method starts by drawing $K$
 * initial (different) samples from the data set $X$ and associates a
 * \emph{weight} vector $w$ to it, which is initialized to the value of the
 * corresponding random sample. Then it randomly draws one previously unused
 * sample $x$ from $X$ and adds it to the weight vector $w\opt$ which
 * satisfies $$(w\opt)\transp x/\|w\opt\|_1 \geq w\transp x/\|w\|_1.$$ This
 * procedure continues until all data samples in $X$ have been used. The
 * final $k$-th atom is the Hamming average of all the samples added to its
 * corresponding weight vector, that is, $d_k[i] = \lfloor w_k[i] / n_k
 * \rceil$, where $n_k$ is the number of samples from $X$ which were added to
 * $w_k$. This is a \emph{soft} version of the original graph grow algorithm,
 * where the criterion used to assign a sample to a given atom $k$ is based
 * on the Hamming distance between $w$ and a sample $x$.
 *
 */
void initialize_dictionary_graph_grow(const binary_matrix& E, 
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
/**
 * Each atom is initialized with samples from a Bernoulli(1/2) IID process
 */
void initialize_dictionary_random(const binary_matrix& E,
				  const binary_matrix& H,
				  binary_matrix& D,
				  binary_matrix& A);

//#define initialize_dictionary initialize_dictionary_partition //  GOOD
//#define initialize_dictionary initialize_dictionary_neighbor  // BEST SO FAR
//#define initialize_dictionary initialize_dictionary_graph_grow VERY SLOW AND DOES NOT WORK WELL, DONT KNOW WHY
//#define initialize_dictionary initialize_dictionary_random_centroids WORST

#endif
