#ifndef INITIALIZE_DICTIONARY_H
#define INITIALIZE_DICTIONARY_H

#include "binmat.h"

void initialize_dictionary_neighbor(const binary_matrix& E, 
				    const binary_matrix& H,
				    binary_matrix& D, 
				    binary_matrix& A);

void initialize_dictionary_random_centroids(const binary_matrix& E, 
					    const binary_matrix& H,
					    binary_matrix& D, 
					    binary_matrix& A);

void initialize_dictionary_random_centroids_xor(const binary_matrix& E, 
						const binary_matrix& H,
						binary_matrix& D, 
						binary_matrix& A);

void initialize_dictionary_graph_grow(const binary_matrix& E, 
				      const binary_matrix& H,
				      binary_matrix& D, 
				      binary_matrix& A);

void initialize_dictionary_partition(const binary_matrix& E, 
				     const binary_matrix& H,
				     binary_matrix& D, 
				     binary_matrix& A);

void initialize_dictionary_random(const binary_matrix& E,
				  const binary_matrix& H,
				  binary_matrix& D,
				  binary_matrix& A);

//#define initialize_dictionary initialize_dictionary_partition //  GOOD
//#define initialize_dictionary initialize_dictionary_neighbor  // BEST SO FAR
//#define initialize_dictionary initialize_dictionary_graph_grow VERY SLOW AND DOES NOT WORK WELL, DONT KNOW WHY
//#define initialize_dictionary initialize_dictionary_random_centroids WORST

#endif
