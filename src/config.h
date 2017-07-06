#ifndef CONFIG_H
#define CONFIG_H
#include "binmat.h"

typedef void (*mi_algorithm_t)(const binary_matrix& E, 
			       const binary_matrix& H,
			       binary_matrix& D, 
			       binary_matrix& A);

typedef idx_t (*es_algorithm_t)(binary_matrix& E, 
				const binary_matrix& H,
				const binary_matrix& D, 
				binary_matrix& A,
				const idx_t max_a_weight,
				const idx_t max_e_weight);

typedef idx_t (*du_algorithm_t)(binary_matrix& E, 
				const binary_matrix& H,
				binary_matrix& D, 
				binary_matrix& A);

typedef idx_t (*ml_algorithm_t)(binary_matrix& X,
				const binary_matrix& H,
				binary_matrix& E, 
				binary_matrix& D, 
				binary_matrix& A);

extern mi_algorithm_t initialize_dictionary;
extern es_algorithm_t encode_samples;
extern du_algorithm_t update_dictionary;
extern ml_algorithm_t learn_model;
extern ml_algorithm_t learn_model_inner;

extern const char* mi_algorithm_names[];
extern const char* es_algorithm_names[];
extern const char* du_algorithm_names[];
extern const char* lm_algorithm_names[];


void learn_model_setup(int mi_algo, int es_algo, int du_algo, int lm_algo, int lmi_algo);

#endif
