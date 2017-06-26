#ifndef BSVD_H
#define BSVD_H

idx_t learn_model_traditional(binary_matrix& X,
			      binary_matrix& E, 
			      binary_matrix& D, 
			      binary_matrix& A);

idx_t learn_model_alter1(binary_matrix& X,
			 binary_matrix& E, 
			 binary_matrix& D, 
			 binary_matrix& A);

idx_t learn_model_alter2(binary_matrix& X,
			 binary_matrix& E, 
			 binary_matrix& D, 
			 binary_matrix& A);

idx_t learn_model_alter3(binary_matrix& X,
			 binary_matrix& E, 
			 binary_matrix& D, 
			 binary_matrix& A);

idx_t learn_model_mdl_forward_selection(binary_matrix& X,
					binary_matrix& E, 
					binary_matrix& D, 
					binary_matrix& A);

idx_t learn_model_mdl_backward_selection(binary_matrix& X,
					 binary_matrix& E, 
					 binary_matrix& D, 
					 binary_matrix& A);

idx_t learn_model_mdl_full_search(binary_matrix& X,
				  binary_matrix& E, 
				  binary_matrix& D, 
				  binary_matrix& A);

typedef void (*mi_algorithm_t)(const binary_matrix& E, 
				binary_matrix& D, 
				binary_matrix& A);

typedef idx_t (*cu_algorithm_t)(binary_matrix& E, 
				const binary_matrix& D, 
				binary_matrix& A);

typedef idx_t (*du_algorithm_t)(binary_matrix& E, 
				binary_matrix& D, 
				binary_matrix& A);

typedef idx_t (*ml_algorithm_t)(binary_matrix& X,
				binary_matrix& E, 
				binary_matrix& D, 
				binary_matrix& A);

extern mi_algorithm_t initialize_dictionary;
extern cu_algorithm_t encode_samples;
extern du_algorithm_t update_dictionary;
extern ml_algorithm_t learn_model;
extern ml_algorithm_t learn_model_inner;

extern const char* mi_algorithm_names[];
extern const char* cu_algorithm_names[];
extern const char* du_algorithm_names[];
extern const char* lm_algorithm_names[];


void learn_model_setup(int mi_algo, int cu_algo, int du_algo, int lm_algo, int lmi_algo);

#endif
