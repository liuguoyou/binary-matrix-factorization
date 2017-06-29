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

#endif
