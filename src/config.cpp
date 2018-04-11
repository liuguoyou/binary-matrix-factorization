#include "config.h"
#include <cstdlib>
#include "bsvd.h"
#include "initialize_dictionary.h"
#include "coefficients_update.h"
#include "update_dictionary.h"

mi_algorithm_t initialize_dictionary = initialize_dictionary_neighbor;
es_algorithm_t coefficients_update = coefficients_update_omp;
du_algorithm_t update_dictionary = update_dictionary;
ml_algorithm_t learn_model = learn_model_traditional;
ml_algorithm_t learn_model_inner = learn_model_traditional;

mi_algorithm_t mi_algorithm_catalog[] = {initialize_dictionary_neighbor,
					 initialize_dictionary_partition,
					 initialize_dictionary_random_centroids,
					 initialize_dictionary_random_centroids_xor,
					 initialize_dictionary_graph_grow,
					 initialize_dictionary_random,
					 0};

const char* mi_algorithm_names[] = {"Neighbor initialization",
				    "Partition initialization",
				    "Random centroids initialization",
				    "Random centroids (in mod-2 algebra) initialization",
				    "Graph growing initialization",
				    "purely random dictionary",
				    0
};

es_algorithm_t es_algorithm_catalog[] = {coefficients_update_omp,
					 coefficients_update_basic,
					 coefficients_update_corr,
					 0};

const char* es_algorithm_names[] = {"OpenMP basic coefficients update",
				    "Basic coefficients update",
				    "Coefficients update using correlation",
				    0
};

du_algorithm_t du_algorithm_catalog[] = {update_dictionary_steepest,
					 update_dictionary_proximus,
					 update_dictionary_steepest_omp,
					 update_dictionary_proximus_omp,
					 0};

const char* du_algorithm_names[] = {"Steepest descent (a la MOD)  dictionary update",
				    "Proximus-like dictionary update",
				    "Steepest descent (a la MOD)  dictionary update (OMP)",
				    "Proximus-like dictionary update (OMP)",0
};

ml_algorithm_t learn_model_algorithm_catalog[] = {learn_model_traditional,
						  learn_model_alter1,
						  learn_model_alter2,
						  learn_model_alter3,
						  learn_model_mdl_forward_selection,
						  learn_model_mdl_backward_selection,
						  learn_model_mdl_full_search,
						  0};
const char* lm_algorithm_names[] = {"Model learning by traditional alternate descent",
				    "Role-switching learning 1: at each iteration, the role of A and D are switched",
				    "Role-switched learning 2: after convergence, the role of A and D are switched and traditional model is applied again",
				    "Role switched learning 3: like RS1 but only update_dictionary is applied (for use with Proximus",
				    "MDL/forward selection",
				    "MDO/backward selection",
				    "MDL/full search"
};


void learn_model_setup(int mi_algo, int es_algo, int du_algo, int lm_algo, int lmi_algo) {
  if (mi_algo > 5) { std::cerr << "Invalid model initialization algorithm (0-" << 5 << ')' << std::endl; exit(-1); }
  if (es_algo > 2) { std::cerr << "Invalid coefficients update algorithm (0-" << 2 << ')' << std::endl; exit(-1); }
  if (du_algo > 3) { std::cerr << "Invalid dictionary update algorithm (0-" << 3 << ')' << std::endl; exit(-1); }
  if (lm_algo > 6) { std::cerr << "Invalid model learning algorithm (0-" << 6 << ')' << std::endl; exit(-1); }
  if (lmi_algo > 3) { std::cerr << "Invalid inner model learning algorithm (0-" << 3 << ')' << std::endl; exit(-1); }

  initialize_dictionary = mi_algorithm_catalog[mi_algo];
  std::cout << "Using " << mi_algorithm_names[mi_algo] << std::endl;
  coefficients_update = es_algorithm_catalog[es_algo];
  std::cout << "Using " << es_algorithm_names[es_algo] << std::endl;
  update_dictionary = du_algorithm_catalog[du_algo];
  std::cout << "Using " << du_algorithm_names[du_algo] << std::endl;
  learn_model = learn_model_algorithm_catalog[lm_algo];
  std::cout << "Using " << lm_algorithm_names[lm_algo] << " for outer learning loop." << std::endl;
  learn_model_inner = learn_model_algorithm_catalog[lmi_algo];
  std::cout << "Using " << lm_algorithm_names[lmi_algo] << " for inner learning." << std::endl;
}
