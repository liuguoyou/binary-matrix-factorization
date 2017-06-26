#include <algorithm>
#include "gsl/gsl_randist.h"
#include <omp.h>
#include "util.h"
#include <iomanip>
#include "bsvd.h"

static gsl_rng* get_rng() {
  static gsl_rng* rng = 0;
  if (rng == 0) { 
    rng = gsl_rng_alloc (gsl_rng_rand48);
    gsl_rng_set (rng, random_seed); // set random_seed
  } 
  return rng;
}

mi_algorithm_t initialize_dictionary = initialize_dictionary_neighbor;
cu_algorithm_t encode_samples = update_coefficients_omp;
du_algorithm_t update_dictionary = update_dictionary;
ml_algorithm_t learn_model = learn_model_traditional;
ml_algorithm_t learn_model_inner = learn_model_traditional;

long random_seed = 34503498;

mi_algorithm_t mi_algorithm_catalog[] = {initialize_dictionary_neighbor,
					 initialize_dictionary_partition,
					 initialize_dictionary_random_centroids,
					 initialize_dictionary_random_centroids_xor,
					 initialize_dictionary_graph_grow,
					 0};

const char* mi_algorithm_names[] = {"Neighbor initialization",
				    "Partition initialization",
				    "Random centroids initialization",
				    "Random centroids (in mod-2 algebra) initialization",
				    "Graph growing initialization",0
};

cu_algorithm_t cu_algorithm_catalog[] = {encode_samples_omp,
					 encode_samples,
					 encode_samples_fast, // broken
					 0};

const char* cu_algorithm_names[] = {"OpenMP basic coefficients update",
				    "Basic coefficients update",
				    "Fast coefficients update (broken!)",0
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


void learn_model_setup(int mi_algo, int cu_algo, int du_algo, int lm_algo, int lmi_algo) {
  if (mi_algo > 4) { std::cerr << "Invalid model initialization algorithm (0-" << 4 << ')' << std::endl; exit(-1); }
  if (cu_algo > 2) { std::cerr << "Invalid coefficients update algorithm (0-" << 2 << ')' << std::endl; exit(-1); }
  if (du_algo > 3) { std::cerr << "Invalid dictionary update algorithm (0-" << 3 << ')' << std::endl; exit(-1); }
  if (lm_algo > 6) { std::cerr << "Invalid model learning algorithm (0-" << 6 << ')' << std::endl; exit(-1); }
  if (lmi_algo > 3) { std::cerr << "Invalid inner model learning algorithm (0-" << 3 << ')' << std::endl; exit(-1); }

  initialize_dictionary = mi_algorithm_catalog[mi_algo];
  std::cout << "Using " << mi_algorithm_names[mi_algo] << std::endl;
  encode_samples = cu_algorithm_catalog[cu_algo];
  std::cout << "Using " << cu_algorithm_names[cu_algo] << std::endl;
  update_dictionary = du_algorithm_catalog[du_algo];
  std::cout << "Using " << du_algorithm_names[du_algo] << std::endl;
  learn_model = learn_model_algorithm_catalog[lm_algo];
  std::cout << "Using " << lm_algorithm_names[lm_algo] << " for outer learning loop." << std::endl;
  learn_model_inner = learn_model_algorithm_catalog[lmi_algo];
  std::cout << "Using " << lm_algorithm_names[lmi_algo] << " for inner learning." << std::endl;
}

idx_t learn_model_traditional(binary_matrix& X,
			      binary_matrix& E, 
			      binary_matrix& D, 
			      binary_matrix& A) {
  mul(A,false,D,false,E);
  add(E,X,E);
  //  std::cout << "trad" << std::endl;
  idx_t changed = 1;
  idx_t iter = 0;
  //std::cout << "iter=" << std::setw(8) << iter 
//	    << "\t||E||=" << std::setw(8) << E.weight()
//	    << "\t||D||=" << std::setw(8) << D.weight()
//	    << "\t||A||=" << std::setw(8) << A.weight() << std::endl;
  while (changed > 0) {    
    iter++;
    idx_t changed_coefs = encode_samples(E,D,A);
  //  std::cout << "iter=" << std::setw(8) << iter 
//	      << "\t||E||=" << std::setw(8) << E.weight()
//	      << "\t||D||=" << std::setw(8) << D.weight()
//	      << "\t||A||=" << std::setw(8) << A.weight()
//	      << "\tchanged coefs=" << std::setw(8) << changed_coefs << std::endl;
    changed = changed_coefs + update_dictionary(E,D,A);
//    std::cout << "iter=" << std::setw(8) << iter 
//	      << "\t||E||=" << std::setw(8) << E.weight()
//	      << "\t||D||=" << std::setw(8) << D.weight()
//	      << "\t||A||=" << std::setw(8) << A.weight()
//	      << "\tchanged atoms=" << std::setw(8) << changed << std::endl;
  }
  return iter;
}


idx_t learn_model_alter1(binary_matrix& X,
			 binary_matrix& E, 
			 binary_matrix& D, 
			 binary_matrix& A) {
  //  std::cout << "alter1" << std::endl;
  const idx_t N = E.get_rows();
  const idx_t M = E.get_cols();
  const idx_t K = D.get_rows();

  mul(A,false,D,false,E);
  add(E,X,E);
  binary_matrix Dt(M,K);
  binary_matrix At(K,N);
  binary_matrix Et(M,N);

  //
  // RUN BSVD
  //
  idx_t changed = 1;
  idx_t iter = 0;
//  std::cout << "iter=" << std::setw(8) << iter 
//	    << "\t||E||=" << std::setw(8) << E.weight()
//	    << "\t||D||=" << std::setw(8) << D.weight()
//	    << "\t||A||=" << std::setw(8) << A.weight() << std::endl;
  while (changed > 0) {    
    iter++;
    idx_t changed_coefs = encode_samples(E,D,A);
//    std::cout << "DIRECT: iter=" << std::setw(8) << iter 
//	      << "\t||E||=" << std::setw(8) << E.weight()
//	      << "\t||D||=" << std::setw(8) << D.weight()
//	      << "\t||A||=" << std::setw(8) << A.weight()
//	      << "\tchanged coefs=" << std::setw(8) << changed_coefs << std::endl;
    changed = changed_coefs + update_dictionary(E,D,A);
  //  std::cout << "DIRECT: iter=" << std::setw(8) << iter 
//	      << "\t||E||=" << std::setw(8) << E.weight()
//	      << "\t||D||=" << std::setw(8) << D.weight()
//	      << "\t||A||=" << std::setw(8) << A.weight()
//	      << "\tchanged atoms=" << std::setw(8) << changed << std::endl;

    A.transpose_to(At);
    D.transpose_to(Dt);
    E.transpose_to(Et);

    changed_coefs = encode_samples(Et,At,Dt);
  //  std::cout << "TRANSP: iter=" << std::setw(8) << iter 
//	      << "\t||E||=" << std::setw(8) << Et.weight()
//	      << "\t||D||=" << std::setw(8) << Dt.weight()
//	      << "\t||A||=" << std::setw(8) << At.weight()
//	      << "\tchanged coefs=" << std::setw(8) << changed_coefs << std::endl;

    changed = update_dictionary(Et,At,Dt);
  //  std::cout << "TRANSP: iter=" << std::setw(8) << iter 
//	      << "\t||E||=" << std::setw(8) << Et.weight()
//	      << "\t||D||=" << std::setw(8) << Dt.weight()
//	      << "\t||A||=" << std::setw(8) << At.weight()
//	      << "\tchanged atoms=" << std::setw(8) << changed << std::endl;
    
    At.transpose_to(A);
    Dt.transpose_to(D);
    Et.transpose_to(E);
  }
  At.destroy();
  Dt.destroy();
  Et.destroy();
  return iter;
}


idx_t learn_model_alter2(binary_matrix& X,
			 binary_matrix& E, 
			 binary_matrix& D, 
			 binary_matrix& A) {
  //  std::cout << "alter2" << std::endl;
  const idx_t N = E.get_rows();
  const idx_t M = E.get_cols();
  const idx_t K = D.get_rows();

  mul(A,false,D,false,E);
  add(E,X,E);
  binary_matrix Dt(M,K);
  binary_matrix At(K,N);
  binary_matrix Et(M,N);
  idx_t changed = 1;
  idx_t iter = 0;
//  std::cout << "iter=" << std::setw(8) << iter 
//	    << "\t||E||=" << std::setw(8) << E.weight()
//	    << "\t||D||=" << std::setw(8) << D.weight()
//	    << "\t||A||=" << std::setw(8) << A.weight() << std::endl;

  idx_t outer_changed = 1;
  while (outer_changed > 0) {
    outer_changed = 0;
    while (changed > 0) {    
      iter++;
      idx_t changed_coefs = encode_samples(E,D,A);
//      std::cout << "DIRECT: iter=" << std::setw(8) << iter 
//		<< "\t||E||=" << std::setw(8) << E.weight()
//		<< "\t||D||=" << std::setw(8) << D.weight()
//		<< "\t||A||=" << std::setw(8) << A.weight()
//		<< "\tchanged coefs=" << std::setw(8) << changed_coefs << std::endl;
      changed = changed_coefs + update_dictionary(E,D,A);
//      std::cout << "DIRECT: iter=" << std::setw(8) << iter 
//		<< "\t||E||=" << std::setw(8) << E.weight()
//		<< "\t||D||=" << std::setw(8) << D.weight()
//		<< "\t||A||=" << std::setw(8) << A.weight()
//		<< "\tchanged atoms=" << std::setw(8) << changed << std::endl;
      outer_changed += changed;
    }
    A.transpose_to(At);
    D.transpose_to(Dt);
    E.transpose_to(Et);
    changed = 1;
    iter = 0;
    while (changed > 0) {    
      iter++;
      idx_t changed_coefs = encode_samples(Et,At,Dt);
//      std::cout << "TRANSPOSED: iter=" << std::setw(8) << iter 
//		<< "\t||E||=" << std::setw(8) << Et.weight()
//		<< "\t||D||=" << std::setw(8) << Dt.weight()
//		<< "\t||A||=" << std::setw(8) << At.weight()
//		<< "\tchanged coefs=" << std::setw(8) << changed_coefs << std::endl;

      changed = changed_coefs + update_dictionary(Et,At,Dt);
//      std::cout << "TRANSPOSED: iter=" << std::setw(8) << iter 
//		<< "\t||E||=" << std::setw(8) << Et.weight()
//		<< "\t||D||=" << std::setw(8) << Dt.weight()
//		<< "\t||A||=" << std::setw(8) << At.weight()
//		<< "\tchanged atoms=" << std::setw(8) << changed << std::endl; 
      outer_changed += changed;   
    }
    At.transpose_to(A);
    Dt.transpose_to(D);
    Et.transpose_to(E);
  }
  At.destroy();
  Dt.destroy();
  Et.destroy();
  return iter;
}


idx_t learn_model_alter3(binary_matrix& X,
			 binary_matrix& E, 
			 binary_matrix& D, 
			 binary_matrix& A) {
  //  std::cout << "alter3" << std::endl;
  const idx_t N = E.get_rows();
  const idx_t M = E.get_cols();
  const idx_t K = D.get_rows();
  mul(A,false,D,false,E);
  add(E,X,E);

  binary_matrix Dt(M,K);
  binary_matrix At(K,N);
  binary_matrix Et(M,N);
  idx_t changed = K+1;
  idx_t iter = 0;
//  std::cout << "iter=" << std::setw(8) << iter 
//	    << "\t||E||=" << std::setw(8) << E.weight()
//	    << "\t||D||=" << std::setw(8) << D.weight()
//	    << "\t||A||=" << std::setw(8) << A.weight() << std::endl;
  while (changed > 0) {    
    iter++;
    A.transpose_to(At);
    D.transpose_to(Dt);
    E.transpose_to(Et);
    changed = update_dictionary(Et,At,Dt);
//    std::cout << "iter=" << std::setw(8) << iter 
//	      << "\t||E||=" << std::setw(8) << Et.weight()
//	      << "\t||D||=" << std::setw(8) << Dt.weight()
//	      << "\t||A||=" << std::setw(8) << At.weight()
//	      << "\tchanged atoms=" << std::setw(8) << changed << std::endl;

    At.transpose_to(A);
    Dt.transpose_to(D);
    Et.transpose_to(E);
    changed = update_dictionary(E,D,A);
//    std::cout << "iter=" << std::setw(8) << iter 
//	      << "\t||E||=" << std::setw(8) << E.weight()
//	      << "\t||D||=" << std::setw(8) << D.weight()
//	      << "\t||A||=" << std::setw(8) << A.weight()
//	      << "\tchanged atoms=" << std::setw(8) << changed << std::endl;
  }
  At.destroy();
  Dt.destroy();
  Et.destroy();
return iter;
}

#include "coding.h"

idx_t model_codelength(const binary_matrix& E, 
		       const binary_matrix& D, 
		       const binary_matrix& A) {

  const idx_t M = E.get_cols();
  const idx_t N = E.get_rows();
  const idx_t K = D.get_rows();


  binary_matrix Dk(1,M);
  binary_matrix Ak(1,N);

  idx_t LE = universal_codelength(E.get_rows()*E.get_cols(),E.weight());
  idx_t LD = 0, LA = 0;
  for (idx_t k = 0; k < K; k++) {
    D.copy_row_to(k,Dk);
    A.copy_col_to(k,Ak);
    LD += universal_codelength(M,Dk.weight());
    LA += universal_codelength(N,Ak.weight());
  }
  Dk.destroy();
  Ak.destroy();
  return LE+LD+LA;
}

idx_t learn_model_mdl_forward_selection(binary_matrix& X,
					binary_matrix& E, 
					binary_matrix& D, 
					binary_matrix& A) {
  const idx_t M = E.get_cols();
  const idx_t N = E.get_rows();
  idx_t K = D.get_rows();
  learn_model_inner(X,E,D,A);
  binary_matrix nextAtom(1,M);
  binary_matrix nextCoefs(N,1);
  binary_matrix currD(D),currA(A),currE(E);  
  binary_matrix nextD,nextA,nextE(N,M);
  idx_t bestK = K;
  idx_t bestL = model_codelength(E,D,A);
  idx_t stuck = 0;
  idx_t sumStuck = 0;
  idx_t allStuck = 0;
  do {
    //
    // initialize curr atom and associated coefs.
    //
    idx_t currL = model_codelength(currE,currD,currA);
    int dif = int(currL) - int(bestL);
    int dev = allStuck > 0 ? (sumStuck/allStuck) : 0;
    std::cout << "currK=" << K << " currL=" << currL << " bestK=" << bestK << " bestL=" << bestL << " stuck=" << stuck << " dif=" << dif << " dev=" << dev << std::endl;
    initialize_dictionary(currE,nextAtom,nextCoefs);
    //learn_model_inner(currE,nextE,nextAtom,nextCoefs);
    //std::cout << nextAtom << std::endl;
    //std::cout << nextAtom.weight() << std::endl;
    //
    // create curr dictionary with old atoms plus new one
    // idem with coefs
    //
    // this is painful
    nextD.destroy();
    nextD.allocate(K+1,M);
    nextD.set_submatrix(0,0,currD);
    nextD.set_submatrix(K,0,nextAtom);
    currD.destroy();
    currD.allocate(K+1,M);
    nextD.copy_to(currD);
    nextD.destroy();

    nextA.destroy();
    nextA.allocate(N,K+1);
    nextA.set_submatrix(0,0,currA);
    nextA.set_submatrix(0,K,nextCoefs);
    currA.destroy();
    currA.allocate(N,K+1);
    nextA.copy_to(currA);
    nextA.destroy();

    learn_model_inner(X,currE,currD,currA);
    currL = model_codelength(currE,currD,currA);
    if ((currL + dev) < bestL) {
      stuck = 0;
      bestL = currL;
      D.destroy();
      D.allocate(K+1,M);
      currD.copy_to(D);
      A.destroy();
      A.allocate(N,K+1);
      currA.copy_to(A);
      currE.copy_to(E);
      bestK = K+1;
    } else { 
      stuck++;
      allStuck++;
      sumStuck += (currL-bestL);
      if (stuck >= 10) {
	std::cout << "No further improvement." << std::endl;
	break;
      }
    }
    K++;
    } while (stuck < 10);
  currD.destroy();
  currA.destroy();
  currE.destroy();
  nextE.destroy();
  nextAtom.destroy();
  nextCoefs.destroy();
  return bestL;
}

idx_t learn_model_mdl_backward_selection(binary_matrix& X,
					 binary_matrix& E, 
					 binary_matrix& D, 
					 binary_matrix& A) {
  const idx_t M = E.get_cols();
  const idx_t N = E.get_rows();
  idx_t K = D.get_rows();
  learn_model_inner(X,E,D,A);
  //mul(A,false,D,false,E);
  //add(E,X,E);
  idx_t bestL = model_codelength(E,D,A);
  idx_t bestK = K;
  idx_t currL = bestL;
  binary_matrix Dk(1,M);
  binary_matrix Ak(1,N);
  binary_matrix nextD,currD(D);
  binary_matrix nextA,currA(A);
  binary_matrix nextE(N,M);
  binary_matrix AkDk(N,M);
  
  idx_t stuck = 0;
  idx_t sumStuck = 0;
  idx_t allStuck = 0;
  for (; K > 0; K--) {
    int dif = int(currL) - int(bestL);
    int dev = allStuck > 0 ? (sumStuck/allStuck) : 0;
    std::cout << "currK=" << K << " currL=" << currL << " bestK=" << bestK << " bestL=" << bestL << " stuck=" << stuck << " dif=" << dif << " dev=" << dev << std::endl;
    idx_t nextk = 0;
    idx_t nextL = ~(1UL<<(sizeof(idx_t)-1));
    for (idx_t k = 0; k < K; k++) {
      currD.copy_row_to(k,Dk);
      currA.copy_col_to(k,Ak);
      mul(Ak,true,Dk,false,AkDk);
      // codelength of nex
      add(AkDk,E,nextE); // residual associated to removing Dk
      idx_t tmpL = model_codelength(nextE,currD,currA);
      tmpL -= universal_codelength(M,Dk.weight());
      tmpL -= universal_codelength(N,Ak.weight());
      if (tmpL < nextL) {
	nextL = tmpL;
	nextk = k;
      }
    }
    //
    // create new dictionary and coefficients without discarded atom bestk
    //    
    nextD.destroy();
    nextA.destroy();
    if (K > 1) {
      nextD.allocate(K-1,M);
      nextA.allocate(N,K-1);
      for (idx_t k = 0; k < nextk; k++) {
	currD.copy_row_to(k,Dk);
	nextD.set_row(k,Dk);
	currA.copy_col_to(k,Ak);
	nextA.set_col(k,Ak);
      }
      for (idx_t k = nextk+1; k < K; k++) {
	currD.copy_row_to(k,Dk);
	nextD.set_row(k-1,Dk);
	currA.copy_col_to(k,Ak);
	nextA.set_col(k-1,Ak);
      }
      learn_model_inner(X,nextE,nextD,nextA);
      nextL = model_codelength(nextE,nextD,nextA);
    } else {
      nextL = model_codelength(nextE,nextD,nextA);
    }

    if (nextL + dev < bestL)  {
      // keep best solution so far
      if (K == 1) {
	std::cout << "Resulted in empty model!" << std::endl;
	D.destroy();
	A.destroy();
	X.copy_to(E);
	break;
      }
      stuck = 0;
      bestK = (K-1);
      bestL = nextL;
      D.destroy();
      D.allocate(K-1,M);
      nextD.copy_to(D);
      A.destroy();
      A.allocate(N,K-1);
      nextA.copy_to(A);
      nextE.copy_to(E);
    } else { // not better, do not update solution
      stuck++;
      allStuck++;
      sumStuck += (nextL-bestL);
      if (stuck >= 10) {
	std::cout << "No further improvement." << std::endl;
	break;
      }
    }
    currD.destroy();
    currD.allocate(K-1,M);
    nextD.copy_to(currD);
    
    currA.destroy();
    currA.allocate(N,K-1);
    nextA.copy_to(currA);

    currL = nextL;
  }

  currD.destroy(); nextD.destroy();
  currA.destroy(); nextA.destroy();
  nextE.destroy();
  AkDk.destroy();
  Ak.destroy();
  Dk.destroy();
  return bestL;
}

idx_t learn_model_mdl_full_search(binary_matrix& X,
				  binary_matrix& E, 
				  binary_matrix& D, 
				  binary_matrix& A) {
  const idx_t M = E.get_cols();
  const idx_t N = E.get_rows();
  idx_t K = D.get_rows();
  
  binary_matrix  candE(N,M);
  idx_t bestL = 1UL<<30; // will not work on 16 bit machines
  idx_t bestk = 0;
  for (idx_t k = 20; k <= K; k+=20) {
  //for (idx_t k = 1; k < K; k++) {
    binary_matrix candD(k,M),candA(N,k);
    initialize_dictionary(X,candD,candA);
    learn_model_inner(X,candE,candD,candA);
#define REPS 10
#if 1
    idx_t aux[REPS];
    std::cout << "K=" << k;
    for (idx_t I = 0; I < REPS; I++) {
      random_seed = (random_seed*31 ) % 17;
      initialize_dictionary(X,candD,candA);
      learn_model_inner(X,candE,candD,candA);
      aux[I] = model_codelength(candE,candD,candA);
      std::cout << " " << aux[I];
    }
    idx_t candL = *std::min_element(aux+0,aux+REPS);
    std::cout << " " << candL << std::endl;
#else
    idx_t candL = model_codelength(candE,candD,candA);
    std::cout << "K=" << k << " candL=" << candL << std::endl;
#endif
    if (candL < bestL) {
      bestL = candL;
      bestk = k;
      candE.copy_to(E);

      D.destroy();
      D.allocate(k,M);
      candD.copy_to(D);

      A.destroy();
      A.allocate(N,k);
      candA.copy_to(A);
    }
    candD.destroy();
    candA.destroy();
  }
  candE.destroy();
  std::cout << "bestL=" << bestL << " bestk=" << bestk << std::endl;
  return bestL;
}
