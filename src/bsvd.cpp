#include <omp.h>
#include "util.h"
#include "bsvd.h"
#include "random_number_generation.h"
#include "update_dictionary.h"
#include "initialize_dictionary.h"
#include "encode_samples.h"
#include "entropy_coding.h"
#include <iomanip>
#include <cstdlib>
#include <algorithm>
#include "config.h"

idx_t learn_model_traditional(binary_matrix& X,
			      const binary_matrix& H,
			      binary_matrix& E,
			      binary_matrix& D, 
			      binary_matrix& A) {
  const idx_t K = D.get_rows();
  const idx_t ma  = K;
  const idx_t me  = 0;
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
    idx_t changed_coefs = encode_samples(E,H,D,A,ma,me);
  //  std::cout << "iter=" << std::setw(8) << iter 
//	      << "\t||E||=" << std::setw(8) << E.weight()
//	      << "\t||D||=" << std::setw(8) << D.weight()
//	      << "\t||A||=" << std::setw(8) << A.weight()
//	      << "\tchanged coefs=" << std::setw(8) << changed_coefs << std::endl;
    changed = changed_coefs + update_dictionary(E,H,D,A);
//    std::cout << "iter=" << std::setw(8) << iter 
//	      << "\t||E||=" << std::setw(8) << E.weight()
//	      << "\t||D||=" << std::setw(8) << D.weight()
//	      << "\t||A||=" << std::setw(8) << A.weight()
//	      << "\tchanged atoms=" << std::setw(8) << changed << std::endl;
  }
  return iter;
}


idx_t learn_model_alter1(binary_matrix& X,
			 const binary_matrix& H,
			 binary_matrix& E, 
			 binary_matrix& D, 
			 binary_matrix& A) {
  //  std::cout << "alter1" << std::endl;
  const idx_t N = E.get_rows();
  const idx_t M = E.get_cols();
  const idx_t K = D.get_rows();
  const idx_t ma  = K;
  const idx_t me  = 0;

  mul(A,false,D,false,E);
  add(E,X,E);
  binary_matrix Dt(M,K);
  binary_matrix At(K,N);
  binary_matrix Et(M,N);
  binary_matrix Ht(M,N);
  H.transpose_to(Ht);
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
    idx_t changed_coefs = encode_samples(E,H,D,A,ma,me);
//    std::cout << "DIRECT: iter=" << std::setw(8) << iter 
//	      << "\t||E||=" << std::setw(8) << E.weight()
//	      << "\t||D||=" << std::setw(8) << D.weight()
//	      << "\t||A||=" << std::setw(8) << A.weight()
//	      << "\tchanged coefs=" << std::setw(8) << changed_coefs << std::endl;
    changed = changed_coefs + update_dictionary(E,H,D,A);
  //  std::cout << "DIRECT: iter=" << std::setw(8) << iter 
//	      << "\t||E||=" << std::setw(8) << E.weight()
//	      << "\t||D||=" << std::setw(8) << D.weight()
//	      << "\t||A||=" << std::setw(8) << A.weight()
//	      << "\tchanged atoms=" << std::setw(8) << changed << std::endl;

    A.transpose_to(At);
    D.transpose_to(Dt);
    E.transpose_to(Et);

    changed_coefs = encode_samples(Et,Ht,At,Dt,ma,me);
  //  std::cout << "TRANSP: iter=" << std::setw(8) << iter 
//	      << "\t||E||=" << std::setw(8) << Et.weight()
//	      << "\t||D||=" << std::setw(8) << Dt.weight()
//	      << "\t||A||=" << std::setw(8) << At.weight()
//	      << "\tchanged coefs=" << std::setw(8) << changed_coefs << std::endl;

    changed = update_dictionary(Et,Ht,At,Dt);
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
			 const binary_matrix& H,
			 binary_matrix& E, 
			 binary_matrix& D, 
			 binary_matrix& A) {
  //  std::cout << "alter2" << std::endl;
  const idx_t N = E.get_rows();
  const idx_t M = E.get_cols();
  const idx_t K = D.get_rows();
  const idx_t ma  = K;
  const idx_t me  = 0;

  mul(A,false,D,false,E);
  add(E,X,E);
  binary_matrix Dt(M,K);
  binary_matrix At(K,N);
  binary_matrix Et(M,N);
  binary_matrix Ht(M,N);
  if (!H.empty()) H.transpose_to(Ht);
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
      idx_t changed_coefs = encode_samples(E,H,D,A,ma,me);
//      std::cout << "DIRECT: iter=" << std::setw(8) << iter 
//		<< "\t||E||=" << std::setw(8) << E.weight()
//		<< "\t||D||=" << std::setw(8) << D.weight()
//		<< "\t||A||=" << std::setw(8) << A.weight()
//		<< "\tchanged coefs=" << std::setw(8) << changed_coefs << std::endl;
      changed = changed_coefs + update_dictionary(E,H,D,A);
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
      idx_t changed_coefs = encode_samples(Et,Ht,At,Dt,ma,me);
//      std::cout << "TRANSPOSED: iter=" << std::setw(8) << iter 
//		<< "\t||E||=" << std::setw(8) << Et.weight()
//		<< "\t||D||=" << std::setw(8) << Dt.weight()
//		<< "\t||A||=" << std::setw(8) << At.weight()
//		<< "\tchanged coefs=" << std::setw(8) << changed_coefs << std::endl;

      changed = changed_coefs + update_dictionary(Et,Ht,At,Dt);
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
			 const binary_matrix& H,
			 binary_matrix& E, 
			 binary_matrix& D, 
			 binary_matrix& A) {
  //  std::cout << "alter3" << std::endl;
  const idx_t N = E.get_rows();
  const idx_t M = E.get_cols();
  const idx_t K = D.get_rows();
  const idx_t ma = K;
  const idx_t me  = 0;
  mul(A,false,D,false,E);
  add(E,X,E);

  binary_matrix Dt(M,K);
  binary_matrix At(K,N);
  binary_matrix Et(M,N);
  binary_matrix Ht(M,N);
  if (!H.empty()) H.transpose_to(Ht);
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
    changed = update_dictionary(Et,Ht,At,Dt);
//    std::cout << "iter=" << std::setw(8) << iter 
//	      << "\t||E||=" << std::setw(8) << Et.weight()
//	      << "\t||D||=" << std::setw(8) << Dt.weight()
//	      << "\t||A||=" << std::setw(8) << At.weight()
//	      << "\tchanged atoms=" << std::setw(8) << changed << std::endl;

    At.transpose_to(A);
    Dt.transpose_to(D);
    Et.transpose_to(E);
    changed = update_dictionary(E,H,D,A);
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

#include "entropy_coding.h"

idx_t learn_model_mdl_forward_selection(binary_matrix& X,
					const binary_matrix& H,
					binary_matrix& E, 
					binary_matrix& D, 
					binary_matrix& A) {
  const idx_t M = E.get_cols();
  const idx_t N = E.get_rows();
  idx_t K = D.get_rows();
  learn_model_inner(X,H,E,D,A);
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
    initialize_dictionary(currE,H,nextAtom,nextCoefs);
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

    learn_model_inner(X,H,currE,currD,currA);
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
					 const binary_matrix& H,
					 binary_matrix& E, 
					 binary_matrix& D, 
					 binary_matrix& A) {
  const idx_t M = E.get_cols();
  const idx_t N = E.get_rows();
  idx_t K = D.get_rows();
  learn_model_inner(X,H,E,D,A);
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
      learn_model_inner(X,H,nextE,nextD,nextA);
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
				  const binary_matrix& H,
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
    initialize_dictionary(X,H,candD,candA);
    learn_model_inner(X,H,candE,candD,candA);
#define REPS 10
#if 1
    idx_t aux[REPS];
    std::cout << "K=" << k;
    for (idx_t I = 0; I < REPS; I++) {
      random_seed = (random_seed*31 ) % 17;
      initialize_dictionary(X,H,candD,candA);
      learn_model_inner(X,H,candE,candD,candA);
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
