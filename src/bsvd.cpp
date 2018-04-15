#include <omp.h>
#include "util.h"
#include "bsvd.h"
#include "random_number_generation.h"
#include "update_dictionary.h"
#include "initialize_dictionary.h"
#include "coefficients_update.h"
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
  const idx_t me  = get_max_err_weight();
  mul(A,false,D,false,E);
  add(E,X,E);
  idx_t changed = 1;
  idx_t iter = 0;
  if (get_verbosity() >= 2) {
    std::cout << "trad" << std::endl;
    std::cout << "iter=" << std::setw(8) << iter 
	    << "\t||E||=" << std::setw(8) << E.weight()
	    << "\t||D||=" << std::setw(8) << D.weight()
	    << "\t||A||=" << std::setw(8) << A.weight() << std::endl;
  }
  while (changed > 0) {    
    iter++;
    idx_t changed_coefs = coefficients_update(E,H,D,A,ma,me);
    if (get_verbosity() >= 2) {
      std::cout << "iter=" << std::setw(8) << iter 
		<< "\t||E||=" << std::setw(8) << E.weight()
		<< "\t||D||=" << std::setw(8) << D.weight()
		<< "\t||A||=" << std::setw(8) << A.weight()
		<< "\tchanged coefs=" << std::setw(8) << changed_coefs << std::endl;
    }
    changed = changed_coefs + update_dictionary(E,H,D,A);
    if (get_verbosity() >= 2) {
      codelength L = model_codelength(E,D,A);
      std::cout << "iter=" << std::setw(8) << iter 
		<< "\t||E||=" << std::setw(8) << E.weight()
		<< "\t||D||=" << std::setw(8) << D.weight()
		<< "\t||A||=" << std::setw(8) << A.weight()
		<< "\tL(E)=" << std::setw(8) << L.E
		<< "\tL(D)=" << std::setw(8) << L.D
		<< "\tL(A)=" << std::setw(8) << L.A
		<< "\tL(X)=" << std::setw(8) << L.X
		<< "\tchanged atoms=" << std::setw(8) << changed << std::endl;
    }
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
    idx_t changed_coefs = coefficients_update(E,H,D,A,ma,me);
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

    changed_coefs = coefficients_update(Et,Ht,At,Dt,ma,me);
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
      idx_t changed_coefs = coefficients_update(E,H,D,A,ma,me);
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
      idx_t changed_coefs = coefficients_update(Et,Ht,At,Dt,ma,me);
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
  codelength Lparts = model_codelength(E,D,A); 
  idx_t bestL = Lparts.X; 
  idx_t stuck = 0;
  idx_t sumStuck = 0;
  idx_t allStuck = 0;
  //  Lparts = model_codelength(currE,currD,currA);
  do {
    //
    // initialize curr atom and associated coefs.
    //
    idx_t currL = Lparts.X;
    int dif = int(currL) - int(bestL);
    int dev = allStuck > 0 ? (sumStuck/allStuck) : 0;    
    std::cout << "currK=" << K
	      << " |currE|=" << currE.weight()
	      << " |currA|=" << currA.weight()
	      << " |currD|=" << currD.weight()
	      << " currL=" << currL << " bestK=" << bestK
	      << " bestL=" << bestL << " stuck=" << stuck
	      << " dif=" << dif << " dev=" << dev << std::endl;
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
    Lparts = model_codelength(currE,currD,currA); // RARO: dos veces??
    currL = Lparts.X;
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

#if 0
//
// THIS IMPLEMENTATION IS FASTER BUT IS NOT WORKING WELL
//
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
  codelength Lparts = model_codelength(E,D,A);
  idx_t bestL = Lparts.X;
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
    // need a good rationale for this!
    int dev = allStuck > 0 ? (sumStuck/allStuck) : 0;
    std::cout << "currK=" << K << " currL=" << currL << " bestK=" << bestK << " bestL=" << bestL << " stuck=" << stuck << " dif=" << dif << " dev=" << dev << std::endl;
    idx_t nextk = 0;
    idx_t nextL = ~(1UL<<(sizeof(idx_t)-1));
    for (idx_t r = 0; r < K; k++) {
      currD.copy_row_to(r,Dk);
      currA.copy_col_to(r,Ak);
      mul(Ak,true,Dk,false,AkDk);
      // codelength of nex
      add(AkDk,E,nextE); // residual associated to removing Dk
      const codelength Lparts = model_codelength(nextE,currD,currA);
      const idx_t LDk = universal_codelength(M,Dk.weight());
      const idx_t LAk =  universal_codelength(N,Ak.weight());
      const idx_t tmpL = Lparts.X - LDk - LAk;
      if (get_verbosity() >= 1) {
	std::cout << "backsel: K=" << K << " candidate= " << k
		  << " L(E)=" << Lparts.E
		  << " L(D)=" << Lparts.D
		  << " L(A)=" << Lparts.A
		  << " -L(Dk)=" << LDk
		  << " -L(Ak)=" << LAk
		  << " L=" << tmpL << "  nextk=" << nextk << " nextL=" << nextL << std::endl;
      }
      if (tmpL < nextL) {
	nextL = tmpL;
	nextk = k;
	if (get_verbosity() >= 1) {
	  std::cout << "new best! " << std::endl;
	}
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
      //learn_model_inner(X,H,nextE,nextD,nextA);
      const codelength Lparts = model_codelength(nextE,nextD,nextA);
      nextL = Lparts.X;
    } else {
      const codelength Lparts = model_codelength(nextE,nextD,nextA);
      nextL = Lparts.X;
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

#else

//
// SLOWER IMPLEMENTATION, BUT SAFE
//
idx_t learn_model_mdl_backward_selection(binary_matrix& X,
					 const binary_matrix& H,
					 binary_matrix& E, 
					 binary_matrix& D, 
					 binary_matrix& A) {
  const idx_t M = E.get_cols();
  const idx_t N = E.get_rows();
  idx_t K = D.get_rows();
  learn_model_inner(X,H,E,D,A);
  binary_matrix Dk(1,M);
  binary_matrix Ak(1,N);
  binary_matrix bestD,candD,innerD;
  binary_matrix bestA,candA,innerA;
  binary_matrix innerE(N,M);
  //
  // first stage: get rid of unused atoms
  //
  idx_t nused = 0;
  char iused[K];
  for (size_t k = 0; k < K; k++) {
    if (A.col_weight(k) > 0) {
      iused[k] = 1;
      nused++;
    } else {
      iused[k] = 0;
    }
  }
  if (get_verbosity() >= 1) {
    std::cout << "Removing " << (K-nused) << " unused atoms." << std::endl;
  }
  codelength Lparts = model_codelength(E,D,A);
  idx_t bestL = Lparts.X;
  idx_t bestK = nused;
  idx_t candL = bestL;
  idx_t candK = bestK;
  //
  // initial candidate is best candidate after pruning
  //
  bestD.allocate(bestK,M);
  bestA.allocate(N,bestK);
  candD.allocate(bestK,M);
  candA.allocate(N,bestK);
  for (size_t k = 0, r = 0; r < bestK; k++) {
    if (iused[k]) {
      D.copy_row_to(k,Dk);
      bestD.set_row(r,Dk);
      candD.set_row(r,Dk);
      A.copy_col_to(k,Ak);
      bestA.set_row(r,Ak);
      candA.set_row(r,Ak);
      r++;
    }
  }
  idx_t stuck = 0;
  idx_t sumStuck = 0;
  idx_t allStuck = 0;
  for (candK = bestK-1; candK > 1; candK--) {
    int dif = int(candL) - int(bestL);
    // need a good rationale for this!
    int dev = allStuck > 0 ? (sumStuck/allStuck) : 0;
    std::cout << "candK=" << candK << " candL=" << candL << " bestK=" << bestK << " bestL=" << bestL << " stuck=" << stuck << " dif=" << dif << " dev=" << dev << std::endl;
    idx_t best_inner_k = 0;
    idx_t best_inner_L = ~(1UL<<(sizeof(idx_t)-1));

    innerD.destroy();
    innerA.destroy();
    innerD.allocate(candK,M);
    innerA.allocate(N,candK);
    for (idx_t k = 0; k < candK; k++) {
      //
      // copy all but the k-th row and col of D and A
      //
      for (idx_t r = 0; r < k; r++) {
	candD.copy_row_to(r,Dk);
	innerD.set_row(r,Dk);
	candA.copy_col_to(r,Ak);
	innerA.set_col(r,Ak);
      }
      for (idx_t r = k+1; r < (candK+1); r++) {
	candD.copy_row_to(r,Dk);
	innerD.set_row(r-1,Dk);
	candA.copy_col_to(r,Ak);
	innerA.set_col(r-1,Ak);
      }
      //
      // update model after removing k-th atom
      //
      learn_model_inner(X,H,innerE,innerD,innerA);
      //
      // compute codelength
      // 
      const codelength Lparts = model_codelength(innerE,innerD,innerA);
      if (get_verbosity() >= 1) {
	std::cout << "backsel: candK=" << candK << " candidate= " << k
		  << " L(E)=" << Lparts.E
		  << " L(D)=" << Lparts.D
		  << " L(A)=" << Lparts.A
		  << " L=" << Lparts.X << "  bestk=" << best_inner_k << " bestL=" << best_inner_L << std::endl;
      }
      //
      // compare to best within the current inner iteration (over k)
      //
      if (Lparts.X < best_inner_L) {
	best_inner_L =Lparts.X;
	best_inner_k = k;
      }
    }
    //
    // see if we have improved compared to previous dictionary size 
    //
    if (best_inner_L + dev < bestL) {
      std::cout << "new best K " << candK << std::endl;
      //
      // YES: best <- candidate - k-th atom
      // 
      bestL = best_inner_L;
      bestK = candK;
      bestD.destroy();
      bestD.allocate(candK,M);      
      bestA.destroy();
      bestA.allocate(N,candK);
      for (idx_t r = 0; r < best_inner_k; r++) {
	candD.copy_row_to(r,Dk);
	bestD.set_row(r,Dk);
	candA.copy_col_to(r,Ak);
	bestA.set_col(r,Ak);
      }
      for (idx_t r = best_inner_k+1; r < candK; r++) {
	candD.copy_row_to(r,Dk);
	bestD.set_row(r-1,Dk);
	candA.copy_col_to(r,Ak);
	bestA.set_col(r-1,Ak);
      }
    } else {
      //
      // NO improvement
      //
      stuck++;
      allStuck++;
      sumStuck += (candL-bestL);
      if (stuck >= 10) {
	std::cout << "No further improvement." << std::endl;
	break;
      }
    }
    //
    // In any case, if we are still searching, new candidate 
    // is the best dictionary obtained in this iteration
    // candidate <- candidate - k+th atom
    for (idx_t r = 0; r < best_inner_k; r++) {
      candD.copy_row_to(r,Dk);
      innerD.set_row(r,Dk);
      candA.copy_col_to(r,Ak);
      innerA.set_col(r,Ak);
    }
    for (idx_t r = best_inner_k+1; r < candK; r++) {
      candD.copy_row_to(r,Dk);
      innerD.set_row(r-1,Dk);
      candA.copy_col_to(r,Ak);
      innerA.set_col(r-1,Ak);
    }    
    candD.destroy();
    candA.destroy();
    candD.allocate(candK,M);
    candA.allocate(N,candK);
    candL = best_inner_L;
  }

  candD.destroy(); innerD.destroy();
  candA.destroy(); innerA.destroy();
  innerE.destroy();
  Ak.destroy();
  Dk.destroy();
  if (bestK < D.get_rows()) {
    D.destroy();
    A.destroy();
    D.allocate(bestK,M);
    A.allocate(N,bestK);
    bestD.copy_to(D);
    bestA.copy_to(A);    
  }
  return bestL;
}


#endif

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
#define REPS 1
#if 1
    idx_t aux[REPS];
    std::cout << "K=" << k;
    for (idx_t I = 0; I < REPS; I++) {
      random_seed = (random_seed*31 ) % 17;
      initialize_dictionary(X,H,candD,candA);
      learn_model_inner(X,H,candE,candD,candA);
      const codelength Lparts = model_codelength(candE,candD,candA);
      aux[I] = Lparts.X;
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
