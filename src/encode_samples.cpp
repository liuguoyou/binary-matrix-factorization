#include <algorithm>
#include <omp.h>
#include <iomanip>
#include "encode_samples.h"
#include "intmat.h"
#include <cassert>

//==========================================================================

idx_t encode_samples_corr(binary_matrix& E,
			  const binary_matrix& H,
			  const binary_matrix& D,
			  binary_matrix& A,
			  const idx_t max_a_weight,
			  const idx_t max_e_weight) 
{
  const idx_t m = E.get_cols();
  const idx_t n = E.get_rows();
  const idx_t p = D.get_rows();
  std::cout << "cu/corr" << std::endl;
  //
  // SPARSE CODING STEP
  // Go over each row and encode it using rows from D until the residual weight
  // is no longer diminished.
  //
  idx_t changed = 0;
  //
  // weight of each row of Dw, for normalized correlation
  //
  integer_matrix Dw(1,p);
  for (size_t i = 0; i < p; i++) {
    Dw.set(0,i,D.row_weight(i)); // weight is SQUARED Euclidean norm
  }
  //
  // modulo-2 Gram matrix of D
  //
  binary_matrix G2(p,p);
  mul_ABt(D,D,G2);
  binary_matrix Ei(1,m);
  binary_matrix Ai(1,p);
  binary_matrix Dk(1,m);
  integer_matrix gi(1,p);
  for (idx_t i = 0; i < n; i++) {
    E.copy_row_to(i,Ei);
    A.clear();
    // compute initial (unnormalized) correlation with D and Ei
    mul_ABt(D,Ei,gi); 
    bool improved = false;
    bool ichanged = false;
    idx_t t = 0;
    double max_corr, max_weight;
    idx_t max_k;
    while (improved) {
      //
      // current weight of residual
      // 
      const idx_t ew = Ei.weight();
      if (ew == 0) { continue; };
      //
      // get atom with maximum correlation 
      //
      double gk = gi.get(0,0);
      max_k = 0;
      max_corr = gk*gk;
      max_weight = (double) Dw.get(0,0);
      for (size_t k = 1; k < m; k++) {
	gk = gi.get(0,k);
	// comparison of correlation is squared
	const double corr = gk*gk;
	const double weight = Dw.get(0,k);
	// instead of corr/weight > max_corr/max_weight
	if (max_weight*corr > weight*max_corr) { 
	  max_k = k;
	  max_corr = corr;
	  max_weight = weight;
	}
      }
      if (max_corr == 0) {
	improved = false; break;
      }
      // if there is correlation, go on
      Dk = D.get_row(max_k);
      const idx_t wk = Dw.get(0,max_k);
      //
      // add Dk to Ei mod 2
      //
      add(Ei,Dk,Ei);
      //
      // update correlation: g(t+1) = g(t) - G2_k
      //
      for (size_t i = 0; i < p; i++) {
	if (G2.get(max_k,i)) {
	  gi.dec(0,i);
	}
      }
      // next iteration t <- t+1
      t++;
    } // while there is any improvement in the i-th sample residual
    if (ichanged) {
      changed++;
      E.set_row(i,Ei);
      A.set_row(i,Ai);
    }
    A.set_row(i,Ai);
    std::cout << "i=" << i << " changed=" << ichanged << " |Ei|=" << Ei.weight() << "|Ai|=" << Ai.weight() << std::endl;
  } // for each reow in E
  gi.destroy();
  Ei.destroy();
  Ai.destroy();
  return changed;
}

//==========================================================================

idx_t encode_samples_basic(binary_matrix& E,
			   const binary_matrix& H,
			   const binary_matrix& D,
			   binary_matrix& A,
			   const idx_t max_a_weight,
			   const idx_t max_e_weight) 
{
  const idx_t m = E.get_cols();
  const idx_t n = E.get_rows();
  const idx_t p = D.get_rows();
  std::cout << "cu/basic" << std::endl;
  //
  // SPARSE CODING STEP
  // Go over each row and encode it using rows from D until the residual weight
  // is no longer diminished.
  //
  idx_t changed = 0;
  binary_matrix Ei(1,m);
  binary_matrix Ai(1,p);
  binary_matrix Dk(1,m);
  for (idx_t i = 0; i < n; i++) {
    E.copy_row_to(i,Ei);
    A.copy_row_to(i,Ai);
    bool improved = true;
    bool ichanged = false;
    idx_t iter = 0;
    while (improved) {
      idx_t w = Ei.weight();
      if (w <= max_e_weight) { // reached maximum error goal
	break;
      }
      if (Ai.weight() >= max_a_weight) { // reached maximum allowable weight in A
	break;
      }
      D.copy_row_to(0,Dk);
      idx_t bestk = 0;
      idx_t bestd = dist(Ei,Dk);
      //bool_and(Ei,Dk,Dk);
      //     idx_t r = Dk.weight();
      //      std::cout << "\titer=" << iter << " k=" << 0 << " r=" << r << " d=" << bestd << std::endl;
      for (idx_t k = 1; k < p; k++) {
	D.copy_row_to(k,Dk);
	const idx_t dk = dist(Ei,Dk);
	//bool_and(Ei,Dk,Dk);
	//	idx_t r = Dk.weight();
	//	std::cout << "\titer" << iter << " k=" << k << " r=" << r << " d=" << dk << std::endl;
	if (dk < bestd) {
	  bestd = dk;
	  bestk = k;
	} 
      } // for each candidate atom
      //      std::cout << "i=" << i << " w=" << w << " bestk=" << bestk << " bestd=" << bestd << std::endl;
      if (bestd < w) {
	D.copy_row_to(bestk,Dk);
	Ai.flip(0,bestk);
	add(Ei,Dk,Ei);
	w = bestd;
	ichanged = true;
      } else {
	improved = false;
      }  
      iter++;
    } // while there is any improvement in the i-th sample residual
    if (ichanged) {
      changed++;
      E.set_row(i,Ei);
      A.set_row(i,Ai);
    }
    std::cout << "i=" << i << " changed=" << ichanged << " |Ei|=" << Ei.weight() << "|Ai|=" << Ai.weight() << std::endl;
  } // for each row in E
  Ei.destroy();
  Ai.destroy();
  Dk.destroy();
  return changed;
}

//==========================================================================

idx_t encode_samples_missing_data_omp(binary_matrix& E,
				      const binary_matrix& H,
				      const binary_matrix& D,
				      binary_matrix& A,
				      const idx_t max_a_weight,
				      const idx_t max_e_weight);

//==========================================================================

idx_t encode_samples_omp(binary_matrix& E,
			 const binary_matrix& H,
			 const binary_matrix& D,
			 binary_matrix& A,
			 const idx_t max_a_weight,
			 const idx_t max_e_weight) 
{
  if (!H.empty()) return encode_samples_missing_data_omp(E,H,D,A,max_a_weight,max_e_weight);

  //  std::cout << "cu/omp" << std::endl;
  const idx_t m = E.get_cols();
  const idx_t n = E.get_rows();
  const idx_t p = D.get_rows();
  //
  // SPARSE CODING STEP
  // Go over each row and encode it using rows from D until the residual weight
  // is no longer diminished.
  //
  idx_t changed = 0;
  omp_set_num_threads(omp_get_max_threads());
  idx_t NT;
#pragma omp parallel
  {  
 NT = omp_get_num_threads(); 
  }
  //  std::cout << "THREADS=" << NT << std::endl;
  binary_matrix Ei[NT];
  binary_matrix Ai[NT];
  binary_matrix Dk[NT];
  for (idx_t TT = 0; TT < NT; TT++) {
    Ei[TT].allocate(1,m);
    Ai[TT].allocate(1,p);
    Dk[TT].allocate(1,m);
  }

#pragma omp parallel for schedule(dynamic)
  for (idx_t i = 0; i < n; i++) {
    const idx_t T = omp_get_thread_num();
    E.copy_row_to(i,Ei[T]);
    A.copy_row_to(i,Ai[T]);
    bool improved = true;
    bool ichanged = false;
    idx_t iter = 0;
    while (improved) {
      idx_t w = Ei[T].weight();
      if (w <= max_e_weight) { // reached maximum error goal
#ifdef DEBUG
	std::cout << "STOP: error below maximum allowed error weight. " << std::endl;
#endif
	break;
      }
      if (Ai[T].weight() >= max_a_weight) { // reached maximum allowable weight in A
#ifdef DEBUG
	std::cout << "STOP: coefficients above maximum allowed  weight. " << std::endl;
#endif
	break;
      }
      D.copy_row_to(0,Dk[T]);
      idx_t bestk = 0, bestd = dist(Ei[T],Dk[T]);
      //bool_and(Ei[T],Dk[T],Dk[T]);
      //idx_t r = Dk[T].weight();q
      //std::cout << "\titer=" << iter << " k=" << 0 << " d=" << bestd << std::endl;
      for (idx_t k = 1; k < p; k++) {
	D.copy_row_to(k,Dk[T]);
	const idx_t dk = dist(Ei[T],Dk[T]);
	//bool_and(Ei[T],Dk[T],Dk[T]);
	//	idx_t r = Dk[T].weight();
	//std::cout << "\titer" << iter << " k=" << k << " d=" << dk << std::endl;
	if (dk < bestd) {
          bestd = dk;
          bestk = k;
	} 
      }
#ifdef DEBUG
      //std::cout << "i=" << std::setw(10) << i << " iter="<< std::setw(4) << iter;
      //std::cout << " wk=" << std::setw(3) << Dk[T].weight() << " w=" << w << " bestk=" << std::setw(4) << bestk << " bestd=" << std::setw(4) << bestd << std::endl;
#endif
      if (bestd < w) {
	D.copy_row_to(bestk,Dk[T]);
	Ai[T].flip(0,bestk);
	add(Ei[T],Dk[T],Ei[T]);
	w = bestd;
	ichanged = true;
      } else {
	improved = false;
      }  
      iter++;
    } // while there is any improvement
    if (ichanged) {
      changed++;
      E.set_row(i,Ei[T]);
      A.set_row(i,Ai[T]);
    }
#ifdef DEBUG
    std::cout << "i=" << std::setw(10) << i << " changed=" << std::setw(4) << ichanged << " |Ei|=" << std::setw(4) << Ei[T].weight() << " |Ai|=" << std::setw(4) << Ai[T].weight() << std::endl;
#endif
  }
  for (idx_t T = 0; T < NT; T++) {
    Ei[T].destroy();
    Ai[T].destroy();
    Dk[T].destroy();
  }
  return changed;
}

//==========================================================================

idx_t encode_samples_missing_data_omp(binary_matrix& E,
				      const binary_matrix& H,
				      const binary_matrix& D,
				      binary_matrix& A,
				      const idx_t max_a_weight,
				      const idx_t max_e_weight) 
{
  //  std::cout << "cu/omp" << std::endl;
  const idx_t m = E.get_cols();
  const idx_t n = E.get_rows();
  const idx_t p = D.get_rows();
  //
  // SPARSE CODING STEP
  // Go over each row and encode it using rows from D until the residual weight
  // is no longer diminished.
  //
  idx_t changed = 0;
  omp_set_num_threads(omp_get_max_threads());
  idx_t NT;
#pragma omp parallel
  {  
 NT = omp_get_num_threads(); 
  }
  //  std::cout << "THREADS=" << NT << std::endl;
  binary_matrix Ei[NT];
  binary_matrix Hi[NT];
  binary_matrix Ai[NT];
  binary_matrix Dk[NT];
  for (idx_t TT = 0; TT < NT; TT++) {
    Ei[TT].allocate(1,m);
    Hi[TT].allocate(1,m);
    Ai[TT].allocate(1,p);
    Dk[TT].allocate(1,m);
  }

#pragma omp parallel for schedule(dynamic)
  for (idx_t i = 0; i < n; i++) {
    const idx_t T = omp_get_thread_num();
    E.copy_row_to(i,Ei[T]);
    H.copy_row_to(i,Hi[T]);
    A.copy_row_to(i,Ai[T]);
    bool improved = true;
    bool ichanged = false;
    idx_t iter = 0;
    while (improved) {
      idx_t w = Ei[T].weight();
      if (w <= max_e_weight) { // reached maximum error goal
#ifdef DEBUG
	std::cout << "STOP: error below maximum allowed error weight. " << std::endl;
#endif
	break;
      }
      if (Ai[T].weight() >= max_a_weight) { // reached maximum allowable weight in A
#ifdef DEBUG
	std::cout << "STOP: coefficients above maximum allowed  weight. " << std::endl;
#endif
	break;
      }
      D.copy_row_to(0,Dk[T]);
      idx_t bestk = 0, bestd = weighted_dist(Ei[T],Dk[T],Hi[T]);
      //bool_and(Ei[T],Dk[T],Dk[T]);
      //idx_t r = Dk[T].weight();q
      //std::cout << "\titer=" << iter << " k=" << 0 << " d=" << bestd << std::endl;
      for (idx_t k = 1; k < p; k++) {
	D.copy_row_to(k,Dk[T]);
	const idx_t dk = weighted_dist(Ei[T],Dk[T],Hi[T]);
	//bool_and(Ei[T],Dk[T],Dk[T]);
	//	idx_t r = Dk[T].weight();
	//std::cout << "\titer" << iter << " k=" << k << " d=" << dk << std::endl;
	if (dk < bestd) {
          bestd = dk;
          bestk = k;
	} 
      }
#ifdef DEBUG
      //std::cout << "i=" << std::setw(10) << i << " iter="<< std::setw(4) << iter;
      //std::cout << " wk=" << std::setw(3) << Dk[T].weight() << " w=" << w << " bestk=" << std::setw(4) << bestk << " bestd=" << std::setw(4) << bestd << std::endl;
#endif
      if (bestd < w) {
	D.copy_row_to(bestk,Dk[T]);
	Ai[T].flip(0,bestk);
	add(Ei[T],Dk[T],Ei[T]);
	w = bestd;
	ichanged = true;
      } else {
	improved = false;
      }  
      iter++;
    } // while there is any improvement
    if (ichanged) {
      changed++;
      E.set_row(i,Ei[T]);
      A.set_row(i,Ai[T]);
    }
#ifdef DEBUG
    std::cout << "i=" << std::setw(10) << i << " changed=" << std::setw(4) << ichanged << " |Ei|=" << std::setw(4) << Ei[T].weight() << " |Ai|=" << std::setw(4) << Ai[T].weight() << std::endl;
#endif
  }
  for (idx_t T = 0; T < NT; T++) {
    Ei[T].destroy();
    Ai[T].destroy();
    Dk[T].destroy();
  }
  return changed;
}



//==========================================================================
