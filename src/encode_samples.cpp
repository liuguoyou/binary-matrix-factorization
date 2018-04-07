#include <algorithm>
#include <omp.h>
#include <iomanip>
#include "encode_samples.h"
#include "intmat.h"
#include <cassert>

//==========================================================================
//
// AUXILIARY FUNCTIONS FOR CORRELATION-BASED ENCODING
//

/// compute unnormalized  correlation between dict. D and sample x, D (pxm), x (1xm): g = Dx^t and
/// squared norm of the rows of D, w
static void  compute_correlation(const binary_matrix& D, const binary_matrix& x, integer_matrix& g, integer_matrix& w) {
  const size_t p = D.get_rows();
  const size_t m = D.get_cols();
  assert(x.get_cols() == m);
  assert(g.get_cols() == p);
  assert(w.get_cols() == p);
  mul_ABt(D,x,g); // integer
  // now square every element in g
  // and compute the row weights of D
  for (size_t i = 0; i < p; i++) {
    w.set(0,i,D.row_weight(i));
  }
}


/// update correlation given a 
static void  update_correlation(const binary_matrix& D, const idx_t& i, integer_matrix& g) {
  const size_t p = D.get_rows();
  const size_t m = D.get_cols();
  assert(g.get_cols() == p);
  //
  // this is tricky: we want to update the standard correlation between two Euclidean vectors
  // but the newly added atom is combined using XOR addition, so the result is more convolved
  //
  // Typically, if we have a correlation at step k g^k = D^tr^k and a new atom D_i is added (with coefficient a_i^{k+1}
  // we would have g^{k+1} = D^t(r^k - a_i^{k+1}D_i)= D^tr^k - a_i^{k+1}DD_i^t = Dr^k - a_i^{k+1}G_i
  //
  // Here the thing is very different: a_i doesn't matter, it only matters that it switched, so only the index i
  // matters. We have
  // g^{k+1} = D( r^k XOR {k+1}D_i )
  // but g^{k+1} is still the traditional (unnormalized) correlation, so we have to update it using the result
  // of the XOR operation for each element j in r^k and D_i:
  // furthermore, we need to take into account the norm of the atoms, but using integer only arithmethic, so
  // we cannot use the Euclidean norm, but the squared euclidean norm
  // theremore we must work with squared correlations: g^k_j = (D_j*r^k)^2 / h(D_j) = 
  //       r^k_j  D_ij    r^{k+1}_j
  // j= 0  0      0       0         nothing changes
  // j= 1  0      1       1         the j-th row  of D is ADDED (in integer arithmethic) to the unnormalized g^k
  //       1      0       1         nothing changes
  // j= K  1      1       0         the j-th row of D is SUBSTRACTED from  g^k
  //
  // NOTE: the operations below are transposed to what the text says, so rows are columns, columns  are rows, etc.
  binary_matrix D_i = D.get_row(i);
  binary_matrix Dcol(1,p);
  for (size_t j = 0; j < m; j++) {
    if (D.get(0,j)) {
      D.copy_col_to(j,Dcol);
      if (g.get(0,j)) { // substract column
	for (size_t k = 0; k < m; k++) {
	  if (Dcol.get(0,k))
	    g.dec(0,k);
	}
      } else { // add column to g
	for (size_t k = 0; k < m; k++) {
	  if (Dcol.get(0,k))
	    g.inc(0,k);
	}
      }
    }
  }
  D_i.destroy();
  Dcol.destroy();
}

idx_t get_most_correlated(const integer_matrix& g, const integer_matrix& w) {
  const size_t m = g.get_rows();
  assert(w.get_rows() == m);
  double gk = g.get(0,0);
  idx_t max_k = 0;
  double max_corr = gk*gk / (double) w.get(0,0);
  for (size_t k = 1; k < m; k++) {
    gk = g.get(0,k);
    double corr = gk*gk / w.get(0,k);
    if (corr > max_corr) {
      max_k = k;
      max_corr = corr;
    }
  }
  return max_k;
}
		  
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
  binary_matrix Ei(1,m);
  binary_matrix Ai(1,p);
  for (idx_t i = 0; i < n; i++) {
    E.copy_row_to(i,Ei);
    A.copy_row_to(i,Ai);
    bool improved = true;
    bool ichanged = false;
    idx_t iter = 0;
    while (improved) {
      // PENDING!!!
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
