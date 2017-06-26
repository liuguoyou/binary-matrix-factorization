#include "encode_samples.h"
#include <algorithm>

idx_t encode_samples_basic(binary_matrix& E, const binary_matrix& D,  binary_matrix& A) 
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
      D.copy_row_to(0,Dk);
      idx_t bestk = 0, bestd = dist(Ei,Dk);
      bool_and(Ei,Dk,Dk);
      //     idx_t r = Dk.weight();
      //      std::cout << "\titer=" << iter << " k=" << 0 << " r=" << r << " d=" << bestd << std::endl;
      for (idx_t k = 1; k < p; k++) {
	D.copy_row_to(k,Dk);
	const idx_t dk = dist(Ei,Dk);
	bool_and(Ei,Dk,Dk);
	//	idx_t r = Dk.weight();
	//	std::cout << "\titer" << iter << " k=" << k << " r=" << r << " d=" << dk << std::endl;
	if (dk < bestd) {
	  bestd = dk;
	  bestk = k;
	} 
      }
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
    } // while there is any improvement
    if (ichanged) {
      changed++;
      E.set_row(i,Ei);
      A.set_row(i,Ai);
    }
  }
  Ei.destroy();
  Ai.destroy();
  Dk.destroy();
  return changed;
}



idx_t encode_samples_omp(binary_matrix& E, const binary_matrix& D,  binary_matrix& A) 
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
      D.copy_row_to(0,Dk[T]);
      idx_t bestk = 0, bestd = dist(Ei[T],Dk[T]);
      bool_and(Ei[T],Dk[T],Dk[T]);
      //     idx_t r = Dk[T].weight();
      //      std::cout << "\titer=" << iter << " k=" << 0 << " r=" << r << " d=" << bestd << std::endl;
      for (idx_t k = 1; k < p; k++) {
	D.copy_row_to(k,Dk[T]);
	const idx_t dk = dist(Ei[T],Dk[T]);
	bool_and(Ei[T],Dk[T],Dk[T]);
	//	idx_t r = Dk[T].weight();
	//	std::cout << "\titer" << iter << " k=" << k << " r=" << r << " d=" << dk << std::endl;
	if (dk < bestd) {
          bestd = dk;
          bestk = k;
	} 
      }
      //      std::cout << "i=" << i << " w=" << w << " bestk=" << bestk << " bestd=" << bestd << std::endl;
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
  }
  for (idx_t T = 0; T < NT; T++) {
    Ei[T].destroy();
    Ai[T].destroy();
    Dk[T].destroy();
  }
  return changed;
}


  //
  // DOES NOT WORK WELL!
  // I need to find proper surrogates to the distance between E and D that
  // can actually be stored.
  // But speed is of no concern now.
  //
idx_t encode_samples_fast(binary_matrix& E, const binary_matrix& D,  binary_matrix& A) 
{
  //  std::cout << "cu/fast" << std::endl;
  const idx_t m = E.get_cols();
  const idx_t n = E.get_rows();
  const idx_t p = D.get_rows();
  //
  // SPARSE CODING STEP
  // Go over each row and encode it using rows from D until the residual weight
  // is no longer diminished.
  //
  idx_t changed = 0;
  idx_t G[p][p];
  idx_t r[p];
  idx_t h[p];
  binary_matrix dk(1,m),eik(1,m);
  binary_matrix& d1 = dk; // aliases
  binary_matrix& d2 = eik;
  for (idx_t k1 = 0; k1 < p; ++k1) {
    D.copy_row_to(k1,d1);
    for (idx_t k2 = k1; k2 < p; ++k2) {
      D.copy_row_to(k2,d2);
      bool_and(d1,d2,d2);
      G[k1][k2] = G[k2][k1] = d2.weight();
    }
  }
#if 0
  for (idx_t k1 = 0; k1 < p; ++k1) {
    for (idx_t k2 = 0; k2 < p; ++k2) {
      std::cout << " k1=" << k1 << " k2=" << k2 << " G=" << G[k1][k2] << std::endl;
    }
  }
#endif
  binary_matrix Ei(1,m);
  binary_matrix Ai(1,p); 
  for (idx_t i = 0; i < n; i++) {
    E.copy_row_to(i,Ei);
    A.copy_row_to(i,Ai);
    idx_t w = Ei.weight();
    // EitD^(0)
    // h(0)
    bool ichanged = false;
    idx_t iter = 0;
    for (idx_t k = 0; k < p; ++k) {
      D.copy_row_to(k,dk);
      bool_and(dk,Ei,eik);
      r[k] = eik.weight();
      h[k] = G[k][k] - 2*r[k] + w;
      //std::cout << "\titer=" << iter << " k=" << k1 << " r=" << r[k1] << " d=" << h[k1] << std::endl;
    }
    // a(0) = 0
    bool improved = true;
    while (improved) {
      const idx_t bestk = std::min_element(h,h+p) - h;
      const idx_t bestd = h[bestk];
      //std::cout << "i=" << i << " w=" << w << " bestk=" << bestk << " bestd=" << bestd << std::endl;
      if (bestd < w) {
	idx_t olda = Ai.get(0,bestk);
	Ai.flip(0,bestk);
	D.copy_row_to(bestk,dk);
	add(Ei,dk,Ei);
	w = bestd;
	idx_t oldrbestk = r[bestk];
	iter++;
	if (olda) { // remove atom from solution
	  for (idx_t k = 0; k < p; ++k) {
	    //std::cout << "G[k][bestk]=" << G[k][bestk] << std::endl;
	    r[k] = r[k] + G[k][bestk];  
	    h[k] = h[k] + G[bestk][bestk] + 2*(((int)oldrbestk-G[k][bestk]));  
	    D.copy_row_to(k,dk); // DEBUG!
	    //std::cout << "\titer=" << iter << " k=" << k << " r=" << r[k] << " d=" << h[k] << " d=" << dist(Ei,dk) << std::endl;
	  }
	} else  { // add atom to solution
	  for (idx_t k = 0; k < p; ++k) {
	    //std::cout << "G[k][bestk]=" << G[k][bestk] << std::endl;
	    r[k] = r[k] - G[k][bestk];  
	    h[k] = h[k] + G[bestk][bestk] - 2*(((int)oldrbestk-G[k][bestk]));  
	    D.copy_row_to(k,dk); // DEBUG!
	    //std::cout << "\titer=" << iter << " k=" << k << " r=" << r[k] << " d=" << h[k] << " d=" << dist(Ei,dk) << std::endl;
	  }
	}
	ichanged = true;
      } else {
	improved = false;
      } 
    } // while there is any improvement
    if (ichanged) {
      changed++;
      E.set_row(i,Ei);
      A.set_row(i,Ai);
    }
  }
  dk.destroy();
  eik.destroy();
  Ei.destroy();
  Ai.destroy();
  return changed;
}
