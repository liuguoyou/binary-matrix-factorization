#include "bsvd.h"
#include <algorithm>
#include "gsl/gsl_randist.h"
#include <omp.h>


idx_t update_coefficients_omp(binary_matrix& E, const binary_matrix& D,  binary_matrix& A) 
{
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


idx_t update_dictionary_omp(binary_matrix& E, binary_matrix& D, const binary_matrix& A)
{
  //
  // DICTIONARY UPDATE STEP
  // For each atom,
  // take the samples that use it
  // if for the j-th dimension the weight is larger than m/2, turn j-th
  // atom element on, else turn it off
  //
  const idx_t m = E.get_cols();
  const idx_t n = E.get_rows();
  const idx_t p = D.get_rows();
  binary_matrix Dk(1,m);
  binary_matrix newDk(1,m);
  idx_t NT;
#pragma omp parallel
  {  
 NT = omp_get_num_threads(); 
  }
  //  std::cout << "THREADS=" << NT << std::endl;
  binary_matrix Ei[NT];
  for (idx_t TT = 0; TT < NT; TT++) {
    Ei[TT].allocate(1,m);
  }

  idx_t changed = 0;
  for (idx_t k = 0; k < p; k++) {
    D.copy_row_to(k,Dk);
    idx_t weights[(int)m];
    idx_t atom_usage = 0;
    std::fill(weights,weights+m,0);
    //#pragma omp parallel for schedule(dynamic)
    for (idx_t i = 0; i < n; i++) {
      const idx_t T = omp_get_thread_num();
      idx_t aik = A.get(i,k);
      if (aik) { // if atom k was used in i-th samples
	atom_usage++;
	E.copy_row_to(i,Ei[T]);
	add(Ei[T],Dk,Ei[T]); // add-back old atom, we do not take it into account
	for (idx_t j = 0; j < m ; j++) {
	  if (Ei[T].get(0,j)) {
	    //#pragma omp atomic
	    weights[j]++;
	  }
	} // add i-th sample to weights vector
      } // if atom was used
    } // for each data sample
    if (!atom_usage) 
      continue;
    // now update atom based on weights
    const idx_t u = atom_usage/2;
    Dk.copy_to(newDk);
    for (idx_t j = 0; j < m ; j++) {
      newDk.set(0,j, weights[j] > u);
    }
    if (dist(newDk,Dk) > 0) { 
      // there was a change in the atom
      changed++;
      D.set_row(k,newDk);
      // update residual E
#pragma omp parallel for schedule(dynamic)
      for (idx_t i = 0; i < n; i++) {
	const idx_t T = omp_get_thread_num();
	idx_t aik = A.get(i,k);
	if (!aik) 
	  continue;
	E.copy_row_to(i,Ei[T]);
	add(Ei[T],Dk,Ei[T]); // add-back old atom
	add(Ei[T],newDk,Ei[T]); // substract (same as add here) new atom
	E.set_row(i,Ei[T]);
      }
    }
  } // for each atom
  newDk.destroy();
  Dk.destroy();
  for (idx_t T = 0; T < NT; T++) {
    Ei[T].destroy();
  }
  return changed;
} // end


idx_t update_dictionary_proximus_omp(binary_matrix& E, binary_matrix& D, binary_matrix& A)
{
  //
  // DICTIONARY UPDATE STEP
  // For each atom dk,
  // take the subset of samples that use it
  // remove dk from them
  // apply the proximus rank-one approximation algorithm to dk and a^k
  // 
  const idx_t m = E.get_cols();
  const idx_t n = E.get_rows();
  const idx_t p = D.get_rows();
  binary_matrix Ei(1,m); // rows of E
  binary_matrix Ej(1,n); // column of E
  binary_matrix Dk(1,m); // column of D
  binary_matrix Ak(1,n); // row of coefs
  binary_matrix newDk(1,m); // column of D
  binary_matrix newAk(1,n); // row of coefs


  idx_t changed = 0;
  //
  // FOR EACH ATOM
  //
  for (idx_t k = 0; k < p; k++) {
    //
    // PROXIMUS-LIKE ITERATION, UNTIL NO CHANGE OCCURS
    // 
    bool kchanged = false;
    bool converged;
    //    std::cout << "PROXIMUS: k=" << k << std::endl;

    do { // until converged is true
      converged = true;
      //
      // Update atom Dk
      //
      idx_t u = 0;
      idx_t Dw[(int)m];
      std::fill(Dw,Dw+m,0);
      D.copy_row_to(k,Dk);

      for (idx_t i = 0; i < n; i++) {
	idx_t aik = A.get(i,k);
	if (aik) { // if atom k was used in i-th samples
	  u++;
	  E.copy_row_to(i,Ei);
	  add(Ei,Dk,Ei); // add-back old atom, we do not take it into account
	  for (idx_t j = 0; j < m ; j++) {
	    if (Ei.get(0,j)) {
	      Dw[j]++;
	    }
	  } // add i-th sample to weights vector
	} // if atom was used
      } // for each data sample
      idx_t dd = 0;
      if (u) { 
	// now update atom based on weights
	u /= 2; 
	Dk.copy_to(newDk);
	for (idx_t j = 0; j < m ; j++) {
	  newDk.set(0,j, Dw[j] > u);
	}
	dd = dist(newDk,Dk);
	if (dd > 0) { 
	  // there was a change in the atom
	  D.set_row(k,newDk);
	  converged = false;
	  kchanged = true;
	  // update residual matrix E
	  for (idx_t i = 0; i < n; i++) {
	    idx_t aik = A.get(i,k);
	    if (!aik) 
	      continue;
	    E.copy_row_to(i,Ei);
	    add(Ei,Dk,Ei); // add-back old atom
	    add(Ei,newDk,Ei); // substract (same as add here) new atom
	    E.set_row(i,Ei);
	  }
	}
      }
      //      std::cout << "\tu=" << u << "\tdd=" << dd << "\t|E|=" << E.weight() << std::endl;      
      //
      // UPDATE COEFs, Ak
      //
      u = 0;
      A.copy_col_to(k,Ak);
      idx_t Aw[(int)n];
      std::fill(Aw,Aw+n,0);
      idx_t da = 0;
      for (idx_t j = 0; j < m; j++) {
	idx_t dkj = D.get(k,j);
	if (dkj) { // if atom k was used in i-th samples
	  u++;
	  E.copy_col_to(j,Ej);
	  add(Ej,Ak,Ej); // add-back old coefs, we do not take it into account
	  for (idx_t i = 0; i < n ; i++) {
	    if (Ej.get(0,i)) {
	      Aw[i]++;
	    }
	  } // add i-th sample to weights vector
	} // if atom was used
      } // for each data sample
      if (u) {
	// now update atom based on weights
	u /= 2;
	Ak.copy_to(newAk);
	for (idx_t i = 0; i < n ; i++) {
	  newAk.set(0,i, Aw[i] > u);
	}
	da = dist(newAk,Ak);
	if (da > 0) { 
	  // there was a change 
	  A.set_col(k,newAk);
	  converged = false;
	  kchanged = true;
	  // update residual E
	  for (idx_t j = 0; j < m; j++) {
	    idx_t dkj = D.get(k,j);
	    if (!dkj) 
	      continue;
	    E.copy_col_to(j,Ej);
	    add(Ej,Ak,Ej); // add-back old atom
	    add(Ej,newAk,Ej); // substract (same as add here) new atom
	    E.set_col(j,Ej);
	  }
	}
	//	std::cout << "\tu=" << u << "\tda=" << da << "\t|E|=" << E.weight() << std::endl;      
      }
    } while (!converged);
    //
    // the pair Dk, Ak has changed during this run
    //
    if (kchanged)
      changed++;
  } // for each atom
    // finished
  Ei.destroy();
  Ej.destroy();
  newDk.destroy();
  newAk.destroy();
  Dk.destroy();
  Ak.destroy();
  return changed;
} // end


