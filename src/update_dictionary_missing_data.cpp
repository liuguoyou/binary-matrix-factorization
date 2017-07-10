#include "update_dictionary.h"
#include <omp.h>

idx_t update_dictionary_steepest_missing_data(binary_matrix& E,
					      const binary_matrix& H,
					      binary_matrix& D,
					      binary_matrix& A)
{
  //
  // DICTIONARY UPDATE STEP
  // For each atom,
  // take the samples that use it
  // if for the j-th dimension the weight is larger than m/2, turn j-th
  // atom element on, else turn it off
  //
  //  std::cout << "du/stee" << std::endl;
  const idx_t m = E.get_cols();
  const idx_t n = E.get_rows();
  const idx_t p = D.get_rows();
  binary_matrix Dk(1,m);
  binary_matrix newDk(1,m);
  binary_matrix Ei(1,m);

  idx_t changed = 0;
  for (idx_t k = 0; k < p; k++) {
    D.copy_row_to(k,Dk);
    idx_t weights[(int)m];
    idx_t atom_usage = 0;
    std::fill(weights,weights+m,0);
    for (idx_t i = 0; i < n; i++) {
      idx_t aik = A.get(i,k);
      if (aik) { // if atom k was used in i-th samples
	atom_usage++;
	E.copy_row_to(i,Ei);
	add(Ei,Dk,Ei); // add-back old atom, we do not take it into account
	for (idx_t j = 0; j < m ; j++) {
	  if (Ei.get(0,j)) {
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
  } // for each atom
  newDk.destroy();
  Dk.destroy();
  Ei.destroy();
  return changed;
} // end


idx_t update_dictionary_proximus_missing_data(binary_matrix& E,
					      const binary_matrix& H,
					      binary_matrix& D,
					      binary_matrix& A)
{
  //
  // DICTIONARY UPDATE STEP
  // For each atom dk,
  // take the subset of samples that use it
  // remove dk from them
  // apply the proximus rank-one approximation algorithm to dk and a^k
  // 
  //  std::cout << "du/prox" << std::endl;
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
    //
    // PROXIMUS INITIALIZATION of Ak: maximize correlation, not distance
    //
    // ----------------
#if 0
    A.copy_col_to(k,Ak);
    idx_t da = 0;
    //
    // i. compute s = <E,D_k^t>
    //
    std::pair<idx_t,idx_t> s[(int)n];
    for (idx_t i = 0 ; i < n; i++) {
      s[i].first = 0; s[i].second = i;
    }
    for (idx_t j = 0; j < m; j++) {
      idx_t dkj = D.get(k,j);
      if (dkj) { // if atom k was used in i-th samples
	E.copy_col_to(j,Ej);
	add(Ej,Ak,Ej); // add-back old coefs, we do not take it into account
	for (idx_t i = 0; i < n ; i++) {
	  if (Ej.get(0,i)) {
	    s[i].first++;
	  }
	} // add i-th sample to weights vector
      } // if atom was used
    } // for each data sample
      //
      // ii. find A_k that maximizes |A_k^t s| / |A_k|
      // this is simply setting the p largest elements of E D_k^t, where p is such that 
      // z_(p+1) < \sum_{k=1}^{p} z_k / p, where z is the sorted version of s
      // THIS SHOULD BE REPLACED WITH COUNTING SORT, WHICH IS O(n)
      // O(nlog n) kills the algorithm
    std::sort(s,s+n);    
    std::cout << s[n-1].first << std::endl;
    //counting_sort(s,n);
    newAk.clear();
    for (int i = (n-1), sp = 0; (i>=0) && (s[i].first > sp); i--) {
      newAk.set(0,s[i].second);
      sp += s[i].first;
    }
    da = dist(newAk,Ak);
    if (da > 0) { 
      // there was a change 
      A.set_col(k,newAk);
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
    } else {
      continue;
    }
    
#endif
    // ----------------

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
	  //	  kchanged = true;
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


idx_t update_dictionary_steepest_missing_data_omp(binary_matrix& E,
						  const binary_matrix& H,
						  binary_matrix& D,
						  binary_matrix& A)
{
  //
  // DICTIONARY UPDATE STEP
  // For each atom,
  // take the samples that use it
  // if for the j-th dimension the weight is larger than m/2, turn j-th
  // atom element on, else turn it off
  //
  //  std::cout << "du/stee" << std::endl;
  const idx_t m = E.get_cols();
  const idx_t n = E.get_rows();
  const idx_t p = D.get_rows();
  binary_matrix Dk(1,m);
  binary_matrix newDk(1,m);

  omp_set_num_threads(omp_get_max_threads());
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
  const binary_matrix cA = A;
  //  binary_matrix Ei(1,m);

  idx_t changed = 0;
  for (idx_t k = 0; k < p; k++) {
    const idx_t T = omp_get_thread_num();    
    D.copy_row_to(k,Dk);
    idx_t weights[(int)m];
    idx_t atom_usage = 0;
    std::fill(weights,weights+m,0);
#pragma omp parallel for 
    for (idx_t i = 0; i < n; i++) {
      if (cA.get(i,k)) { // if atom k was used in i-th samples
	atom_usage++;
	E.copy_row_to(i,Ei[T]);
	add(Ei[T],Dk,Ei[T]); // add-back old atom, we do not take it into account
	for (idx_t j = 0; j < m ; j++) {
	  if (Ei[T].get(0,j)) {
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
#pragma omp parallel for 
      for (idx_t i = 0; i < n; i++) {
	if (!cA.get(i,k)) 
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


idx_t update_dictionary_proximus_missing_data_omp(binary_matrix& E,
						  const binary_matrix& H,
						  binary_matrix& D,
						  binary_matrix& A)
{
  //
  // DICTIONARY UPDATE STEP
  // For each atom dk,
  // take the subset of samples that use it
  // remove dk from them
  // apply the proximus rank-one approximation algorithm to dk and a^k
  // 
  //  std::cout << "du/prox" << std::endl;
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
    //
    // PROXIMUS INITIALIZATION of Ak: maximize correlation, not distance
    //
    // ----------------
#if 0
    A.copy_col_to(k,Ak);
    idx_t da = 0;
    //
    // i. compute s = <E,D_k^t>
    //
    std::pair<idx_t,idx_t> s[(int)n];
    for (idx_t i = 0 ; i < n; i++) {
      s[i].first = 0; s[i].second = i;
    }
    for (idx_t j = 0; j < m; j++) {
      idx_t dkj = D.get(k,j);
      if (dkj) { // if atom k was used in i-th samples
	E.copy_col_to(j,Ej);
	add(Ej,Ak,Ej); // add-back old coefs, we do not take it into account
	for (idx_t i = 0; i < n ; i++) {
	  if (Ej.get(0,i)) {
	    s[i].first++;
	  }
	} // add i-th sample to weights vector
      } // if atom was used
    } // for each data sample
      //
      // ii. find A_k that maximizes |A_k^t s| / |A_k|
      // this is simply setting the p largest elements of E D_k^t, where p is such that 
      // z_(p+1) < \sum_{k=1}^{p} z_k / p, where z is the sorted version of s
      // THIS SHOULD BE REPLACED WITH COUNTING SORT, WHICH IS O(n)
      // O(nlog n) kills the algorithm
    std::sort(s,s+n);    
    std::cout << s[n-1].first << std::endl;
    //counting_sort(s,n);
    newAk.clear();
    for (int i = (n-1), sp = 0; (i>=0) && (s[i].first > sp); i--) {
      newAk.set(0,s[i].second);
      sp += s[i].first;
    }
    da = dist(newAk,Ak);
    if (da > 0) { 
      // there was a change 
      A.set_col(k,newAk);
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
    } else {
      continue;
    }
    
#endif
    // ----------------

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
	  //	  kchanged = true;
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
