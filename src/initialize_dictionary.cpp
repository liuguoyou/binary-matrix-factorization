#include "initialize_dictionary.h"

void initialize_dictionary_random_centroids_xor(const binary_matrix& E, 
					   binary_matrix& D, 
					   binary_matrix& A) {
  //
  // Initialize using random clusters
  //
  //
  // Initialize dictionary
  //
  const idx_t m = E.get_cols();
  const idx_t n = E.get_rows();
  const idx_t p = D.get_rows();
  binary_matrix Dk(1,m);
  binary_matrix Ei(1,m);
  gsl_rng *rng = get_rng();  // random number generator
  A.clear();
  D.clear();
  for (idx_t i = 0; i < n; ++i) {
    idx_t k = (idx_t) gsl_rng_uniform_int(rng,p);
    A.set(i,k);
    D.copy_row_to(k,Dk);
    E.copy_row_to(i,Ei);
    add(Dk,Ei,Dk);
    D.set_row(k,Dk);
  }
  Ei.destroy();
  Dk.destroy();
}

void initialize_dictionary_random_centroids(const binary_matrix& E, 
				       binary_matrix& D, 
				       binary_matrix& A) {
  //
  // Initialize using random clusters
  //
  //
  // Initialize dictionary
  //
  const idx_t m = E.get_cols();
  const idx_t n = E.get_rows();
  const idx_t p = D.get_rows();
  binary_matrix Dk(1,m);
  binary_matrix Ei(1,m);
  gsl_rng *rng = get_rng();  // random number generator
  A.clear();
  D.clear();
  idx_t s[p][m];
  idx_t u[p];
  std::fill(&s[0][0],&s[0][0]+m*p,0);
  std::fill(u,u+p,0);
  idx_t S = 0;
  for (idx_t i = 0; i < n; ++i) {
    idx_t k = (idx_t) gsl_rng_uniform_int(rng,p);
    A.set(i,k);
    u[k]++;
    E.copy_row_to(i,Ei);
    for (idx_t j = 0; j < m; j++)
      if (Ei.get(0,j)) { S++; s[k][j]++; }
  }
  S /= (m*n);
  for (idx_t k = 0; k < p; k++) {
    for (idx_t j = 0; j < m; j++) {
      D.set(k,j, 2*s[k][j] >= u[k]);
    }
  }
  Ei.destroy();
  Dk.destroy();
}

/**
 * Works with up to m atoms in principle,
 * Each atom Dk is initialized to the centroid of all rows for which E
 * has ones on its k-th column
 */
void initialize_dictionary_partition(const binary_matrix& E, 
				binary_matrix& D, 
				binary_matrix& A) {
  //
  // Initialize using random clusters
  //
  //
  // Initialize dictionary
  //
  const idx_t m = E.get_cols();
  const idx_t n = E.get_rows();
  const idx_t p = D.get_rows();
  binary_matrix Dk(1,m);
  binary_matrix Ei(1,m);
  A.clear();
  D.clear();
  idx_t s[m];
  aux_t ranking[m];
  for (idx_t k = 0; k < m ; k++) {
    ranking[k].first = E.col_weight(k);
    ranking[k].second = k;
  }
  counting_sort(ranking,m);

  idx_t u;
  for (idx_t k = 0; (k < p) && (k < m); k++) {
    u = 0;
    std::fill(s,s+m,0);
    idx_t pivot = ranking[m-k-1].second;
    //    std::cout << "k=" << k << " pivot=" << pivot << " score=" << ranking[m-k-1].first << std::endl;
    for (idx_t i = 0; i < n; ++i) {
      if (E.get(i,pivot)) {
	u++;
	E.copy_row_to(i,Ei);
	for (idx_t j = 0; j < m; ++j) {
	  if (Ei.get(0,j))
	    s[j]++;
	}
      } // Ei had a row       
    }
    for (idx_t j = 0; j < m; ++j) 
      D.set(k,j,s[j] >= u/2);
  }
  // if p > m, the other atoms are left uninitialized!
  Ei.destroy();
  Dk.destroy();
}

/**
 * For each atom Dk we choose a random sample from E, Ej
 * and compute Dk as the centroid of all samples Ei that are neighbors to Ej in some way
 * in Boolean algebra, and PROXIMUS, this is interpreted as all other samples Ej that have a non-empty
 * intersection with Ei, that is, such that Ej AND Ei is non-zero. 
 */
void initialize_dictionary_neighbor(const binary_matrix& E, 
			       binary_matrix& D, 
			       binary_matrix& A) {
  const idx_t m = E.get_cols();
  const idx_t n = E.get_rows();
  const idx_t p = D.get_rows();
  binary_matrix Dk(1,m);
  binary_matrix Ei(1,m);
  binary_matrix Ej(1,m);
  gsl_rng *rng = get_rng();  // random number generator
  A.clear();
  D.clear();
  for (idx_t k = 0; k < p; ) {
    // pick a random row
    idx_t i = gsl_rng_uniform_int(rng,n);
    E.copy_row_to(i,Ei);
    if (Ei.weight()==0) continue;
    idx_t s[m];
    idx_t u = 0;
    std::fill(s,s+m,0);
    for (idx_t j = 0; j < n; ++j) {
      E.copy_row_to(j,Ej);
      bool_and(Ej,Ei,Ej);
      if (Ej.weight() > 0) {
	u++;
	for (idx_t j = 0; j < m; ++j) {
	  if (Ej.get(0,j))
	    s[j]++;
	}
      }
    }
    if (u > 0) { // the chosen has neighbors
      for (idx_t j = 0; j < m; ++j) 
	D.set(k,j,s[j] >= u/2);
      k++;
    }
  }
  Ei.destroy();
  Ej.destroy();
  Dk.destroy();
}

idx_t accumulate_to(binary_matrix& v, idx_t* s) {
  const idx_t m = v.get_cols();
  idx_t S=0;
  for (idx_t i = 0; i < m; i++) {
    if (v.get(0,i)) { 
      s[i]++;
      S++;
    }
  }
  return S;
}
/**
 * For each atom Dk we define a 'subgraph' and initialize it to a random row Ei, 
 * Then grow the subgraph by adding those rows Ej which share support with any element in the subgraph
 * I guess a good idea is to add not all, but
 * Finally, set Dk as the centroid of all samples in its corresponding part.
 */
void initialize_dictionary_graph_grow(const binary_matrix& E, 
				 binary_matrix& D, 
				 binary_matrix& A) {
  const idx_t m = E.get_cols();
  const idx_t n = E.get_rows();
  const idx_t p = D.get_rows();
  binary_matrix Dk(1,m);
  binary_matrix Ei(1,m);
  gsl_rng *rng = get_rng();  // random number generator
  A.clear();
  D.clear();
  idx_t s[p][m]; // part representation
  idx_t u[p]; // elements in each part
  bool t[n]; //mark if a given sample was selected
  std::fill(&s[0][0],&s[0][0]+m*p,0);
  std::fill(u,u+p,0);
  std::fill(t,t+n,false);
  //
  // initialize p 'subgraphs'. graphs are represented by a counts vector, which adds up 
  // the vectors in the graph.
  //
  int left = n;
  for (idx_t k = 0; (left >= 0) && (k < p); ) {
    // pick a random row
    idx_t i;
    do { 
      i = gsl_rng_uniform_int(rng,n);
    } while (t[i]);
    E.copy_row_to(i,Ei);
    for (idx_t j = 0; j <m; j++) {
      if (E.get(i,j)) s[k][j]=1; 
    }
    t[i] = true;
    left--;
    u[k] = 1;
    k++;
  }
  //
  // for each part, add to it the previously unused pattern that is closest
  // to its support. The original algorithm means ANY row that shares ANY non-zero 
  // with the part. HEre we pivot by weight: we prefer those that share the most used
  // atoms within that part.
  //
  while (left > 0) {
    for (idx_t k = 0; k < p; k++) {
      idx_t maxscore = 0;
      idx_t maxi = 0;
      idx_t score = 0;
      //
      // choose best newcomer
      //
      if (left <= 0) 
	break;
      for (idx_t i = 0; left && (i < n); i++) {
	if (t[i]) continue;
	for (idx_t j = 0; j <m; j++) {
	  //	  if (E.get(i,j)) score += s[k][j];  // NOT BAD< NEED TO TEST
	  if (E.get(i,j)) score  = 1; // do not accumulate overlapped supports, just support
	}
	if (score > maxscore) {
	  maxscore = score;
	  maxi = i;
	}
      }
      //std::cout << "left=" << left << " maxscore=" << maxscore << " maxi=" << maxi << " k=" << k << std::endl;
      if (maxscore == 0) { // reset part!
	idx_t i;
	do {
	  i = gsl_rng_uniform_int(rng,n);
	} while (t[i]);
	E.copy_row_to(i,Ei);
	for (idx_t j = 0; j < m; j++) {
	  s[k][j] = Ei.get(0,j) ? 1: 0;
	}
	t[i] = true;
	u[k] = 1;      
	left--;
      } else {
	//
	// add to part
	//
	t[maxi] = true;
	for (idx_t j = 0; j <m; j++) {
	  if (E.get(maxi,j)) s[k][j]++;
	}
	left--;
	u[k]++;
      }
    }
  }
  for (idx_t k = 0; k < p; k++) {
    for (idx_t j = 0; j < m; ++j) 
      D.set(k,j,s[k][j] >= u[k]/2);
  }
  
  Ei.destroy();
  Dk.destroy();
}
  
void initialize_dictionary_random(const binary_matrix& E, 
			     binary_matrix& D, 
			     binary_matrix& A) {
  gsl_rng *rng = get_rng();  // random number generator
  const idx_t K = D.get_rows();
  const idx_t M = D.get_cols();
  for (idx_t k = 0; k < K; k++) {
    for (idx_t j = 0; j < M; j++) {
      D.set(k,j,gsl_ran_bernoulli (rng,0.5));
    }
  }
  A.clear();
}
