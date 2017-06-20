#include "binmat.h"
#include "pbm.h"
#include <cstdio>
#include <cstdlib>
#include "bsvd.h"
#include <iomanip>
#include <cmath>
#include "util.h"

/**
 * Basic MOD-like binary dictionary learning algorithm applied to
 * a binary matrix (stored as a PBM image file, but treated as a plain matrix)
 */
int main(int argc, char **argv) {   
  idx_t rows,cols;
  int res;
  FILE* fimg;
  fimg = fopen(argc > 1 ? argv[1]: "data/marvelt.pbm","r");
  if (!fimg) return -1;
  res = read_pbm_header(fimg,rows,cols);
  std::cout << "N=" << rows << " M=" << cols << std::endl;
  binary_matrix X(rows,cols);
  read_pbm_data(fimg,X);
  if (res !=PBM_OK) { std::cerr << "Error " << res << " reading image."  << std::endl; std::exit(1); }
  fclose(fimg);
  //std::cout << "==== ORIGINAL =====\n" << std::endl;
  const idx_t N = rows;
  const idx_t M = cols;
  const idx_t K = argc > 2 ? atoi(argv[2]) : 32;
  binary_matrix D(K,M);
  binary_matrix A(N,K);
  binary_matrix E(N,M);
  binary_matrix V(1,M), P(sqrt(M),sqrt(M));
  mul(A,false,D,false,E);
  add(E,X,E);

  initialize_model(X,D,A);
  std::cout << "M=" << M << " N=" << N << " K=" << K << std::endl;
  //
  // RUN BSVD
  //
  idx_t changed = K+1;
  idx_t iter = 0;
  std::cout << "iter=" << std::setw(8) << iter 
	    << "\t||E||=" << std::setw(8) << E.weight()
	    << "\t||D||=" << std::setw(8) << D.weight()
	    << "\t||A||=" << std::setw(8) << A.weight() << std::endl;
  while (changed > 0) {    
    iter++;
    idx_t changed_coefs = update_coefficients_omp(E,D,A);
    std::cout << "iter=" << std::setw(8) << iter 
	      << "\t||E||=" << std::setw(8) << E.weight()
	      << "\t||D||=" << std::setw(8) << D.weight()
	      << "\t||A||=" << std::setw(8) << A.weight()
	      << "\tchanged coefs=" << std::setw(8) << changed_coefs << std::endl;
    changed = update_dictionary_proximus(E,D,A);
    std::cout << "iter=" << std::setw(8) << iter 
	      << "\t||E||=" << std::setw(8) << E.weight()
	      << "\t||D||=" << std::setw(8) << D.weight()
	      << "\t||A||=" << std::setw(8) << A.weight()
	      << "\tchanged atoms=" << std::setw(8) << changed << std::endl;
  }
  render_mosaic(D,"digits_prox_mosaic.pgm");
  

  E.destroy();
  V.destroy();
  D.destroy();
  X.destroy();
  P.destroy();
  A.destroy();
  return 0;
}
