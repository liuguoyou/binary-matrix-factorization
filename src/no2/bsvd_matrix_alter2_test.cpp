#include "binmat.h"
#include "pbm.h"
#include <cstdio>
#include <cstdlib>
#include "bsvd.h"
#include <iomanip>
#include "util.h"

/**
 * Alternant MOD-like binary dictionary learning algorithm applied to
 * a binary matrix (stored as a PBM image file, but treated as a plain matrix)
 *
 * Roles of A and D are switched after each iteration.
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
  binary_matrix P,V;
  const idx_t N = rows;
  const idx_t M = cols;
  const idx_t K = argc > 2 ? atoi(argv[2]) : 32;
  binary_matrix D(K,M);
  binary_matrix A(N,K);
  binary_matrix E(N,M);
  mul(A,false,D,false,E);
  add(E,X,E);

  binary_matrix Dt(M,K);
  binary_matrix At(K,N);
  binary_matrix Et(M,N);

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

  idx_t outer_changed = 1;
  while (outer_changed > 0) {
    outer_changed = 0;
    while (changed > 0) {    
      iter++;
      //      idx_t changed_coefs = update_coefficients(E,D,A);
      idx_t changed_coefs = update_coefficients_omp(E,D,A);
      std::cout << "DIRECT: iter=" << std::setw(8) << iter 
		<< "\t||E||=" << std::setw(8) << E.weight()
		<< "\t||D||=" << std::setw(8) << D.weight()
		<< "\t||A||=" << std::setw(8) << A.weight()
		<< "\tchanged coefs=" << std::setw(8) << changed_coefs << std::endl;
      changed = update_dictionary(E,D,A);
      std::cout << "DIRECT: iter=" << std::setw(8) << iter 
		<< "\t||E||=" << std::setw(8) << E.weight()
		<< "\t||D||=" << std::setw(8) << D.weight()
		<< "\t||A||=" << std::setw(8) << A.weight()
		<< "\tchanged atoms=" << std::setw(8) << changed << std::endl;
      outer_changed += changed;
    }
    A.transpose_to(At);
    D.transpose_to(Dt);
    E.transpose_to(Et);
    changed = 1;
    iter = 0;
    while (changed > 0) {    
      iter++;
      idx_t changed_coefs = update_coefficients_omp(Et,At,Dt);
      std::cout << "TRANSPOSED: iter=" << std::setw(8) << iter 
		<< "\t||E||=" << std::setw(8) << Et.weight()
		<< "\t||D||=" << std::setw(8) << Dt.weight()
		<< "\t||A||=" << std::setw(8) << At.weight()
		<< "\tchanged coefs=" << std::setw(8) << changed_coefs << std::endl;

      changed = update_dictionary(Et,At,Dt);
      std::cout << "TRANSPOSED: iter=" << std::setw(8) << iter 
		<< "\t||E||=" << std::setw(8) << Et.weight()
		<< "\t||D||=" << std::setw(8) << Dt.weight()
		<< "\t||A||=" << std::setw(8) << At.weight()
		<< "\tchanged atoms=" << std::setw(8) << changed << std::endl; 
      outer_changed += changed;   
    }
    At.transpose_to(A);
    Dt.transpose_to(D);
    Et.transpose_to(E);
  }

  write_pbm(D,"D.pbm");
  write_pbm(A,"A.pbm");
  write_pbm(E,"E.pbm");

  Et.destroy();
  At.destroy();
  Dt.destroy();
  A.destroy();
  D.destroy();
  E.destroy();
  V.destroy();
  X.destroy();
  P.destroy();
  return 0;
}
