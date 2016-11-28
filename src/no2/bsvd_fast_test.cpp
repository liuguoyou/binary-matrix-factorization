#include "binmat.h"
#include "pbm.h"
#include <cstdio>
#include <cstdlib>
#include "bsvd.h"
#include "gsl/gsl_randist.h"
#include <iomanip>

//#define W 7 // diabolic

int main(int argc, char **argv) {   
  idx_t rows,cols;
  int res;
  FILE* fimg;
  fimg = fopen(argc > 1 ? argv[1]: "data/test.pbm","r");
  const idx_t W = argc > 2 ? atoi(argv[2]) : 16;
  set_grid_width(W);
  if (!fimg) return -1;
  res = read_pbm_header(fimg,rows,cols);
  std::cout << "rows=" << rows << " cols=" << cols << std::endl;
  binary_matrix I(rows,cols);
  read_pbm_data(fimg,I);
  if (res !=PBM_OK) { std::cerr << "Error " << res << " reading image."  << std::endl; std::exit(1); }
  fclose(fimg);
  //std::cout << "==== ORIGINAL =====\n" << std::endl;
  //  std::cout << std::endl << I << std::endl;
  binary_matrix P,V;
  idx_t Ny = (W-1+rows)/W;
  idx_t Nx = (W-1+cols)/W;
  idx_t M = W*W;
  std::cout << "Nx=" << Nx << " Ny=" << Ny << std::endl;
  idx_t N = Nx*Ny;
  binary_matrix X(N,M);
  const idx_t K = argc > 3 ? atoi(argv[3]) : 2*M;
  binary_matrix D(K,M);
  binary_matrix A(N,K);
  std::cout << "M=" << M << " N=" << N << " K=" << K << std::endl;
  idx_t li = 0;
  binary_matrix I2(rows,cols);
  I2.clear();
  //
  // Initialize data
  //
  for (idx_t i = 0; i < Ny; i++) {
    for (idx_t j = 0; j < Nx; j++,li++) {
      //     std::cout << "n=" << li << std::endl;
      P = I.get_submatrix(i*W,(i+1)*W,j*W,(j+1)*W);
      V = P.get_vectorized();
      X.set_row(li,V);
    }
  }
  //
  // Initialize dictionary
  //
  initialize_model(X,D,A);
#if 0
  long seed = 34503498;
  gsl_rng *rng;  // random number generator
  rng = gsl_rng_alloc (gsl_rng_rand48);     // pick random number generator
  gsl_rng_set (rng, seed);                  // set seed
  for (idx_t k = 0; k < K; k++) {
    for (idx_t j = 0; j < M; j++) {
      D.set(k,j,gsl_ran_bernoulli (rng,0.5));
    }
  }
  gsl_rng_free (rng);                       // dealloc the rn
#endif
  //
  // RUN BSVD
  //
  idx_t changed = 513;
  binary_matrix E = X.get_copy();
  idx_t iter = 0;
  std::cout << "iter=" << std::setw(8) << iter 
	    << "\t||E||=" << std::setw(8) << E.weight()
	    << "\t||D||=" << std::setw(8) << D.weight()
	    << "\t||A||=" << std::setw(8) << A.weight() << std::endl;
  while (changed > 0) {    
    iter++;
    //    idx_t changed_coefs = update_coefficients(E,D,A);
    idx_t changed_coefs = update_coefficients_fast(E,D,A);
    std::cout << "iter=" << std::setw(8) << iter 
	      << "\t||E||=" << std::setw(8) << E.weight()
	      << "\t||D||=" << std::setw(8) << D.weight()
	      << "\t||A||=" << std::setw(8) << A.weight()
	      << "\tchanged coefs=" << std::setw(8) << changed_coefs << std::endl;
    changed = update_dictionary(E,D,A);
    std::cout << "iter=" << std::setw(8) << iter 
	      << "\t||E||=" << std::setw(8) << E.weight()
	      << "\t||D||=" << std::setw(8) << D.weight()
	      << "\t||A||=" << std::setw(8) << A.weight()
	      << "\tchanged atoms=" << std::setw(8) << changed << std::endl; 
    //    X.copy_to(E);
 }
  // E::cout << "==== RECONSTRUCTED =====\n" << std::endl;
  // std::cout << std::endl << I << std::endl;
  // fimg = fopen("reconstructed.pbm","w");
  // if (!fimg) return -2;
  // write_pbm(I,fimg);
  // fclose(fimg);
  // if (dist(I,I2)) {
  //   std::cout << "DIFFER after first pass!" << std::endl;
  // }
  for (idx_t k = 0; k < K; k++) {
      V = D.get_row(k);
      P.set_vectorized(V);
      std::cout << "atom" << k << P << std::endl;
  }
  li = 0;
  for (idx_t i = 0; i < Ny; i++) {
    for (idx_t j = 0; j < Nx; j++,li++) {
      //     std::cout << "n=" << li << std::endl;
      V = E.get_row(li);
      P.set_vectorized(V);
      I2.set_submatrix(i*W,j*W,P);
    }
  }
  std::cout << "==== RESIDUAL =====\n" << std::endl;
  //std::cout << std::endl << E << std::endl;
  fimg = fopen("residual.pbm","w");
  if (!fimg) return -2;
  write_pbm(I2,fimg);
  fclose(fimg);

  E.destroy();
  V.destroy();
  D.destroy();
  X.destroy();
  P.destroy();
  I.destroy();
  I2.destroy();
  return 0;
}
