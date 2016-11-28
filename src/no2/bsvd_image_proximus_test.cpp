#include "binmat.h"
#include "pbm.h"
#include <cstdio>
#include <cstdlib>
#include "bsvd.h"
#include "gsl/gsl_randist.h"
#include <iomanip>
#include "util.h"

int mi_algo = 0;
int cu_algo = 0;
int du_algo = 0;
int lm_algo = 0;
idx_t W = 16;
idx_t K = 512;

const char* iname = "data/test.pbm";

void parse_args(int argc, char **argv) {
  for (int i = 0; i < argc; ++i) {		       
    if (argv[i][0] == '-') { 
      if (i == (argc-1)) {
	std::cerr << "Missing argument for " << argv[i] << std::endl; exit(-1); 
      }
      const char* val = argv[i+1];
      std::cout << "Parameter " << argv[i] << " value " << val << std::endl;
      switch (argv[i][1]) {
      case 'i': case 'I': mi_algo = atoi(val); break;
      case 'c': case 'C': cu_algo = atoi(val); break;
      case 'd': case 'D': du_algo = atoi(val); break;
      case 'l': case 'L': lm_algo = atoi(val); break;
      case 'w': case 'W': W = (idx_t) atoi(val); break;
      case 'k': case 'K': K = (idx_t) atoi(val); break;
      case 'r': case 'R': random_seed = atol(val); break;
      default: std::cerr << "Invalid option " << argv[i] << std::endl; exit(-1);
      }
      i++;
    } else {
      iname = argv[i];
    }
  }
}

/**
 * KSVD-like binary dictionary learning algorithm applied to
 * image patches.
 */
int main(int argc, char **argv) {   
  idx_t rows,cols;
  int res;
  FILE* fimg;
  parse_args(argc,argv);
  learn_model_setup(mi_algo,cu_algo,du_algo,lm_algo);
  fimg = fopen(iname,"r");
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
  binary_matrix E(N,M);
  mul(A,false,D,false,E);
  add(E,X,E);
#if 0
  //
  // RUN BSVD
  //
  idx_t iter = 0;
  std::cout << "iter=" << std::setw(8) << iter 
	    << "\t||E||=" << std::setw(8) << E.weight()
	    << "\t||D||=" << std::setw(8) << D.weight()
	    << "\t||A||=" << std::setw(8) << A.weight() << std::endl;
  while (changed > 0) {    
    iter++;
    idx_t changed_coefs = update_coefficients(E,D,A);
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
#else
  learn_model(E,D,A);
#endif
#if 0
  FILE* fatom;
  for (idx_t k = 0; k < K; k++) {
    char aux[128];
      V = D.get_row(k);
      P.set_vectorized(V);
      snprintf(aux,128,"atom_%04lu.pbm",k);
      fatom = fopen(aux,"w");
      write_pbm(P,fatom);
      fclose(fatom);
  }
#else
  render_mosaic(D,"mosaic.pbm");
#endif
  li = 0;
  for (idx_t i = 0; i < Ny; i++) {
    for (idx_t j = 0; j < Nx; j++,li++) {
      //     std::cout << "n=" << li << std::endl;
      V = E.get_row(li);
      P.set_vectorized(V);
      I2.set_submatrix(i*W,j*W,P);
    }
  }  
  fimg = fopen("residual.pbm","w");
  if (!fimg) return -2;
  write_pbm(I2,fimg);
  fclose(fimg);
  mul(A,false,D,false,E);
  add(E,X,E);
  std::cout << "|E|" << E.weight() << std::endl;
  E.destroy();
  V.destroy();
  D.destroy();
  X.destroy();
  P.destroy();
  I.destroy();
  I2.destroy();
  return 0;
}
