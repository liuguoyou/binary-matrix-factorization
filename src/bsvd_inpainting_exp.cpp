/**
 * B-SVD inpainting experiment
 * Takes a binary samples matrix, a binary mask of the same size, and an initial  dictionary as input
 * Adapt the dictionary using the input samples, each time using only those samples which are available and the corresponding subdictionary
 * Encode the known parts of each sample using the corresponding subdictionary to obtain encoding coefficients
 * Estimate the missing samples using the corresponding complementary set of rows of the dictionary and the coefficients computed as above
 * patches to which it belongs (usually w*w patches)
 */
 
#include "binmat.h"
#include "pbm.h"
#include <cstdio>
#include <cstdlib>
#include "bsvd.h"
#include "random_number_generation.h"
#include <iomanip>
#include "util.h"
#include "config.h"

int mi_algo = 0;
int es_algo = 0;
int du_algo = 0;
int lm_algo = 0;
int lmi_algo = 0;

idx_t W = 16;
idx_t K = 512;
bool image_mode = false;
bool force_mosaic = true;
bool force_residual_mosaic = true;
const char* iname = "data/test.pbm";

void parse_args(int argc, char **argv) {
  for (int i = 0; i < argc; ++i) {		       
    if (argv[i][0] == '-') { 
      if (i == (argc-1)) {
	std::cerr << "Missing argument for " << argv[i] << std::endl; exit(-1); 
      }
      const char* val = argv[i+1];
//      std::cout << "Parameter " << argv[i] << " value " << val << std::endl;
      switch (argv[i][1]) {
      case 'i': mi_algo = atoi(val); break;
      case 'c': es_algo = atoi(val); break;
      case 'd': du_algo = atoi(val); break;
      case 'l': lm_algo = atoi(val); break;
      case 'L': lmi_algo = atoi(val); break;
      case 'w': W = (idx_t) atoi(val); break;
      case 'k': K = (idx_t) atoi(val); break;
      case 'r': random_seed = atol(val); break;
      case 'I': image_mode = (atoi(val) > 0); break;
      case 'm': force_mosaic = (atoi(val) > 0); break;
      case 'M': force_residual_mosaic = (atoi(val) > 0); break;
      default: std::cerr << "Invalid option " << argv[i] << std::endl; exit(-1);
      }
      i++;
    } else {
      iname = argv[i];
    }
  }
}

int main(int argc, char **argv) {   
  idx_t rows,cols;
  parse_args(argc,argv);
  learn_model_setup(mi_algo,es_algo,du_algo,lm_algo,lmi_algo);
  //
  // input data
  // 
  FILE* fX;
  fX = fopen(iname,"r");
  if (!fX) return -1;
  int res;
  res = read_pbm_header(fX,rows,cols);
  std::cout << "rows=" << rows << " cols=" << cols << std::endl;
  binary_matrix X(rows,cols);
  read_pbm_data(fX,X);
  if (res !=PBM_OK) { std::cerr << "Error " << res << " reading data."  << std::endl; std::exit(1); }
  fclose(fX);

  const idx_t M = X.get_cols();
  const idx_t N = X.get_rows();

  binary_matrix D(K,M);
  binary_matrix A(N,K);
  std::cout << "M=" << M << " N=" << N << " K=" << K << std::endl;

  //
  // Initialize dictionary
  //
  initialize_dictionary(X,D,A);
  binary_matrix E(N,M);
  //
  //  2. learn model
  //
  learn_model(X,E,D,A);
  //
  // 3. write output
  //
  write_pbm(D,"dictionary.pbm");
  write_pbm(A,"coefficients.pbm");
  write_pbm(E,"residual.pbm");
  if (image_mode) {
    render_mosaic(D,"atoms_mosaic.pbm");
    idx_t Ny = (W-1+rows)/W;
    idx_t Nx = (W-1+cols)/W;
    idx_t li = 0;
    binary_matrix P(W,W),V(1,W*W);
    for (idx_t i = 0; i < Ny; i++) {
      for (idx_t j = 0; j < Nx; j++,li++) {
	//     std::cout << "n=" << li << std::endl;
        E.copy_row_to(li,V);
        P.set_vectorized(V);
        X.set_submatrix(i*W,j*W,P);
      }
    }  
    P.destroy();
    V.destroy();
    fX = fopen("residual.pbm","w");
    if (!fX) return -2;
    write_pbm(X,fX);
    fclose(fX);
  } else {
    if (force_mosaic)
      render_mosaic(D,"atoms_mosaic.pbm");
  }
  if (force_residual_mosaic) {
    render_mosaic(E,"residual_mosaic.pbm");
  }
  mul(A,false,D,false,E);
  add(E,X,E);
  std::cout << "|E|" << E.weight() << std::endl;
  A.destroy();
  E.destroy();
  D.destroy();
  X.destroy();
  return 0;
}
