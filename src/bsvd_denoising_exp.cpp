/**
 * B-SVD denoising experiment
 * Takes a binary image and an initial wxw image patches dictionary as input
 * Adds Bernoulli noise to the image of known probability of error p (use GSL RNG functions together with get_rng() defined here in random_number_generation.h)
 * Split the image into wxw overlapping patches (use binmat functions for this)
 * Further adapt the initial dictionary using the noisy patches
 * Obtain clean patch estimates by encoding each noisy sample so that the Hamming error falls below w*w*p
 * Estimate the clean image as a superposition of the estimated patches; here the reasonable operation for
 * deciding upon the final value of a pixel is majority vote: use the value that appears most frequently among the estimates for that pixel coming from different
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
idx_t K = 16;
bool image_mode = false;
bool force_mosaic = true;
bool force_residual_mosaic = true;
const char* iname = "data/test.pbm";
const char* dname = 0;
const char* oname = "denoised.pbm";
double error_probability = 0.1;

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
      case 'k': K = (idx_t) atoi(val); break;
      case 'D': dname = val; break;
      case 'o': oname = val; break;
      case 'p': error_probability = atof(val); break;
      case 'r': random_seed = atol(val); break;
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
  idx_t M,N;
  int res;
  FILE* fX;
  parse_args(argc,argv);
  learn_model_setup(mi_algo,es_algo,du_algo,lm_algo,lmi_algo);
  fX = fopen(iname,"r");
  if (!fX) return -1;
  res = read_pbm_header(fX,N,M);
  std::cout << "M=" << M << " N=" << N << std::endl;

  //
  // input data
  // 
  binary_matrix X(N,M);
  read_pbm_data(fX,X);
  if (res !=PBM_OK) {
    std::cerr << "Error " << res << " reading image."  << std::endl; std::exit(1);
    fclose(fX);
    exit(-1);
  }
  fclose(fX);
  

  //
  // Initial dictionary
  //
  binary_matrix D,A;
  if (dname) {
    fX = fopen(dname,"r");
    if (!fX) return -1;
    idx_t Md;
    read_pbm_header(fX,K,Md);
    if (Md != M) {
      std::cerr << "Dictionary dimension " << Md << " does not match data dimension " << M << "." << std::endl;
      fclose(fX);
      exit(-1);
    }
    D.allocate(K,Md);
    read_pbm_data(fX,D);
    if (res !=PBM_OK) {
      std::cerr << "Error " << res << " reading image."  << std::endl; std::exit(1);
      fclose(fX);
      exit(-1);
    }
    fclose(fX);
  } else {
    D.allocate(K,M);
    initialize_dictionary(X,D,A);
  }
  A.allocate(N,K);
  std::cout << "M=" << M << " N=" << N << " K=" << K << std::endl;
  binary_matrix E(N,M);
  //
  //  2. further update dictionary
  //
  learn_model(X,E,D,A);
  //
  // 3. denoise: average number of errors in Bernoulli(p) on a
  //    vector of length M is Mp
  //
  const idx_t me = M*error_probability;
  std::cout << "Average number of errors per row " << me << std::endl;
  encode_samples(X,D,A,K,me);
  mul(A,false,D,false,X); // estimated denoised signal
  //
  // 3. write output
  //
  write_pbm(D,"dictionary.pbm");
  write_pbm(A,"coefficients.pbm");
  write_pbm(E,"residual.pbm");
  write_pbm(X,oname);
  if (force_mosaic)
    render_mosaic(D,"atoms_mosaic.pbm");
  if (force_residual_mosaic) {
    render_mosaic(E,"residual_mosaic.pbm");
  }
  A.destroy();
  E.destroy();
  D.destroy();
  X.destroy();
  return 0;
}
