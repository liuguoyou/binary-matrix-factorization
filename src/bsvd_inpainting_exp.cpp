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

idx_t K = 512;
bool image_mode = false;
bool force_mosaic = true;
bool force_residual_mosaic = true;
const char* iname = "data/einstein.pbm";
const char* hname = "data/einstein_erasures1.pbm";
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
      case 'r': random_seed = atol(val); break;
      case 'm': force_mosaic = (atoi(val) > 0); break;
      case 'M': force_residual_mosaic = (atoi(val) > 0); break;
      default: std::cerr << "Invalid option " << argv[i] << std::endl; exit(-1);
      }
      i++;
    } else {
      iname = argv[i];
      hname = argv[i+1]; // erasure mask name
    }
  }
}

int main(int argc, char **argv) {   
  idx_t rows,cols;
  parse_args(argc,argv);
  learn_model_setup(mi_algo,es_algo,du_algo,lm_algo,lmi_algo);
  FILE* file = fopen(iname,"r");
  if (!file) return -1;
  int res = read_pbm_header(file,rows,cols);
  std::cout << "cols=" << cols << " rows=" << rows << std::endl;

  //
  // input data
  // 
  binary_matrix X(rows,cols);
  read_pbm_data(file,X);
  if (res !=PBM_OK) {
    std::cerr << "Error " << res << " reading data matrix."  << std::endl; std::exit(1);
    fclose(file);
    exit(-1);
  }
  fclose(file);
  //
  // erasure mask
  //
  idx_t rows2,cols2;
  res = read_pbm_header(file,rows2,cols2);
  std::cout << "cols2=" << cols2 << " rows2=" << rows2 << std::endl;
  if (cols2 != cols) {
	std::cerr << "wrong erasure mask columns: " << cols2 << " should be " << cols << std::endl;
  }
  if (rows2 != rows) {
	std::cerr << "wrong erasure mask rows: " << rows2 << " should be " << rows << std::endl;
  }
  binary_matrix H(rows,cols);
  read_pbm_data(file,H);
  if (res !=PBM_OK) {
    std::cerr << "Error " << res << " reading erasure mask matrix."  << std::endl; std::exit(1);
    fclose(file);
    exit(-1);
  }
  fclose(file);
  

  //
  // Initialize dictionary
  //
  binary_matrix D,A;
  if (dname) {
    file = fopen(dname,"r");
    if (!file) return -1;
    idx_t colsd;
    read_pbm_header(file,K,colsd);
    if (colsd != cols) {
      std::cerr << "Dictionary dimension " << colsd << " does not match data dimension " << cols << "." << std::endl;
      fclose(file);
      exit(-1);
    }
    D.allocate(K,colsd);
    read_pbm_data(file,D);
    if (res !=PBM_OK) {
      std::cerr << "Error " << res << " reading image."  << std::endl; std::exit(1);
      fclose(file);
      exit(-1);
    }
    fclose(file);
  } else {
    D.allocate(K,cols);
    initialize_dictionary(X,H,D,A);
  }
  A.allocate(rows,K);
  A.clear();
  std::cout << "cols=" << cols << " rows=" << rows << " K=" << K << std::endl;
  binary_matrix E(rows,cols);
  if (force_mosaic)
    render_mosaic(D,"denoising_initial_dictionary.pbm");
  //
  //  2. further update dictionary
  //
  std::cout << "Further adapting dictionary to data." << std::endl;
  learn_model(X,H,E,D,A);
  render_mosaic(D,"denoising_adapted_dictionary.pbm");
  //
  // 3. denoise: average number of errors in Bernoulli(p) on a
  //    vector of length cols is colsp
  //
  const idx_t me = cols*error_probability;
  std::cout << "Denoising. <<" << std::endl;
  std::cout << "Average number of errors per row " << me << std::endl;
  A.clear();
  coefficients_update(X,H,D,A,K,me);  
  X.copy_to(E); // at this point X contains the residual
  std::cout << "Average residual weight=" << (double)E.weight()/(double)rows << std::endl;
  std::cout << "Average coefficients weight=" << (double)A.weight()/(double)rows << std::endl;
  mul(A,false,D,false,X); // now X is the estimated denoised signal
  //
  // 3. write output
  //
  write_pbm(D,"dictionary.pbm");
  write_pbm(A,"coefficients.pbm");
  write_pbm(X,oname);
  A.destroy();
  E.destroy();
  D.destroy();
  X.destroy();
  return 0;
}
