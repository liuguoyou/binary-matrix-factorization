/**
 * Classification experiment.
 * Input pbm "image" with binary samples as columns and an ASCII space-separated file with class labels as numbers
 * Procedure:
 * 1) split dataset into training and testing subsets
 * 2) learns a dictionary for  the samples of class c=1,...,C within the training subset; this results in C dictionaries
 * 3) classifies each sample x0 in the testing subset using each of the C dictionaries; this yields a set of scores {l_1,l_2,...,l_C}
 * 4) declares  sample x0 to belong to class c* if c* = argmin {l_c: c=1,...,C}
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
const char* Xname = "data/mnist_data.pbm";
const char* Lname = "data/mnist_labels.ascii";

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
      Xname = argv[i];
      Lname = argv[i+1];
    }
  }
}

int main(int argc, char **argv) {   
  idx_t rows,cols;
  int res;
  FILE* fX, *fL;
  parse_args(argc,argv);
  learn_model_setup(mi_algo,es_algo,du_algo,lm_algo,lmi_algo);
  fX = fopen(Xname,"r");
  set_grid_width(W);
  if (!fX) return -1;
  res = read_pbm_header(fX,rows,cols);
  std::cout << "rows=" << rows << " cols=" << cols << std::endl;

  //
  // input data
  // 
  binary_matrix X(rows,cols);
  read_pbm_data(fX,X);
  if (res !=PBM_OK) { std::cerr << "Error " << res << " reading image."  << std::endl; std::exit(1); }
  fclose(fX);
  char* L = new char[cols];
  int l;
  unsigned i = 0;
  fL = fopen(Lname,"r");
  if (!fL) { std::cerr << "Error  reading labels file " << Lname << std::endl; std::exit(1); }
  while (fscanf(fL,"%d ",&l)) {
    L[i++] = l;
  } 
  fclose(fL);

  //
  // classify
  //

  //
  // cleanup
  //
  X.destroy();
  delete[] L;
  return 0;
}
