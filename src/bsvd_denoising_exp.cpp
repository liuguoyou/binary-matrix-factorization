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

idx_t K = 512;
bool image_mode = false;
bool force_mosaic = true;
bool force_residual_mosaic = true;
const char* iname = "data/test.pbm";
const char* dname = 0;
const char* oname = "denoised.pbm";
double error_probability = 0.1;

void list_choices(const char* prefix, const char* opts[]) {
  size_t i = 0;
  while (opts[i]) {
    std::cout << prefix << i << " -> " << opts[i] << std::endl;
    i++;
  }
}


void show_help(const char* pname) {
  std::cout << "USAGE: " << pname << " [options] <data_file>" << std::endl;
  std::cout << "where [options can be any combination of: " << std::endl;
  std::cout << "\t-i initialization algorithm. Choices are:" << std::endl;
  list_choices("\t\t",mi_algorithm_names);
  std::cout << "\t-c coefficients update algorithm. Choices are:" << std::endl;
  list_choices("\t\t",es_algorithm_names);
  std::cout << "\t-d dictionary update algorithm. Choices are:" << std::endl;
  list_choices("\t\t",du_algorithm_names);
  std::cout << "\t-l model selection algorithm. Choices are:" << std::endl;
  list_choices("\t\t",lm_algorithm_names);
  std::cout << "\t-L model learning algorithm. Choices are:" << std::endl;
  list_choices("\t\t",lm_algorithm_names);
  std::cout << "\t-w patch width, for image analysis." << std::endl;
  std::cout << "\t-k model size/rank (bool" << std::endl;
  std::cout << "\t-r random number generator seed (int)" << std::endl;
  std::cout << "\t-I model size/rank (bool)" << std::endl;
  std::cout << "\t-m force generation of atom mosaic image (bool)." << std::endl;
  std::cout << "\t-M force generation of residual samples mosaic image (bool)." << std::endl;
  std::cout << "\t-v increase verbosity." << std::endl;
  std::cout << "\t-h show this message." << std::endl;
}

void parse_args(int argc, char **argv) {
  for (int i = 0; i < argc; ++i) {		       
    if (argv[i][0] == '-') { 
      if (argv[i][1] == 'h') {
        show_help(argv[0]); exit(0);
      } else if (argv[i][1] == 'v') {
	inc_verbosity();
	continue;
      } else if (i == (argc-1)) {
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
  idx_t rows,cols;
  FILE* file;
  parse_args(argc,argv);
  learn_model_setup(mi_algo,es_algo,du_algo,lm_algo,lmi_algo);
  file = fopen(iname,"r");
  if (!file) return -1;
  int res = read_pbm_header(file,rows,cols);
  std::cout << "cols=" << cols << " rows=" << rows << std::endl;
  
  //
  // input data
  // 
  binary_matrix X(rows,cols);
  read_pbm_data(file,X);
  if (res !=PBM_OK) {
    std::cerr << "Error " << res << " reading image."  << std::endl; std::exit(1);
    fclose(file);
    exit(-1);
  }
  fclose(file);
  

  //
  // Initialize dictionary
 //
  binary_matrix D,A;
  binary_matrix H;
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
  //    vector of length M is colsp
  //
  const idx_t me = cols*error_probability;
  std::cout << "Denoising. <<" << std::endl;
  std::cout << "Average number of errors per row " << me << std::endl;
  //A.clear();
  //coefficients_update(X,H,D,A,K,me);  
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
