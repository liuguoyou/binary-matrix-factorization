#include "binmat.h"
#include "pbm.h"
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include "bsvd.h"
#include "random_number_generation.h"
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

/**
 * KSVD-like binary dictionary learning algorithm applied to
 * image patches.
 */
int main(int argc, char **argv) {   
  idx_t rows,cols;
  int res;
  FILE* fimg;
  parse_args(argc,argv);
  learn_model_setup(mi_algo,es_algo,du_algo,lm_algo,lmi_algo);
  fimg = fopen(iname,"r");
  set_grid_width(W);
  if (!fimg) return -1;
  res = read_pbm_header(fimg,rows,cols);
  std::cout << "rows=" << rows << " cols=" << cols << std::endl;

  //
  // input data
  // 
  binary_matrix I(rows,cols);
  read_pbm_data(fimg,I);
  if (res !=PBM_OK) { std::cerr << "Error " << res << " reading image."  << std::endl; std::exit(1); }
  fclose(fimg);

  idx_t M,N;
  binary_matrix X,H;
  if (image_mode) {
    std::cout << "==== DATA TREATED AS IMAGE, VECTORS ARE PATCHES =====\n" << std::endl;
    idx_t Ny = (W-1+rows)/W;
    idx_t Nx = (W-1+cols)/W;
    M = W*W;
    std::cout << "Nx=" << Nx << " Ny=" << Ny << std::endl;
    N = Nx*Ny;
    X.allocate(N,M);
    //
    // Initialize data
    //
    idx_t li = 0;
    binary_matrix P(W,W),V(1,W*W);
    for (idx_t i = 0; i < Ny; i++) {
      for (idx_t j = 0; j < Nx; j++,li++) {
	I.copy_submatrix_to(i*W,(i+1)*W,j*W,(j+1)*W,P);
	P.copy_vectorized_to(V);
	X.set_row(li,V);
      }
    }
    N = li;
    P.destroy();
    V.destroy();
  } else {
    std::cout << "==== DATA TREATED AS MATRIX, VECTORS ARE ROWS =====\n" << std::endl;
    X = I.get_copy();
    M = I.get_cols();
    N = I.get_rows();
  }
  binary_matrix D(K,M);
  binary_matrix A(N,K);
  std::cout << "M=" << M << " N=" << N << " K=" << K << std::endl;

  //
  // Initialize dictionary
  //
  initialize_dictionary(X,H,D,A);
  binary_matrix E(N,M);
  //
  //  2. learn model
  //
  idx_t L = learn_model(X,H,E,D,A);
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
	I.set_submatrix(i*W,j*W,P);
      }
    }  
    P.destroy();
    V.destroy();
    fimg = fopen("residual.pbm","w");
    if (!fimg) return -2;
    write_pbm(I,fimg);
    fclose(fimg);
  } else {
    if (force_mosaic)
      render_mosaic(D,"atoms_mosaic.pbm");
  }
  if (force_residual_mosaic) {
    render_mosaic(E,"residual_mosaic.pbm");
  }
  mul(A,false,D,false,E);
  add(E,X,E);
  std::cout << "FINAL: |E|=" << E.weight()
	    << " |A|=" << A.weight() 
	    << " |D|=" << D.weight() 
	    << " L(X)=" << L
	    << " K=" << D.get_rows() << std::endl;
  A.destroy();
  I.destroy();
  E.destroy();
  D.destroy();
  X.destroy();
  return 0;
}
