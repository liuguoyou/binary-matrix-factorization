/**
 * Take an image, a patch width, and save a matrix whose columns are the image patches
 * as columns.
 */
 
#include "binmat.h"
#include "pbm.h"
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include "util.h"

idx_t W = 16;
const char* iname = "data/einstein.pbm";
const char* oname = "res/einstein_patches.pbm";

void parse_args(int argc, char **argv) {
  if (argc < 4) { 
    std::cerr << "USAGE: " << argv[0] << " image.pbm patch_width patches.pbm" << std::endl;
    exit(-1); 
  }
  iname = argv[1];
  W = atoi(argv[2]);
  oname = argv[3];
  
}

int main(int argc, char **argv) {   
  parse_args(argc,argv);
  //
  // input data
  // 
  idx_t rows,cols;
  int res;
  FILE* fimg;
  fimg = fopen(iname,"r");
  if (!fimg) return -1;
  res = read_pbm_header(fimg,rows,cols);
  std::cout << "Input: " << iname << " rows=" << rows << " cols=" << cols << std::endl;
  binary_matrix I(rows,cols);
  read_pbm_data(fimg,I);
  if (res !=PBM_OK) { std::cerr << "Error " << res << " reading image."  << std::endl; std::exit(1); }
  fclose(fimg);

  binary_matrix X;
  const idx_t Ny = (rows-W+1);
  const idx_t Nx = (cols-W+1);
  const idx_t M = W*W;
  const idx_t N = Nx*Ny;
  std::cout << "Output: " << oname << " M=" << M << "Nx=" << Nx << " Ny=" << Ny << " N=" << N << std::endl;
  X.allocate(N,M);
  //
  // Initialize data
  //
  idx_t li = 0;
  binary_matrix P(W,W),V(1,W*W);
  for (idx_t i = 0; i < Ny; i++) {
    for (idx_t j = 0; j < Nx; j++,li++) {
      I.copy_submatrix_to(i,i+W,j,j+W,P);
      P.copy_vectorized_to(V);
      X.set_row(li,V);
    }
  }
  P.destroy();
  V.destroy();
  I.destroy();
  write_pbm(X,oname);
  X.destroy();
  return 0;
}
