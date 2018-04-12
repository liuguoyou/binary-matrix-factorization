/**
 * Take an image, a patch width, and save a matrix whose columns are the image patches
 * as columns.
 */
 
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include "util.h"
#include "binmat.h"
#include "pbm.h"
#include "binimage.h"

idx_t W = 16;
idx_t s = 1;

const char* iname = "data/einstein.pbm";
const char* oname = "res/einstein_patches.pbm";

void parse_args(int argc, char **argv) {
  if (argc < 4) { 
    std::cerr << "USAGE: " << argv[0] << " image.pbm patch_width stride patches.pbm" << std::endl;
    exit(-1); 
  }
  iname = argv[1];
  W = atoi(argv[2]);
  s = atoi(argv[3]);
  oname = argv[4];
  
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
  std::cout << "Parameters: width=" << W << " stride=" << s << std::endl;
  binary_matrix I(rows,cols);
  read_pbm_data(fimg,I);
  if (res !=PBM_OK) { std::cerr << "Error " << res << " reading image."  << std::endl; std::exit(1); }
  fclose(fimg);

  binary_matrix X;
  const idx_t Ny = compute_grid_size(rows,W,s,EXTRACT_EXACT);
  const idx_t Nx = compute_grid_size(cols,W,s,EXTRACT_EXACT);
  const idx_t M = W*W;
  const idx_t N = Nx*Ny;
  std::cout << "Output: " << oname << " M=" << M << "N=" << N << " Ny=" << Ny << " Nx=" << Nx << std::endl;
  X.allocate(N,M);
  extract_patches(I,W,s,EXTRACT_EXACT,X);
  I.destroy();
  write_pbm(X,oname);
  X.destroy();
  return 0;
}

