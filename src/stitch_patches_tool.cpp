/**
 * Takes a matrix whose rows are vectorized square patches  from an image,
 * the image width and height, and reconstructs the image.
 */
 
#include "binmat.h"
#include "pbm.h"
#include "util.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>

const char* iname = "data/patches.pbm";
const char* oname = "res/stitched.pbm";
idx_t m,n;

void parse_args(int argc, char **argv) {
  if (argc < 5) { 
    std::cerr << "Usage: " << argv[0] << " patches.pbm m n image.pbm" << std::endl; 
    exit(-1); 
  }
  iname = argv[1];
  m = atoi(argv[2]);
  n = atoi(argv[3]);
  oname = argv[1];
}

int main(int argc, char **argv) {   
  idx_t rows,cols;
  int res;
  //
  // input data: patches matrix, rows are patches
  // 
  FILE* fX;
  parse_args(argc,argv);
  fX = fopen(iname,"r");
  if (!fX) return -1;
  res = read_pbm_header(fX,rows,cols);
  const idx_t W = (idx_t) sqrt(cols);
  std::cout << "Input: " << iname << " << rows=" << rows << " cols=" << cols << " patch width=" << W << std::endl;
  binary_matrix X(rows,cols);
  read_pbm_data(fX,X);
  if (res !=PBM_OK) { std::cerr << "Error " << res << " reading image."  << std::endl; std::exit(1);  }
  fclose(fX);
  //
  // output data: image matrix
  //
  std::cout << "Output: " << oname << " << rows=" << m << " cols=" << n << std::endl;
  
  // we need to create an auxiliary integer matrix of size m x n 
  // to accumulate the occurences of ones at each pixel
  binary_matrix I;
  I.allocate(m,n);
  binary_matrix P(W,W),V(1,W*W);
  int* R = new int[m*n];
  std::fill(R,R+m*n,0);
  const idx_t Ny = m-W+1;
  const idx_t Nx = n-W+1;
  for (idx_t i = 0, li = 0; i < Ny ; i++) {
    for (idx_t j = 0; j < Ny ; j++,li++) {
        X.copy_row_to(li,V);
	for (idx_t i2 = 0; i2 < W; i2++) {
	}
    }
  }
  //
  // finish up
  //
  fX = fopen(oname,"w");
  if (!fX) return -2;
  write_pbm(I,fX);
  fclose(fX);
  P.destroy();
  V.destroy();
  I.destroy();
  X.destroy();
  return 0;
}
