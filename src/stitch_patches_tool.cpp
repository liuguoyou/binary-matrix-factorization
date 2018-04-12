/**
 * Takes a matrix whose rows are vectorized square patches  from an image,
 * the image width and height, and reconstructs the image.
 */
 
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "binmat.h"
#include "intmat.h"
#include "pbm.h"
#include "util.h"
#include "binimage.h"

const char* iname = "data/patches.pbm";
const char* oname = "res/stitched.pbm";
idx_t m,n,s;

void parse_args(int argc, char **argv) {
  if (argc < 6) { 
    std::cerr << "Usage: " << argv[0] << " patches.pbm rows cols stride image.pbm" << std::endl; 
    exit(-1);
  }
  iname = argv[1];
  n = atoi(argv[2]);
  m = atoi(argv[3]);
  s = atoi(argv[4]);
  oname = argv[5];
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
  std::cout << "Input: " << iname << " rows=" << rows << " cols=" << cols << " patch width=" << W << std::endl;
  std::cout << "Parameters: width=" << W << " stride=" << s << std::endl;
  binary_matrix X(rows,cols);
  read_pbm_data(fX,X);
  if (res !=PBM_OK) { std::cerr << "Error " << res << " reading patches."  << std::endl; std::exit(1);  }
  fclose(fX);
  //
  // output data: image matrix
  //
  std::cout << "Output: " << oname << " rows=" << m << " cols=" << n << std::endl;
  
  // we need to create an auxiliary integer matrix of size m x n 
  // to accumulate the occurences of ones at each pixel
  integer_matrix A(m,n);
  integer_matrix C(m,n);
  binary_matrix  I(m,n);
  stitch_patches(X,m,n,s,EXTRACT_EXACT,I,A,C);
  A.destroy();
  C.destroy();
  X.destroy();
  fX = fopen(oname,"w");
  if (!fX) {
    std::cerr << "Error writing image " << oname << std::endl; return -2;
  }
  write_pbm(I,fX);
  fclose(fX);
  I.destroy();
  return 0;
}
