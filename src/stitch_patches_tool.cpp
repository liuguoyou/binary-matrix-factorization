/**
 * Takes a matrix whose rows are vectorized square patches  from an image,
 * the image width and height, and reconstructs the image.
 */
 
#include "binmat.h"
#include "pbm.h"
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include "util.h"


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
  std::cout << "rows=" << rows << " cols=" << cols << std::endl;
  binary_matrix X(rows,cols);
  read_pbm_data(fX,X);
  if (res !=PBM_OK) { std::cerr << "Error " << res << " reading image."  << std::endl; std::exit(1); }
  fclose(fX);
  //
  // output data: image matrix
  //
  // we need to create an auxiliary integer matrix of size m x n 
  // to accumulate the occurences of ones at each pixel
  binary_matrix I;
  I.allocate(m,n);
  //
  // finish up
  //
  fX = fopen(oname,"w");
  if (!fX) return -2;
  write_pbm(I,fX);
  fclose(fX);
  I.destroy();
  X.destroy();
  return 0;
}
