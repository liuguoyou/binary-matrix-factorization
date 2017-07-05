/** 
 * Adds binary memoryless noise to an input matrix.
 * Parameters: input.pbm output.pbm p [q]
 * p is the crossover probability; this parameter is mandatory
 * if q is specified, then p is the probability of turning a 0 into 1
 * and q the probability of turning a 1 into 0. Otherwise the channel
 * is a binary symmetric channel of crossover probability p.
 */
#include <cstdio>
#include <cstdlib>
#include "pbm.h"
#include "random_number_generation.h"

const char* ifname;
const char* ofname;

int main(int argc, char** argv) {
  idx_t type,rows,cols;
  FILE* fi, *fo;
  if (argc < 4) {
    std::cerr << "Wrong number of parameters.\n Usage: input.pgm output.pgm p [q]" << std::endl;
  }
  ifname = argv[1]; 
  fi = fopen(ifname,"r");
  if (!fi) {
    std::cerr << "Error opening " << ifname << " for reading." << std::endl;
    return -1;
  }      
  read_pbm_header(fi,rows,cols);
  std::cout << "rows=" << rows << " cols=" << cols << std::endl;
  ofname = argv[2];
  fo = fopen(ofname,"w");
  if (!fo) { 
    std::cerr << "Error opening " << ofname << " for writing." << std::endl;
    return -2;
  }
  const double p = atof(argv[3]);
  const double q = (argc == 5 ? atof(argv[4]) : p); 
  std::cout << "p=" << p << " q=" << q << std::endl;
  
  binary_matrix X(rows,cols);
  binary_matrix Y(rows,cols);

  read_pbm_data(fi,X);
  fclose(fi);

  for (idx_t i = 0; i < rows; i++) {
    for (idx_t j = 0; j < cols; j++) {
      const bool v = X.get(i,j);
      Y.set( i, j, v ? get_bernoulli_sample(1-q) : get_bernoulli_sample(p));
    }
  }
  write_pbm(Y,fo);
  fclose(fo);
  X.destroy();
  Y.destroy();
}
