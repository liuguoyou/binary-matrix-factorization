#include "gsl/gsl_sf_gamma.h"
#include "entropy_coding.h"
#include <cmath>

void bilinear_predictor(const binary_matrix& P, binary_matrix& pP) {
  // can be very quickly implemented at block level using binary operators,
  // but I have to sit down and work the math. It's very easy actually...
  for (idx_t j = 1; j < P.get_cols();++j) {
    pP.set(0,j,XOR(P.get(0,j-1),P.get(0,j))); // first row is order 0 pred
  }
  for (idx_t i = 1; i < P.get_rows();++i) {
    pP.set(i,0,XOR(P.get(i-1,0),P.get(i,0))); // first row is order 0 pred      
    for (idx_t j = 1; j < P.get_cols();++j) {
      pP.set(i,j,XOR(XOR(P.get(i-1,j-1),P.get(i,j-1)),XOR(P.get(i-1,j),P.get(i,j))) ); // first row is order 0 pred      
    }
  }
}

double enumerative_codelength(const unsigned n,
			      const unsigned r) { // number of nonzeroes
  return r>0 ? gsl_sf_lnchoose ( n,r ) *COSMOS_LOG2E : 0.0;
}

double universal_codelength(const unsigned n, const unsigned r) {
  double p1 = (double)r/(double)n;
  //  std::cout << "n= " << n << "  r=" << r << " p1=" << p1 << std::endl;
  if ((r > 0) && (r < n)) { 
    return double(n)*(-p1*log2(p1)-(1.0-p1)*log2(1.0-p1)) + 0.5*log2(n);
  } else {
    return 0.5*log2(n);
  }
}

idx_t model_codelength(const binary_matrix& E, 
		       const binary_matrix& D, 
		       const binary_matrix& A) {

  const idx_t M = E.get_cols();
  const idx_t N = E.get_rows();
  const idx_t K = D.get_rows();


  binary_matrix Dk(1,M);
  binary_matrix Ak(1,N);

  idx_t LE = universal_codelength(E.get_rows()*E.get_cols(),E.weight());
  idx_t LD = 0, LA = 0;
  for (idx_t k = 0; k < K; k++) {
    D.copy_row_to(k,Dk);
    A.copy_col_to(k,Ak);
    LD += universal_codelength(M,Dk.weight());
    LA += universal_codelength(N,Ak.weight());
  }
  Dk.destroy();
  Ak.destroy();
  return LE+LD+LA;
}