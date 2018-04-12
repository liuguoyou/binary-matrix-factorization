#include "binmat.h"
#include "pbm.h"
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include "coefficients_update.h"
#include "random_number_generation.h"
#include "util.h"
#include "config.h"

int main(int argc, char **argv) {   
  idx_t m,n,p,mx;
  int res;
  FILE* fD, *fx;
  idx_t max_e_weight = 0;
  idx_t max_a_weight = 5;
  // 
  // read dictinary
  //
  fD = fopen("data/test/D.pbm","r");
  if (!fD) return -1;
  res = read_pbm_header(fD,p,m);
  if (res !=PBM_OK) { std::cerr << "Error " << res << " reading D.pbm."  << std::endl; std::exit(1); }
  std::cout << "D: p=" << p << " m=" << m << std::endl;
  binary_matrix D(p,m);
  read_pbm_data(fD,D);
  fclose(fD);
  fD = NULL;
  std::cout << "D: " << D << std::endl;
  //
  // read test samples
  //
  fx = fopen("data/test/x.pbm","r");
  if (!fx) return -1;
  res = read_pbm_header(fx,n,mx);
  if (res !=PBM_OK) { std::cerr << "Error " << res << " reading x.pbm."  << std::endl; std::exit(1); }
  std::cout << "D: n=" << n << " m=" << mx << std::endl;
  binary_matrix x(n,mx);
  read_pbm_data(fx,x);
  fclose(fx);
  fx= NULL;
  std::cout << "x: " << x << std::endl;
  
  //
  // encode
  //
  binary_matrix a(n,p);
  binary_matrix e(n,m);
  binary_matrix h;
  x.copy_to(e);
  coefficients_update_corr(e,h,D,a,max_a_weight,max_e_weight); 
  std::cout << "e: " << e << "a: " << a << std::endl;
  //
  // 3. write output
  //
  a.destroy();
  e.destroy();
  x.destroy();
  D.destroy();
  return 0;
}
