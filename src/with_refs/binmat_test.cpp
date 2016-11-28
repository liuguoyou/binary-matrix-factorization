#include "binmat.h"
#include "pbm.h"
#include <cstdio>

int main(int argc, char **argv) {
  binary_matrix A(8,128);
  std::cout << A;
  A.clear();
  std::cout << A;
  for (idx_t i = 0 ; i < A.get_rows(); i++) {
    idx_t m1 = 1UL << (i+1);
    idx_t m2 = 1UL << i;
    //    std::cout << "i=" << i << " m1=" << std::bitset<sizeof(idx_t)*8>(m1) << "\tm2="<< std::bitset<sizeof(idx_t)*8>(m2) << std::endl;
    for (idx_t j = 0 ; j < A.get_cols(); j++) {
      A.set(i,j, (j & (m1-1)) & m2 ? 1 : 0);
    }
  }
  std::cout << A.get(5,0) << std::endl;
  std::cout << A;
  std::cout <<"Total Sum:" << A.sum() << std::endl;
  std::cout <<"Total Weight:" << A.weight() << std::endl;
  for (idx_t i = 0 ; i < A.get_rows(); i++) {
    std::cout << "ROW " << i << ": w=" << A.row_weight(i) << " s=" << A.row_sum(i) << std::endl;
  }
  for (idx_t i = 0 ; i < A.get_cols(); i++) {
    std::cout << "COL " << i << ": w=" << A.col_weight(i) << " s=" << A.col_sum(i) << std::endl;
  }

  binary_matrix B(3,2);
  binary_matrix C(2,3);
  binary_matrix D(3,3);
  B.clear();
  C.clear();
  D.clear();
  B.set(0,0,1);
  B.set(1,0,1);
  //  B.set(2,0,1);
  B.set(0,1,1);
  B.set(2,1,1);
  std::cout << "B=" << B << std::endl;
  C.set(0,1,1);
  C.set(1,0,1);
  C.set(1,2,1);
  std::cout << "C=" << C << std::endl;
  mul(B,false,C,false,D);
  std::cout << "BC=" << D << std::endl;
  mul(B,false,B,true,D);
  std::cout << "BBt=" << D << std::endl;
  mul(C,true,C,false,D);
  std::cout << "CtC=" << D << std::endl;

  idx_t rows,cols;
  FILE* fimg;
  fimg = fopen("data/camera.pbm","r");
  if (!fimg) return -1;
  read_pbm_header(fimg,rows,cols);
  std::cout << "rows=" << rows << " cols=" << cols << std::endl;
  binary_matrix I(rows,cols);
  read_pbm_data(fimg,I);
  fclose(fimg);
  fimg = fopen("copy.pbm","w");
  if (!fimg) return -2;
  write_pbm(I,fimg);
  fclose(fimg);
  std::cout << "==== ORIGINAL =====\n" << std::endl;
  std::cout << std::endl << I << std::endl;
  binary_matrix Iref;
  Iref.ref(I);
  std::cout << "==== REFERENCE =====\n" << std::endl;
  std::cout << std::endl << Iref << std::endl;
  std::cout << "==== REF_BLOCK I[20:60,10:40] =====\n" << std::endl;
  Iref.ref_block(I,20,60,10,40);
  std::cout << std::endl << Iref << std::endl;
  Iref.set();
  std::cout << "==== REF_BLOCK I[20:60,10:40]=1 =====\n" << std::endl;
  Iref.ref_block(I,0,90,0,90);
  std::cout << std::endl << Iref << std::endl;
  std::cout << "==== REF_BLOCK I[30:40,40:70] =====\n" << std::endl;

  Iref.ref_block(I,30,40,40,70);
  std::cout << std::endl << Iref << std::endl;
  Iref.clear();
  std::cout << "==== REF_BLOCK I[30:40,40:70]=0 =====\n" << std::endl;
  Iref.ref_block(I,0,90,0,90);
  std::cout << std::endl << Iref << std::endl;

  std::cout << "==== REF_COL I[:,2]=1 =====\n" << std::endl;
  Iref.ref_col(I,2);
  Iref.set();
  std::cout << std::endl << I << std::endl;

  std::cout << "==== REF_COL I[:,20]=0 =====\n" << std::endl;
  Iref.ref_col(I,20);
  Iref.clear();
  std::cout << std::endl << I << std::endl;

  std::cout << "==== REF_ROW I[3,:]=1 =====\n" << std::endl;
  Iref.ref_row(I,3);
  Iref.set();
  std::cout << std::endl << I << std::endl;

  std::cout << "==== REF_ROW I[30,:]=0 =====\n" << std::endl;
  Iref.ref_row(I,30);
  Iref.clear();
  std::cout << std::endl << I << std::endl;

  std::cout << "==== FULL COPY =====\n" << std::endl;
  binary_matrix Icopy(rows,cols);
  Icopy.full_copy(I);
  I.set();
  std::cout << std::endl << Icopy << std::endl;

  std::cout << "==== ADD =====\n" << std::endl;
  binary_matrix Isum(rows,cols);
  add(I,Icopy,Isum);
  std::cout << std::endl << Isum << std::endl;


#if 0

OK  ref(orig);
OK ref_block(orig,i0,i1,j0,j1);
OK  aligned();
OK  clear(); 
OK  set();
OK  weight(); 
OK  row_weight(i);
OK  col_weight(j);
OK  sum(); 
OK  row_sum(i); 
OK  col_sum(idx_t j); 
OK  shallow_copy(orig);
OK  ref_col(orig,j);
OK  ref_row(orig,i);
OK  full_copy(orig);
  // =
  add(A, B, C);
  eland(A, B, C);
  elor(A, B, C);
  elnot(A, B);
  mul(A, At, B, Bt, C);
#endif
  return 0;
}
