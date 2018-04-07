#ifndef INTMAT_H
#define INTMAT_H
#include <algorithm>
#include <fstream>
#include <cstring>
#include "base_types.h"
#include "binmat.h"

typedef long int_t;

/** 
 * ROW major (C-style) packed integer matrix 
 * Whenever possible, operations are carried out at block level.
 */
class integer_matrix {
public:

  /**
   * Allocates a integer matrix of the specified dimensions.
   */
  integer_matrix(idx_t _rows, idx_t _cols);

  integer_matrix() { memset(this,0,sizeof(integer_matrix)); }

 /**
  * Deep copy, with allocation
  */
 integer_matrix(const integer_matrix& other);

  /** 
   * Destructor
   */
  ~integer_matrix() {  } // we delete explicitly

 void allocate(idx_t _rows, idx_t _cols);

 integer_matrix get_vectorized() const;

 integer_matrix get_col(const idx_t j) const;

 integer_matrix get_row(const idx_t i) const ;

 integer_matrix get_submatrix(const idx_t i0,const  idx_t i1,const  idx_t j0,const  idx_t j1) const;

 integer_matrix get_copy() const;

 integer_matrix get_transposed() const;

 void copy_vectorized_to(integer_matrix& B) const;

 void copy_col_to(const idx_t j, integer_matrix& B) const;

 void copy_row_to(const idx_t i, integer_matrix& B) const ;

 void copy_submatrix_to(const idx_t i0,const  idx_t i1,const  idx_t j0,const  idx_t j1, integer_matrix& B) const;

 void copy_to(integer_matrix& B) const;

 void transpose_to(integer_matrix& B) const;

 void set_vectorized(const integer_matrix& src);

 void set_col(const idx_t j, const integer_matrix& src);

 void set_row(const idx_t i, const integer_matrix& src);

 void set_submatrix(const idx_t i0,const  idx_t j0, const integer_matrix& src);

 void add_rows(idx_t nrows);

 void remove_rows(idx_t nrows); // lightweight, does not deallocate space!

 /** @return the number of rows of the matrix */
 inline idx_t get_rows() const { return rows; }

 /** @return the number of columns of the matrix */
 inline idx_t get_cols() const { return cols; }

 /** @return the total number of bits of the matrix */
 inline idx_t get_len() const { return len; }

 inline bool empty() const { return get_len() == 0; }
 /**
  * Sets all  to zero. 
  */
 inline void clear() { std::fill(data,data+len,0); }


 inline int_t get(const idx_t i, const idx_t j) const { 
   return data[i*cols+j];
 }

 inline void set(const idx_t i, const idx_t j, const int_t v) {
   data[i*cols+j] = v; 
 }

 inline void inc(const idx_t i, const idx_t j) {
   data[i*cols+j]++;
 }

 inline void dec(const idx_t i, const idx_t j) {
   data[i*cols+j]--;
 }

 int_t sum() const; 

 int_t row_sum(idx_t i) const; 

 int_t  col_sum(idx_t j) const; 

 /** Overloaded print operator */
 friend std::ostream& operator<<(std::ostream& out, const integer_matrix& A);

 /** C = A + B. It is assumed that C contains enough space for the result. */
 friend integer_matrix& add(const integer_matrix& A, const integer_matrix& B, integer_matrix& C);


 /** 
  * Linear matrix product, C = AxB (each possibly transposed according to the flags). 
  * It is assumed that C contains enough space for the result. 
  */
 friend integer_matrix& mul(const integer_matrix& A, const bool At, const integer_matrix& B, const bool Bt, integer_matrix& C);
 
 friend integer_matrix& mul(const binary_matrix& A, const bool At, const binary_matrix& B, const bool Bt, integer_matrix& C);

 void destroy() { delete[] data; memset(this,0,sizeof(integer_matrix)); }

 friend idx_t dist(const integer_matrix& A, const integer_matrix& B);

 integer_matrix& operator=(const integer_matrix& A);

 /** C=A*B */
 friend integer_matrix& mul_AB(const integer_matrix& A, const integer_matrix& B, integer_matrix& C);

 /** C=A^t*B */
 friend integer_matrix& mul_AtB(const integer_matrix& A, const integer_matrix& B, integer_matrix& C);

 /** C=A*B^t */
 friend integer_matrix& mul_ABt(const integer_matrix& A, const integer_matrix& B, integer_matrix& C);

 /** C=A^t*B^t */
 friend integer_matrix& mul_AtBt(const integer_matrix& A, const integer_matrix& B, integer_matrix& C);

 /** C=A*B */
 friend integer_matrix& mul_AB(const binary_matrix& A, const binary_matrix& B, integer_matrix& C);

 /** C=A^t*B */
 friend integer_matrix& mul_AtB(const binary_matrix& A, const binary_matrix& B, integer_matrix& C);

 /** C=A*B^t */
 friend integer_matrix& mul_ABt(const binary_matrix& A, const binary_matrix& B, integer_matrix& C);

 /** C=A^t*B^t */
 friend integer_matrix& mul_AtBt(const binary_matrix& A, const binary_matrix& B, integer_matrix& C);

private:
 
 /** rows of the matrix */
  idx_t rows;

  /** columns of the matrix */
  idx_t cols;

  /** total number of bits. */
  idx_t len;

  int_t* data; 
};

void set_grid_width(idx_t g);

#endif
