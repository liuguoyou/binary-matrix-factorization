#include "intmat.h"

integer_matrix::integer_matrix(idx_t _rows, idx_t _cols):rows(_rows),cols(_cols),len(rows*cols) {
  data = new int_t[len];
}

integer_matrix::integer_matrix(const integer_matrix& other): integer_matrix(other.rows, other.cols) {
  other.copy_to(*this);
}

void integer_matrix::allocate(idx_t _rows, idx_t _cols) {
  if (data) destroy();
  rows = _rows;
  cols = _cols;
  len = rows*cols;
  data = new int_t[len];
}

integer_matrix integer_matrix::get_vectorized() const {
  integer_matrix res(1,rows*cols);
  std::copy(data,data+len,res.data);
  return res;
}

integer_matrix integer_matrix::get_col(const idx_t j) const {
  integer_matrix res(rows,1);
  for (size_t k = j, i=0; k < len; k+= cols, ++i) {
    res.data[i] = this->data[k];
  }
  return res;
}

integer_matrix integer_matrix::get_row(const idx_t i) const {
  integer_matrix res(1,cols);
  std::copy(data+i*cols,data+(i+1)*cols,res.data);
  return res;
}

integer_matrix integer_matrix::get_submatrix(const idx_t i0,const  idx_t i1,const  idx_t j0,const  idx_t j1) const {
  integer_matrix res(i1-i0,j1-j0);
  // PENDING
  return res;
}

integer_matrix integer_matrix::get_copy() const {
  integer_matrix res(rows,cols);
  std::copy(data,data+len,res.data);
  return res;
}

integer_matrix integer_matrix::get_transposed() const {
  integer_matrix res(cols,rows);
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      res.data[j*rows+i] = this->data[i*cols+j];
    }
  }
  return res;
}

void integer_matrix::copy_vectorized_to(integer_matrix& B) const {
  B.cols = 1;
  B.rows = this->rows*this->cols;
  std::copy(data,data+len,B.data);
}

void integer_matrix::copy_col_to(const idx_t j, integer_matrix& B) const {
  // PENDING
}

void integer_matrix::copy_row_to(const idx_t i, integer_matrix& B) const {
  // PENDING
}

void integer_matrix::copy_submatrix_to(const idx_t i0,const  idx_t i1,const  idx_t j0,const  idx_t j1, integer_matrix& B) const {
  // PENDING
}

void integer_matrix::copy_to(integer_matrix& B) const {
  // PENDING
}

void integer_matrix::transpose_to(integer_matrix& B) const {
  // PENDING
}

void integer_matrix::set_vectorized(const integer_matrix& src) {
  // PENDING
}

void integer_matrix::set_col(const idx_t j, const integer_matrix& src) {
  // PENDING
}

void integer_matrix::set_row(const idx_t i, const integer_matrix& src) {
  // PENDING
}

void integer_matrix::set_submatrix(const idx_t i0,const  idx_t j0, const integer_matrix& src) {
  // PENDING
}

void integer_matrix::add_rows(idx_t nrows) {
  // PENDING
}

void integer_matrix::remove_rows(idx_t nrows) {
  // PENDING
}


int_t integer_matrix::sum() const {
  return 0;
}

int_t integer_matrix::row_sum(idx_t i) const {
  // PENDING
  return 0;
}

int_t integer_matrix::col_sum(idx_t j) const {
  // PENDING
  return 0;
}

std::ostream& operator<<(std::ostream& out, const integer_matrix& A) {
  // PENDING
  return out;
}

integer_matrix& add(const integer_matrix& A, const integer_matrix& B, integer_matrix& C) {
  // PENDING
  return C;
}


integer_matrix& mul(const integer_matrix& A, const bool At, const integer_matrix& B, const bool Bt, integer_matrix& C) {
  // PENDING
  return C;
}


idx_t dist(const integer_matrix& A, const integer_matrix& B) {
  // PENDING
  return 0;
}

integer_matrix& integer_matrix::operator=(const integer_matrix& A) {
  this->destroy();
  std::copy(&A,&A+sizeof(A),this);
  return *this;
}

integer_matrix& mul_AB(const integer_matrix& A, const integer_matrix& B, integer_matrix& C) {
  // PENDING
  return C;
}

integer_matrix& mul_AtB(const integer_matrix& A, const integer_matrix& B, integer_matrix& C) {
  // PENDING
  return C;
}

 integer_matrix& mul_ABt(const integer_matrix& A, const integer_matrix& B, integer_matrix& C) {
  // PENDINGF
  return C;
}

 integer_matrix& mul_AtBt(const integer_matrix& A, const integer_matrix& B, integer_matrix& C) {
  // PENDING
  return C;
}

