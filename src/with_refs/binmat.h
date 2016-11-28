#ifndef BINMAT_H
#define BINMAT_H
//#include <ostream>
#include <iostream>
#include <bitset>

typedef unsigned long idx_t;
typedef unsigned long block_t;

#define BITS_PER_BLOCK (sizeof(block_t)*8)
typedef std::bitset<BITS_PER_BLOCK> bm_bitset; // for printing

/** 
 * ROW major (C-style) packed binary matrix 
 * Whenever possible, operations are carried out at block level.
 */
class binary_matrix {
public:

  /** 
   * Default constructor.
   */
 binary_matrix(): rows(), cols(), len(), stride(), first_bit_offset(), last_bit_offset(), data_blocks(), data(0),data_owner(true), head_mask(~block_t(0)), trail_mask(~block_t(0)) {}

  /**
   * Allocates a binary matrix of the specified dimensions.
   */
  binary_matrix(idx_t _rows, idx_t _cols);

  /**
   * Initializes  a binary matrix of the specified dimensions, with external storage.
   */
  binary_matrix(idx_t _rows, idx_t _cols, block_t* _data);

  /** 
   * Destructor
   */
  ~binary_matrix() { if (data_owner) delete[] data; }

  inline binary_matrix& ref(binary_matrix& orig) { 
    shallow_copy(orig);
    return *this;
  }


 /** 
  * Makes this matrix a shallow reference to the j-th column of the original matrix
  * passed as argument.
  */
 binary_matrix& ref_col(binary_matrix& orig, idx_t j);

 /** 
  * Makes this matrix a shallow reference to the i-th row of the original matrix
  * passed as argument.
  */
 binary_matrix& ref_row(binary_matrix& orig, idx_t i);

 /** 
  * Makes this matrix a shallow reference to a given block  of the original matrix
  * passed as argument.
  */
 binary_matrix& ref_block(binary_matrix& orig, idx_t i0, idx_t i1, idx_t j0, idx_t j1);


 /** 
  * Makes this matrix a shallow reference to the original matrix
  * passed as argument.
  */
 void shallow_copy(binary_matrix& orig);

 /** 
  * Makes this matrix a full copy of the original matrix passed as reference
  */
 void full_copy(const binary_matrix& orig);

 /** 
  * Overloaded assignement operator.
  */ 
 inline binary_matrix& operator=(binary_matrix& orig) { shallow_copy(orig); return *this; }

  /**
   * @return true if the first bit starts at the MSB of the first block.
   */
 inline bool aligned() const { return !first_bit_offset; }

 /** @return the number of rows of the matrix */
 inline idx_t get_rows() const { return rows; }

 /** @return the number of columns of the matrix */
 inline idx_t get_cols() const { return cols; }

 /** @return the total number of bits of the matrix */
 inline idx_t get_len() const { return len; }

 /**
  * Sets all bits to zero. FIX: will erase whole block!
  */
 void clear(); 

 /**
  * Sets all bits to one.
  */
 void set();

 /**
  * Flip bits
  */
 void invert();

 /**
  * @return the bit at the specified position
  */
 inline bool get(const idx_t i, const idx_t j) const { 
   return data[i*stride + ((j+first_bit_offset) / BITS_PER_BLOCK)] & (1UL<<(BITS_PER_BLOCK - ((j+first_bit_offset) % BITS_PER_BLOCK) - 1)); 
 }

 /**
  * Sets the bit at the specified position
  */
 inline void set(const idx_t i, const idx_t j, const bool v) {
   if (v) set(i,j); else clear(i,j);
 }

 /**
  * Sets the bit at the specified position to 1
  */
 inline void set(const idx_t i, const idx_t j) {
   data[i*stride + ((j+first_bit_offset) / BITS_PER_BLOCK)] |= (1UL<<(BITS_PER_BLOCK - ((j+first_bit_offset) % BITS_PER_BLOCK) - 1));
 }

 /**
  * Sets the bit at the specified position to 0
  */
 inline void clear(const idx_t i, const idx_t j) {
   data[i*stride + ((j+first_bit_offset) / BITS_PER_BLOCK)] &= ~block_t(0) ^ (1UL<<(BITS_PER_BLOCK - ((j+first_bit_offset) % BITS_PER_BLOCK) - 1)); 
 }


 /** @return the number of ones in the matrix */
 idx_t weight() const; 

 /** @return the number of ones in the given row */
 idx_t row_weight(idx_t i) const; 

 /** @return the number of ones in the given column */
 idx_t col_weight(idx_t j) const; 

 /** @return the boolean sum (XOR), that is, the parity, of the whole matrix */
 bool sum() const; 

 /** @return the boolean sum (XOR), that is, the parity, of a given row */
 bool row_sum(idx_t i) const; 

 /** @return the boolean sum (XOR), that is, the parity, of a given column */
 bool col_sum(idx_t j) const; 

 /** Overloaded print operator */
 friend std::ostream& operator<<(std::ostream& out, const binary_matrix& A);

 /** C = A XOR B. It is assumed that C contains enough space for the result. */
 friend binary_matrix& add(const binary_matrix& A, const binary_matrix& B, binary_matrix& C);

 //#define xor add

 /** C = A AND B. It is assumed that C contains enough space for the result. */
 friend binary_matrix& eland( const binary_matrix& A, const binary_matrix& B, binary_matrix& C);

 /** C = A AOR B. It is assumed that C contains enough space for the result. */
 friend binary_matrix& elor( const binary_matrix& A, const binary_matrix& B, binary_matrix& C);

 /** B = NOT A. It is assumed that B contains enough space for the result. */
 friend binary_matrix& elnot(const binary_matrix& A, binary_matrix& B);

 /** 
  * Linear matrix product, C = AxB (each possibly transposed according to the flags). 
  * It is assumed that C contains enough space for the result. 
  */
 friend binary_matrix& mul(const binary_matrix& A, const bool At, const binary_matrix& B, const bool Bt, binary_matrix& C);

private:

 /** @return the specified data block, aligned case */
 inline block_t get_block_aligned(const idx_t i, const idx_t j) const {
   return (j < stride) ? data[i*stride + j] : data[i*stride + j]; 
 }

 /** @return the specified data block, generic case */
 block_t get_block(const idx_t i, const idx_t j) const;

 /** Sets the specified data block to the given value, aligned case */
 inline void set_block_aligned(const idx_t i, const idx_t j, const block_t b)  { 
   data[i*stride+j] = b;
 }

 /** Sets the specified data block to the given value, generic case */
 void set_block(const idx_t i, const idx_t j, const block_t b);

 /** C=A*B */
 friend binary_matrix& mul_AB(const binary_matrix& A, const binary_matrix& B, binary_matrix& C);

 /** C=A^t*B */
 friend binary_matrix& mul_AtB(const binary_matrix& A, const binary_matrix& B, binary_matrix& C);

 /** C=A*B^t */
 friend binary_matrix& mul_ABt(const binary_matrix& A, const binary_matrix& B, binary_matrix& C);

 /** C=A^t*B^t */
 friend binary_matrix& mul_AtBt(const binary_matrix& A, const binary_matrix& B, binary_matrix& C);

 /** rows of the matrix */
  idx_t rows;

  /** columns of the matrix */
  idx_t cols;

  /** total number of bits. */
  idx_t len;

  /** number of blocks per row */
  idx_t stride;

  /** first bit offset, if 0 this means the matrix is aligned to block boundary, which makes things MUCH easier */
  idx_t first_bit_offset; 

  /** last bit offset, may be non-zero if cols is not a multiple of block size */
  idx_t last_bit_offset; 

  /** number of actual data blocks allocated for storage */
  idx_t data_blocks;

  idx_t blocks_per_row; /** number of blocks required to hold a row; the matrix may actually span one more memory block than blocks_per_row if it is misaligned */

  idx_t spanned_blocks_per_row; /** number of blocks involved in the representation of the bits of one row of this matrix; this is an internal implementation detail. */
  /** allocated memory */
  block_t* data; 

  /** true if this matrix owns the data. False if this matrix refers to external storage or other matrix */
  bool data_owner;

  /** for masking the heading bits of the matrix, in case it is not aligned */
  block_t head_mask; // mask bits before first column, if bit_offset is nonzero

  /** for masking the trailing bits of the matrix */
  block_t trail_mask; // mask trailing bits
};

#endif
