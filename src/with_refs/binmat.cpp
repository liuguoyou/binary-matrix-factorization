#include "binmat.h"
#include <iomanip>
#include <cstring>
#include <bitset>
#include <cassert>
//#include <omp.h>

/* information on how to vectorize these operations using x86 SSE extensions in C is
   described in https://gcc.gnu.org/onlinedocs/gcc/Vector-Extensions.html */

/* https://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetTable */

/*
 * very efficient general purpose method for counting bits w/o requiring
 * a look up table. An alternative is a LUT as a 32 bit block where one has 16 2 bit
 * counts for every of the 16 possible 4 bit nibbles.
 * 
 * There is a POPCNT CPU extension in intel-base 64 bit processors since 2014,
 * but I don't want to be compiler dependent for now.
 */
static idx_t block_weight(const block_t& v)
{
  static const unsigned char byte_weight_table[256] = 
{
#   define B2(n) n,     n+1,     n+1,     n+2
#   define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#   define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
    B6(0), B6(1), B6(1), B6(2)
};
  const unsigned char * p = (const unsigned char *) &v;
  register idx_t w = 0;
  for (idx_t i = 0; i < sizeof(block_t); ++i) {
    w += byte_weight_table[p[i]];
  }
  return w;
}

/* parity computation, much faster than sum */
static bool block_sum(const block_t& v) {
  static const bool byte_parity_table[256] = {
#   define P2(n) n, n^1, n^1, n
#   define P4(n) P2(n), P2(n^1), P2(n^1), P2(n)
#   define P6(n) P4(n), P4(n^1), P4(n^1), P4(n)
    P6(0), P6(1), P6(1), P6(0)
  };
  const unsigned char * p = (unsigned char *) &v;
  idx_t ti = 0;
#pragma omp parallel for
  for (idx_t i = 0; i < sizeof(block_t); ++i) {
    ti ^= p[i];
  }
  return byte_parity_table[ti];
  //  return ParityTable256[tip[0] ^ p[1] ^ p[2] ^ p[3]];
}

idx_t binary_matrix::weight() const {
  if (rows*cols == 0)
    return 0;
  idx_t w = 0;
  const idx_t Nb = blocks_per_row;  
  for (idx_t i = 0; i < rows; ++i) {
    for (idx_t j = 0; j < Nb ; ++j) {
      w += block_weight(get_block(i,j));
    }
  }
  return w;
}

idx_t binary_matrix::row_weight(idx_t i) const {
  if (rows == 0) return 0;
  idx_t w = 0;
  for (idx_t j = 0; j < (stride-1); ++j) {
    w += block_weight(get_block(i,j));
  }
  return w;
}

idx_t binary_matrix::col_weight(idx_t j) const {
  block_t* pd = &data[(j+first_bit_offset) / BITS_PER_BLOCK];
  block_t mask = 1UL << (BITS_PER_BLOCK - ((first_bit_offset+j) % BITS_PER_BLOCK) - 1);
  idx_t w = 0;
  for (idx_t i = 0; i < data_blocks; i+=stride) {
    if (pd[i] & mask)
      w++;
  }
  return w;
}

bool binary_matrix::sum() const {
  if (rows*cols == 0)
    return 0;
  idx_t w = 0;
  for (idx_t i = 0; i < rows; ++i) {
    for (idx_t j = 0; j < blocks_per_row; ++j) {
      w ^= block_sum(get_block(i,j));
    }
  }
  return w;
}


bool binary_matrix::row_sum(idx_t i) const {
  if (rows*cols == 0)
    return 0;
  idx_t w = 0;
  for (idx_t j = 0; j < blocks_per_row; ++j) {
    w ^= block_sum(get_block(i,j));
  } 
  return w;
}

bool binary_matrix::col_sum(idx_t j) const {
  block_t* pd = &data[(j+first_bit_offset) / BITS_PER_BLOCK];
  block_t bitoff = (BITS_PER_BLOCK - (j % BITS_PER_BLOCK) - first_bit_offset - 1);
  idx_t w = 0;
  for (idx_t i = 0; i < data_blocks; i+=stride) {
    w ^= (pd[i] >> bitoff) & 1UL;
  }
  return w;
}

binary_matrix::binary_matrix(idx_t _rows, idx_t _cols)
  :rows(_rows), cols(_cols), len(_rows*_cols)
{ 
  stride = (cols+BITS_PER_BLOCK-1)/BITS_PER_BLOCK; // ceil
  data_blocks = stride*rows;
  data = new block_t[data_blocks];
  first_bit_offset = 0;
  last_bit_offset = cols-(stride-1)*BITS_PER_BLOCK;
  trail_mask = ~block_t(0) << (stride*BITS_PER_BLOCK-cols);
  if (stride == 1) {  // only one block, so both masks coincide
    head_mask = trail_mask;
  } else {
    head_mask = ~block_t(0);
  }
  blocks_per_row = stride; // these always coincide for data owning matrices
}

binary_matrix::binary_matrix(idx_t _rows, idx_t _cols, block_t* _data)
  :rows(_rows), cols(_cols), len(_rows*_cols)
{ 
  stride = (cols+BITS_PER_BLOCK-1)/BITS_PER_BLOCK; // ceil
  data_blocks = stride*rows;
  data = _data;
  data_owner = false;
  first_bit_offset = 0;
  last_bit_offset = cols-(stride-1)*BITS_PER_BLOCK;
  trail_mask = ~block_t(0) << (stride*BITS_PER_BLOCK-cols);
  if (stride == 1) { // only one block per column, so both masks coincide
    head_mask = trail_mask;
  } else {
    head_mask = ~block_t(0);
  }
  blocks_per_row = stride; 
}

void  binary_matrix::shallow_copy(binary_matrix& orig)  {
  if (data_owner && data) {
    delete[] data;
  }
  rows = orig.rows;
  cols = orig.cols;
  stride = orig.stride; 
  data_blocks = orig.data_blocks;
  data = orig.data;
  first_bit_offset = orig.first_bit_offset;
  last_bit_offset  = orig.last_bit_offset;
  data_owner = false;  
  head_mask = orig.head_mask;
  trail_mask = orig.trail_mask;
  blocks_per_row = orig.blocks_per_row;
}

void binary_matrix::full_copy(const binary_matrix& orig) {
  rows = orig.rows;
  cols = orig.cols;
  assert(data_owner); // we must own this matrix for this method to run
  //
  // if original matrix has non-zero bit offset, 
  // the copied one may have a smaller (by one) stride 
  //
  first_bit_offset = 0;
  if (orig.first_bit_offset > orig.last_bit_offset) { 
    stride = orig.stride-1;
    last_bit_offset = orig.last_bit_offset + BITS_PER_BLOCK - orig.first_bit_offset;
  } else {
    stride = orig.stride;
  }
  //
  // in all cases, we avoid allocating new memory if our matrix
  // has enough room for the copied one
  //
  if (stride*rows > this->data_blocks) {
    this->data_blocks = stride*rows;
    delete[] this->data;
    this->data = new block_t[this->data_blocks];
  } else {
    this->data_blocks = stride*rows;
    if (!this->data) {
      this->data = new block_t[this->data_blocks];
    }
  }
  if (orig.aligned()) {
    // 
    // easy 
    //
    memcpy(data,orig.data,sizeof(block_t)*data_blocks);
    head_mask = orig.head_mask;
    trail_mask = orig.trail_mask;
  } else {
    //
    // not-so-easy
    //
    for (idx_t i = 0; i < rows; ++i) {
      for (idx_t j = 0; j < stride; j++) {
	block_t b = orig.get_block(i,j);
	set_block_aligned(i,j,b);
      }
    }
  }
  blocks_per_row = stride;
  data_owner = true;
}

void binary_matrix::clear() { 
  if (rows*cols == 0) return;
  if (aligned() && (data_owner || (trail_mask == ~block_t(0)))) { 
    memset(data,0,sizeof(block_t)*rows*stride); 
  } else {
    const idx_t Nb = blocks_per_row;
    for (idx_t i = 0; i < rows; ++i) {
      for (idx_t j = 0; j < Nb; ++j) {
	set_block(i,j,block_t(0));
      }
    }
  }
}

void binary_matrix::set() {
  if (rows*cols == 0) return;
  if (aligned() && (data_owner || (trail_mask == ~block_t(0)))) { 
    memset(data,0xff,sizeof(block_t)*rows*stride); 
  } else {
    const idx_t Nb = blocks_per_row;
    for (idx_t i = 0; i < rows; ++i) {
      for (idx_t j = 0; j < Nb; ++j) {
	set_block(i,j,~block_t(0));
      }
    }
  }
}

void binary_matrix::invert() {
  if (rows*cols == 0) return;
  if (aligned() && (data_owner || (trail_mask == ~block_t(0)))) { 
    memset(data,0xff,sizeof(block_t)*rows*stride); 
  } else {
    const idx_t Nb = blocks_per_row;
    for (idx_t i = 0; i < rows; ++i) {
      for (idx_t j = 0; j < Nb; ++j) {
	set_block(i,j,get_block(i,j) ^ ~block_t(0));
      }
    }
  }
}


block_t  binary_matrix::get_block(const idx_t i, const idx_t j) const { 
   assert(i >=0);
   assert(i < rows);
   assert(j >=0);
   assert(j < blocks_per_row);
   if (aligned())  return get_block_aligned(i,j);
   if (j < (blocks_per_row-1)) 
     {
       return ((data[i*stride + j] << first_bit_offset) | (data[i*stride + j + 1] >> (BITS_PER_BLOCK - first_bit_offset)));
     } 
   else if (first_bit_offset + last_bit_offset <= BITS_PER_BLOCK) 
     {
       return ((data[i*stride + j] << first_bit_offset) | (data[i*stride + j + 1] >> (BITS_PER_BLOCK - first_bit_offset))) & trail_mask;
     } 
   else 
     { // only in this case, we don't fetch next word
       return  (data[i*stride + j] << first_bit_offset) & trail_mask;
     }
}

void binary_matrix::set_block(const idx_t i, const idx_t j, const block_t b)  {
  if (aligned()) { set_block_aligned(i,j,b); return; }
  data[i*stride+j] &= ~head_mask;  // erase first bits
  data[i*stride+j] |= (b >> first_bit_offset); // replace by corresp. bits in b

  if ((j < (blocks_per_row-1)) || ((first_bit_offset + last_bit_offset) > BITS_PER_BLOCK)) {
    data[i*stride+j+1] &= ~trail_mask;
    data[i*stride+j+1] |= (b << (BITS_PER_BLOCK-first_bit_offset));    
  }
}

binary_matrix& binary_matrix::ref_col(binary_matrix& orig, idx_t j) {
  assert(0 <= j);
  assert(j <= orig.cols);
  if (data_owner && data) {
    delete[] data;
  }
  // BUG! WILL NOT WORK IF ORIGINAL IS ALREADY MISALIGNED!
  rows = orig.rows;
  cols = 1;
  len = rows;
  stride = orig.stride;
  first_bit_offset = j % BITS_PER_BLOCK;
  last_bit_offset = first_bit_offset;
  data_blocks = orig.data_blocks;
  data = orig.data + (j / BITS_PER_BLOCK);
  data_owner = false;
  head_mask = block_t(1) << (BITS_PER_BLOCK-first_bit_offset);
  trail_mask = head_mask;
  return *this;
}

binary_matrix& binary_matrix::ref_row(binary_matrix& orig, idx_t i) {
  assert(0 <= i);
  assert(i <= orig.rows);
  if (data_owner && data) {
    delete[] data;
  }  
  rows = 1;
  cols = orig.cols;
  len = cols;
  stride = orig.stride;
  first_bit_offset = orig.first_bit_offset;
  last_bit_offset = orig.last_bit_offset;
  data_blocks = orig.data_blocks;
  data = orig.data + (i*stride);
  data_owner = false;
  head_mask = orig.head_mask;
  trail_mask = orig.trail_mask;
  return *this;
}

binary_matrix& binary_matrix::ref_block(binary_matrix& orig, idx_t i0, idx_t i1, idx_t j0, idx_t j1) {
  assert(i0 < i1);
  assert(j0 < j1);
  assert(0 <= i0);
  assert(0 <= j0);
  assert(j1 <= orig.cols);
  assert(i1 <= orig.rows);

  if (data_owner && data) {
    delete[] data;
  }
  data_owner = false;
  stride = orig.stride;
  rows = i1 - i0;
  cols = j1 - j0;
  len = rows*cols;
  data = orig.data + stride*i0 + (j0+orig.first_bit_offset) / BITS_PER_BLOCK;
  first_bit_offset = (j0+orig.first_bit_offset) % BITS_PER_BLOCK;
  last_bit_offset = (j1+orig.first_bit_offset) % BITS_PER_BLOCK;
  blocks_per_row = (BITS_PER_BLOCK-1 + cols) / BITS_PER_BLOCK;
  spanned_blocks_per_row = (BITS_PER_BLOCK-1 + first_bit_offset + cols) / BITS_PER_BLOCK;
  if (blocks_per_row > 1) {
    head_mask = ~block_t(0) >> first_bit_offset;
    trail_mask = ~block_t(0) << (BITS_PER_BLOCK-last_bit_offset);
  } else {
    trail_mask = head_mask = (~block_t(0) >> first_bit_offset) & (~block_t(0) << (BITS_PER_BLOCK-last_bit_offset));
  }
  data_blocks = orig.data_blocks; // not sure what to do with this
  return *this;
}


// slow: I don't think anyone cares about fast dumping, since most of the time
// will be I/O anyway.
std::ostream& operator<<(std::ostream& out, const binary_matrix& A)  {
  out << "rows=" << A.rows << "\tcols=" << A.cols << "\tlen=" << A.len 
      << "\tbpw="<< BITS_PER_BLOCK 
      << "\tdw=" << A.data_blocks
      << "\twpr=" << A.blocks_per_row
      << "\towner=" << A.data_owner << '\n'
      << "fbo=" << A.first_bit_offset
      << "\tlbo=" << A.last_bit_offset << '\n' 
      << "hm=" << bm_bitset(A.head_mask) << '\n'
      << "tm=" << bm_bitset(A.trail_mask);
  out << std::endl;
  out << "        ";
  for (idx_t j = 0; j < A.cols; ++j)
    out << ((j % BITS_PER_BLOCK) ? ' ' : '|') << ' ';
  out << std::endl;
  out << "        ";
  for (idx_t j = 0; j < A.cols; ++j)
    out << ((j % 10) ? '-' : '+') << '-';
  out << std::endl;
  for (idx_t i = 0; i < A.rows; ++i) {
    out << std::setw(5)<< i << ((i % 10) ? " | " : " + ");
    for (idx_t j = 0; j < A.cols; ++j) {
      out << std::setw(1) << A.get(i,j) << ' ';
    }
    out << std::endl;
  }
  return out;
}

// C is assumed to have been allocated and have the appropriate dimension
binary_matrix& add(const binary_matrix& A, const binary_matrix& B, binary_matrix& C)
{
  assert(C.data != 0);
  assert(C.rows == A.rows);
  assert(C.rows == B.rows);
  assert(C.cols == A.cols);
  assert(C.cols == B.cols);
  const idx_t M = A.rows;
  const idx_t Nb = C.blocks_per_row;
  for (idx_t i = 0; i < M; ++i) {      
    for (idx_t j = 0; j < Nb; ++j) {
      C.set_block(i,j, A.get_block(i,j) ^ B.get_block(i,j));
    }
  } 
  return C;
}



// C is assumed to have been allocated and have the appropriate dimension
binary_matrix& mul_AB(const binary_matrix& A, const binary_matrix& B, binary_matrix& C)
{
  assert(C.data != 0);
  assert(C.rows == A.rows);
  assert(C.cols == B.cols);
  assert(A.cols == B.rows);
  const idx_t M = A.rows;
  const idx_t Ns = B.stride;
  const idx_t K = B.rows;
  C.clear();
  idx_t mask = 1UL << (BITS_PER_BLOCK-1);
  idx_t koff = 0;
  for (idx_t k = 0; k < K; ++k) {
    for (idx_t i = 0; i < M; ++i) {
      if (A.get_block(i,koff) & mask) { // if nonzero, the i-th row of C is updated by adding the k-th row of B
	for (idx_t j = 0; j < Ns; ++j) {
	  C.set_block(i,j, C.get_block(i,j) ^ B.get_block(k,j));
	}
      } 
    }
    mask >>= 1;
    if (mask == 0) {
      mask = 1UL << (BITS_PER_BLOCK-1);
      koff++; 
    }
  }
  return C;
}

binary_matrix& mul_AtB(const binary_matrix& A, const binary_matrix& B, binary_matrix& C)
{
  assert(C.data != 0);
  assert(C.rows == A.cols);
  assert(C.cols == B.cols);
  assert(A.rows == B.rows);
  const idx_t K = A.rows;
  const idx_t M = A.cols;
  const idx_t Ns = B.stride;
  C.clear();
  for (idx_t k = 0; k < K; ++k) { 
    idx_t mask = 1UL << (BITS_PER_BLOCK-1);
    idx_t ioff = 0;
    for (idx_t i = 0; i < M; ++i) {      
      //      std::cout << "i=" << i << "mask=" << std::bitset<BITS_PER_BLOCK>(mask) << " C=" << C << std::endl;   
      if (A.get_block(k,ioff) & mask) { // if nonzero, the i-th row of C is updated by adding the k-th row of B
	for (idx_t j = 0; j < Ns; ++j) {
	  C.set_block(i,j, C.get_block(i,j) ^ B.get_block(k,j));
	}
      } 
      mask >>= 1;
      if (mask == 0) {
	mask = 1UL << (BITS_PER_BLOCK-1);
	ioff++;
      }
    } // i
  } // k
  return C;
}

binary_matrix& mul_ABt(const binary_matrix& A, const binary_matrix& B, binary_matrix& C)
{
  assert(C.data != 0);
  assert(C.rows == A.rows);
  assert(C.cols == B.rows);
  assert(A.cols == B.cols);
  const idx_t M = A.rows;
  const idx_t N = B.cols;
  const idx_t K = A.stride;
  for (idx_t i = 0; i < M; ++i) {
    for (idx_t j = 0; j < N; ++j) {
      bool a = 0;
      for (idx_t k = 0; k < K; ++k) {
	a ^= block_sum(A.get_block(i,k) & B.get_block(j,k)); // pA[k] & pB[k]);
      }
      C.set(i,j,a);
    }
  }
  return C;
}

binary_matrix& mul_AtBt(const binary_matrix& A, const binary_matrix& B, binary_matrix& C)
{
  assert(C.data != 0);
  assert(C.rows == A.cols);
  assert(C.cols == B.rows);
  assert(A.rows == B.cols);
  // FALTA!
  return C;
}

binary_matrix& mul(const binary_matrix& A, const bool At, const binary_matrix& B, const bool Bt, binary_matrix& C) {
  if (At && Bt) {
    return mul_AtBt(A,B,C);
  } else if (At && !Bt) {
    return mul_AtB(A,B,C);
  } else if (!At && Bt) {
    return mul_ABt(A,B,C);
  } else {
    return mul_AB(A,B,C);
  }
}

