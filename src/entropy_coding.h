#ifndef ENTROPY_CODING_H
#define ENTROPY_CODING_H
/**
 * Functions for computing approximate (ideal) codelenghts for using with
 * MDL-based model selection methods 
 */
#include "binmat.h"
//#include "GolombCoder.h"

#define COSMOS_2E	5.436563656918090181591196596855297684669
#define COSMOS_2PI	6.283185307179586231995926937088370323181
#define COSMOS_2EPI	17.07946844534713193297648103907704353333
#define COSMOS_LOG2E	1.442695040888963387004650940070860087872
#define COSMOS_LOG2PI	1.651496129472318719066947778628673404455
#define COSMOS_LOG2EPI	4.0941911703612818840269937936682254076

/**
 * Approximate codelength of using a two-parts universal enumerative code
 * for describing a vector of length n with r ones on it.
 */
double enumerative_codelength(const unsigned n,
			      const unsigned r);

/**
 * Approximate codelength of using a one-part universal codelength 
 * for describing a vector of length n with r ones on it.
 */
double universal_codelength(const unsigned n,
			      const unsigned r);

/**
 * Overall codelength of describing the decomposition of a matrix
 * X in terms of E,D,A: X = E + DA
 */
idx_t model_codelength(const binary_matrix& E, 
		       const binary_matrix& D, 
		       const binary_matrix& A);
           
#endif
