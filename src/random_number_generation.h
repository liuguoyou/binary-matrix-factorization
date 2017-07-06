#ifndef RANDOM_NUMBER_GENERATION
#define RANDOM_NUMBER_GENERATION

#include "gsl/gsl_randist.h"

/**
 * random seed used
 */
extern long random_seed;

/**
 * obtains the instance of the random number generator
 */
gsl_rng* get_rng();

/**
 * destroys the RNG
 */
void destroy_rng();

/**
 * obtain an integer sample uniformly distributed between 0 and maxval
 */
unsigned long get_uniform_unsigned_sample(unsigned maxval);

/**
 * @return a Bernoulli-distributed pseudo-random number x with probability P(x=1)=p
 */
bool get_bernoulli_sample(double p);

#endif
