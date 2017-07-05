#include "gsl/gsl_randist.h"

long random_seed = 34503498;


gsl_rng* get_rng() {
  static gsl_rng* rng = 0;
  if (rng == 0) { 
    rng = gsl_rng_alloc (gsl_rng_rand48);
    gsl_rng_set (rng, random_seed); // set random_seed
  } 
  return rng;
}

unsigned long get_uniform_unsigned_sample(unsigned maxval) {
  return (unsigned long) gsl_rng_uniform_int(get_rng(),maxval);  
}

bool get_bernoulli_sample(double p) {
  return gsl_ran_bernoulli (get_rng(),p);
}
