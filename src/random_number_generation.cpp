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
