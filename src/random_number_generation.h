#ifndef RANDOM_NUMBER_GENERATION
#define RANDOM_NUMBER_GENERATION

#include "gsl/gsl_randist.h"

extern long random_seed;

gsl_rng* get_rng();

#endif
