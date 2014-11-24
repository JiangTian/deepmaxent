#ifndef __COMMON_HPP__
#define __COMMON_HPP__

#include <vector>
#include <cstdlib>

typedef double real;
typedef std::vector<real> Datapoint;
typedef std::vector<Datapoint> Dataset;

inline real uniform(real a, real b) {
  return ((real)rand() / (real)RAND_MAX) * (b - a) + a;
}

#endif
