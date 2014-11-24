#ifndef __DENSITY_SAMPLE_HPP__
#define __DENSITY_SAMPLE_HPP__

#include <vector>
#include "common.hpp"

class Feature;

class DensitySample {
public:
  std::vector<std::vector<real> > w;
  std::vector<Feature*> phi;
private:
  std::vector<real> factorsS; // = w_t phi(x) on S
  std::vector<real> factorsSp; // = w_t phi(x) on Sp
  std::vector<real> expFactorsS; // = exp(w_t phi(x) ) on S
  std::vector<real> expFactorsSp; // = exp(w_t phi(x) ) on Sp
  real largestFactor; // largest exp(w.phi(x))
  real normalizerS; // = sum_{x\in S} (exp(w_t phi(x) ) ) on S
  real normalizerSp; // = sum_{x\in S} (exp(w_t phi(x) ) ) on Sp
  real lognormalizerS; // = log( sum_{x\in S} (exp(w_t phi(x) ) ) ) on S
  real lognormalizerSp; // = log( sum_{x\in S} (exp(w_t phi(x) ) ) ) on Sp
  Dataset Sp; // uniformly drawn at construction
  void precomputeFactors(const Dataset & S, std::vector<real> & factors,
			 std::vector<real> & expFactors,
			 real & normalizer, real & lognormalizer);
public:
  DensitySample() {}
  DensitySample(const std::vector<std::vector<real> > & w,
		const std::vector<Feature*> & phi,
		const Dataset & S, int SpSize); //TODO: testing time
  real eval(const Datapoint & X) const; // evaluated on Sp[i] but is faster
  real loss(real beta) const;
  std::vector<real> EphiPW(int j) const; // expectation of the feature j under pw
};

#endif
