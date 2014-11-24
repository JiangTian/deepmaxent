#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <cassert>
#include "common.hpp"
#include "density.hpp"
#include "feature.hpp"
using namespace std;

Density::Density(const vector<vector<real> > & w,
		 const vector<Feature*> & phi,
		 const Dataset & S, int SpSize) //TODO: testing time
  :w(w), phi(phi), normalizerS(0.), normalizerSp(0.),
   lognormalizerS(0.), lognormalizerSp(0.),
   factorsS(S.size()), factorsSp(SpSize),
   expFactorsS(S.size()), expFactorsSp(SpSize), Sp(SpSize) {
  // precompute on S
  precomputeFactors(S, factorsS, expFactorsS, normalizerS, lognormalizerS); // we don't need that!

  // generate Sp
  //  get bounds. TODO: get these bounds once and for all
  int inputSize = S[0].size();
  vector<real> xmin(inputSize, numeric_limits<real>::max());
  vector<real> xmax(inputSize, -numeric_limits<real>::max());
  for (int i = 0; i < S.size(); ++i)
    for (int j = 0; j < S[i].size(); ++j) {
      xmin[j] = min(xmin[j], S[i][j]);
      xmax[j] = max(xmax[j], S[i][j]);
    }
  //  sample Sp between bounds
  for (int i = 0; i < SpSize; ++i) {
    Sp[i] = Datapoint(inputSize);
    for (int j = 0; j < inputSize; ++j)
      Sp[i][j] = uniform(xmin[j], xmax[j]);
  }

  // precompute on Sp
  precomputeFactors(Sp, factorsSp, expFactorsSp, normalizerSp, lognormalizerSp);
}

void Density::precomputeFactors(const Dataset & S, vector<real> & factors,
				vector<real> & expFactors,
				real & normalizer, real & lognormalizer) {
  largestFactor = -numeric_limits<real>::max();
  factors.resize(S.size());
  expFactors.resize(S.size());
  for (int i = 0; i < S.size(); ++i) {
    real t = 0.;
    for (int j = 0; j < phi.size(); ++j) {
      vector<real> phix = phi[j]->eval(S[i]);
      for (int k = 0; k < phix.size(); ++k)
	t += w[j][k] * phix[k];
    }
    factors[i] = t;
    if (t > largestFactor)
      largestFactor = t;
  }
  for (int i = 0; i < S.size(); ++i) { 
    expFactors[i] = exp(factors[i]);
    if (expFactors[i] > numeric_limits<real>::max()) {
      cerr << "The factor is too big for machine precision. "
	   << "Likely that w became too large" << endl;
      exit(0);
    }
  }
  lognormalizer = 0.;
  for (int i = 0; i < S.size(); ++i)
    lognormalizer += exp(factors[i] - largestFactor);
  lognormalizer = log(lognormalizer);
  lognormalizer += largestFactor;
  normalizer = exp(lognormalizer);
}

real Density::evalS(const Datapoint & X) const {
  real t = 0.;
  for (int j = 0; j < phi.size(); ++j) {
    vector<real> phix = phi[j]->eval(X);
    for (int k = 0; k < phix.size(); ++k)
      t += w[j][k] * phix[k];
  }
  return exp(t) / normalizerS; // TODO: this is not a probability 
}

real Density::evalS(int i) const {
  return expFactorsS[i] / normalizerS;
}

real Density::evalSp(const Datapoint & X) const {
  real t = 0.;
  for (int j = 0; j < phi.size(); ++j) {
    vector<real> phix = phi[j]->eval(X);
    for (int k = 0; k < phix.size(); ++k)
      t += w[j][k] * phix[k];
  }
  return exp(t) / normalizerSp; // TODO: this is not a probability 
}

real Density::evalSp(int i) const {
  return expFactorsSp[i] / normalizerSp;
}

real Density::lossS(real beta) const {
  real loss = 0.;
  for (int i = 0; i < factorsS.size(); ++i)
    loss -= factorsS[i];
  loss /= factorsS.size();
  loss += lognormalizerSp;
  for (int i = 0; i < w.size(); ++i) {
    real wnorm1 = 0.;
    for (int j = 0; j < w[i].size(); ++j)
      wnorm1 += abs(w[i][j]);
    loss += (beta + 2 * phi[i]->RademacherComplexity()) * wnorm1;
  }
  return loss;
}
