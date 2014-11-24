#include "densityGrid.hpp"
#include <iostream>
#include <cmath>
#include <limits>
#include "common.hpp"
#include "feature.hpp"
using namespace std;

DensityGrid::DensityGrid(const vector<vector<real> > & w,
			 const vector<Feature*> & phi,
			 const Dataset & S, const Dataset & worldGrid)
  :w(w), phi(phi), logNormalizer(0.), factorsS(S.size()),
   meanFactorsS(0.), EphiPWvals(phi.size()) {
  // compute the normalizer
  real largestFactor = -numeric_limits<real>::max();
  vector<real> factors(worldGrid.size());
  for (int i = 0; i < worldGrid.size(); ++i) {
    real t = 0.;
    for (int j = 0; j < phi.size(); ++j) {
      vector<real> phix = phi[j]->eval(worldGrid[i]);
      for (int k = 0; k < phix.size(); ++k)
	t += w[j][k] * phix[k];
    }
    factors[i] = t;
    if (t > largestFactor)
      largestFactor = t;
  }
  real normalizer = 0.;
  for (int i = 0; i < worldGrid.size(); ++i)
    normalizer += exp(factors[i] - largestFactor);
  logNormalizer = log(normalizer) + largestFactor;

  // compute EphiPWvals
  for (int j = 0; j < phi.size(); ++j)
    EphiPWvals[j] = vector<real>(phi[j]->size());
  for (int i = 0; i < worldGrid.size(); ++i) {
    real pwsi = exp(factors[i] - logNormalizer);
    for (int j = 0; j < phi.size(); ++j) {
      vector<real> phix = phi[j]->eval(worldGrid[i]);
      for (int k = 0; k < phi[j]->size(); ++k)
	EphiPWvals[j][k] += pwsi * phix[k];
    }
  }
  
  // precompute the mean of factors on S (for the loss)
  for (int i = 0; i < S.size(); ++i) {
    real t = 0.;
    for (int j = 0; j < phi.size(); ++j) {
      vector<real> phix = phi[j]->eval(S[i]);
      for (int k = 0; k < phix.size(); ++k)
	meanFactorsS += w[j][k] * phix[k];
    }
  }
  meanFactorsS /= S.size();
}

real DensityGrid::eval(const Datapoint & X) const {
  real t = 0.;
  for (int j = 0; j < phi.size(); ++j) {
    vector<real> phix = phi[j]->eval(X);
    for (int k = 0; k < phix.size(); ++k)
      t += w[j][k] * phix[k];
  }
  return exp(t - logNormalizer);
}

real DensityGrid::loss(real beta) const {
  real loss = -meanFactorsS;
  loss += logNormalizer;
  for (int i = 0; i < w.size(); ++i) {
    real wnorm1 = 0.;
    for (int j = 0; j < w[i].size(); ++j)
      wnorm1 += abs(w[i][j]);
    loss += (beta + 2 * phi[i]->RademacherComplexity()) * wnorm1;
  }
  return loss;
}

vector<real> DensityGrid::EphiPW(int j) const {
  return EphiPWvals[j];
}
