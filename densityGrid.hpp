#ifndef __DENSITY_GRID_HPP__
#define __DENSITY_GRID_HPP__

#include <vector>
#include "common.hpp"

class Feature;

class DensityGrid {
public:
  std::vector<std::vector<real> > w;
  std::vector<Feature*> phi;
  std::vector<real> factorsS;
  std::vector<std::vector<real> > EphiPWvals;
  real logNormalizer;
  real meanFactorsS;
public:
  DensityGrid() {}
  // note that S should be included in worldGrid
  DensityGrid(const std::vector<std::vector<real> > & w,
	      const std::vector<Feature*> & phi,
	      const Dataset & S, const Dataset & worldGrid);
  real eval(const Datapoint & X) const;
  real loss(real beta) const;
  std::vector<real> EphiPW(int j) const; // expectation of the feature j under pw
};

#endif
