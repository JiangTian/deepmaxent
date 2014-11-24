#ifndef __DENSITY_HPP__
#define __DENSITY_HPP__

//#define USE_DENSITY_SAMPLE
#define USE_DENSITY_GRID

#ifdef USE_DENSITY_SAMPLE
#include "densitySample.hpp"
typedef DensitySample Density;
#endif

#ifdef USE_DENSITY_GRID
#include "densityGrid.hpp"
typedef DensityGrid Density;
#endif

#endif
