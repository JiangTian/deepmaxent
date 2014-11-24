#include <iostream>
#include <cmath>
#include <cstdio>
#include <vector>
#include <cstdlib>
#include <limits>
#include <cfloat>
#include <cassert>
#include "deepmaxent.hpp"
#include "data.hpp"
#include "plot.hpp"
using namespace std;

//#define USE_IRIS
#define USE_DUMMY

real uniform() {
  return (real)rand() / (real)RAND_MAX;
}

real normal() {
  real x = 0.;
  int n = 1000;
  for (int i = 1; i < n; ++i)
    x += 2. * uniform() - 1.;
  return x / sqrt(n);
}

enum dataset_type {
  DATASET_DUMMY_UNIFORM,
  DATASET_DUMMY_GAUSSIAN,
  DATASET_IRIS,
  DATASET_ELNINO,
};

Dataset getDataset(dataset_type datatype) {
  Dataset S;
  switch (datatype) {
  case(DATASET_DUMMY_UNIFORM):
    {
      int nDatapoints = 400;
      int inputDim = 2;
      for (int i = 0; i < nDatapoints; ++i) {
	Datapoint p(inputDim);
	for (int j = 0; j < inputDim; ++j) {
	  p[j] = uniform();
	}
	S.push_back(p);
      }
      S[0][0] = -1;
      S[0][1] = -1;
      break;
    }
  case(DATASET_DUMMY_GAUSSIAN):
    {
      int nDatapoints = 400;
      int inputDim = 2;
      for (int i = 0; i < nDatapoints; ++i) {
	Datapoint p(inputDim);
	for (int j = 0; j < inputDim; ++j) {
	  p[j] = normal() + 1;
	  if (p[j] > 2)
	    p[j] = 2;
	  if (p[j] < -2)
	    p[j] = -2;
	}
	S.push_back(p);
      }
      break;
    }
  case(DATASET_IRIS):
    {
      S = readIris("data/iris.data");
      break;
    }
  case(DATASET_ELNINO):
    {
      Dataset S = readElNino("data/tao-all2.dat");
      break;
    }
  }
  return S;
}

Density* pointer_to_pw = NULL;

real eval_pw_at_point_2(real x, real y) {
  Datapoint p(2);
  p[0] = x;
  p[1] = y;
  return pointer_to_pw->eval(p);
}

Dataset* pointer_to_S;
real eval_raw_iris(real x, real y) {
  real n = 0;
  real eps = 0.1;
  Dataset & S = *pointer_to_S;
  for (int i = 0; i < S.size(); i++) {
    if ((0 < S[i][2] - x) && (S[i][2] - x < eps) && (0 < S[i][3] - y) && (S[i][3] - y < eps))
      ++n;
  }
  return n;
}

real eval_pw_at_point_iris(real x, real y) {
  Datapoint p(4);
  p[2] = x;
  p[3] = y;
 
  real sum = 0.;
  for (real u = 0; u < 10; u += 0.2) {
    p[0] = u;
    for (real v = 0; v < 10; v += 0.2) {
      p[1] = v;
      sum += pointer_to_pw->eval(p);
    }
  }
  return sum;
 
  //p[2] = 4.;
  //p[3] = 2.;
  //return pointer_to_pw->eval(p);
}

int main(int argc, char* argv[]) {

  dataset_type datatype = DATASET_IRIS;
  Dataset S = getDataset(datatype);
  pointer_to_S = &S;
  
  // Normalize data to be between -1 and 1
  {
    int inputSize = S[0].size();
    vector<real> xmin(inputSize, numeric_limits<real>::max());
    vector<real> xmax(inputSize, -numeric_limits<real>::max());
    for (int i = 0; i < S.size(); ++i)
      for (int j = 0; j < S[i].size(); ++j) {
	xmin[j] = min(xmin[j], S[i][j]);
	xmax[j] = max(xmax[j], S[i][j]);
      }
    for (int i = 0; i < S.size(); ++i)
      for (int j = 0; j < S[i].size(); ++j)
	S[i][j] = 2. * (S[i][j] - xmin[j]) / (xmax[j] - xmin[j]) - 1;
  }

  //gridPlot(&eval_raw_iris, 3.5, 8.5, 0.5, 5.5);
  gridPlot(&eval_raw_iris, -1, 1, -1, 1, 0);
  
  Density pw = DeepMaxent(S, 1000, 8);
  pointer_to_pw = &pw;
  
  cout << "W=" << endl;
  for (int i = 0; i < pw.w.size(); ++i)
    cout << pw.w[i][0] << " ";
  cout << endl;
  
  switch(datatype) {
  case(DATASET_DUMMY_UNIFORM):
    gridPlot(&eval_pw_at_point_2, -1, 1, -1, 1, 1);
    break;
  case(DATASET_DUMMY_GAUSSIAN):
    gridPlot(&eval_pw_at_point_2, -1, 1, -1, 1, 1);
    break;
  case(DATASET_IRIS):
    gridPlot(&eval_pw_at_point_iris, -1, 1, -1, 1, 1);
    break;
  case(DATASET_ELNINO):
    assert(0); //Not implemented yet
    break;
  }
  
  
  /*
  Datapoint mean(inputDim);
  Datapoint peri(inputDim);
  for (int i = 0; i < inputDim; ++i) {
    mean[i] = 0.;
    peri[i] = 1.;
  }
  cout << "P(mode) = " << pw.evalS(mean);
  cout << " P(peri) = " << pw.evalS(peri) << endl;
 
#if 1
  for (int i = 0; i < pw.w.size(); ++i)
    for (int j = 0; j < pw.w[i].size(); ++j)
      cout << pw.w[i][j] << " ";
  cout << endl;
#endif
  */

  return 0;
}
