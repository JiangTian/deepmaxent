#include <iostream>
#include <cstdio>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <cassert>
#include "common.hpp"
#include "data.hpp"
using namespace std;

Dataset readElNino(const string & filepath) {
  FILE* file = fopen(filepath.c_str(), "r");

  Dataset dataset;
  const int N_MAX = 32;
  char c[12][N_MAX];
  
  int n_missing = 0;

  while(fscanf(file, "%s %s %s %s %s %s %s %s %s %s %s %s",
	       c[0], c[1], c[2], c[3], c[4], c[5], c[6],
	       c[7], c[8], c[9], c[10], c[11]) == 12) {
    Datapoint p(9);
    int k = 0;
    bool has_missing = false;
    for (int i = 0; i < 12; ++i) {
      if ((i == 0) || (i == 3) || (i == 4))
	continue;
      if (strcmp(c[i], ".") != 0) {
	p[k++] = atof(c[i]);
      } else {
	//cout << "Missing point " << dataset.size() << " " << i << endl;
	has_missing = true;
      }
    }
    if (has_missing)
      ++n_missing;
    else
      dataset.push_back(p);
  }

  int inputSize = dataset[0].size();
  for (int i = 0; i < inputSize; ++i) {
    real xmin = dataset[0][i];
    real xmax = dataset[0][i];
    for (int j = 0; j < dataset.size(); ++j) {
      xmin = min(xmin, dataset[j][i]);
      xmax = max(xmax, dataset[j][i]);
    }
    assert(xmin != xmax);
    for (int j = 0; j < dataset.size(); ++j)
      dataset[j][i] = 2. * (dataset[j][i] - xmin) / (xmax - xmin) - 1;
  }
  
  cout << n_missing << " missing over " << (n_missing+dataset.size()) << endl;

  fclose(file);
  return dataset;
}

Dataset readIris(const string & filepath) {
  FILE* file = fopen(filepath.c_str(), "r");

  Dataset dataset;
  const int N_MAX = 128;
  char c[N_MAX];
  char buffer[N_MAX];
  
  while(fscanf(file, "%s", c) == 1) {
    Datapoint p(4);
    
    int i_p = 0;
    int i_buffer = 0;
    for (int i = 0; i < N_MAX; ++i) {
      if (c[i] == '\0')
	break;
      if (c[i] == ',') {
	buffer[i_buffer] = '\0';
	i_buffer = 0;
	p[i_p++] = atof(buffer);
	assert(i_p <= p.size());
      } else {
	buffer[i_buffer++] = c[i];
      }
    }
    dataset.push_back(p);
  }
  
  /*
  int inputSize = dataset[0].size();
  for (int i = 0; i < inputSize; ++i) {
    real xmin = dataset[0][i];
    real xmax = dataset[0][i];
    for (int j = 0; j < dataset.size(); ++j) {
      xmin = min(xmin, dataset[j][i]);
      xmax = max(xmax, dataset[j][i]);
    }
    assert(xmin != xmax);
    for (int j = 0; j < dataset.size(); ++j)
      dataset[j][i] = 2. * (dataset[j][i] - xmin) / (xmax - xmin) - 1;
  }
  */

  fclose(file);
  return dataset;
}
