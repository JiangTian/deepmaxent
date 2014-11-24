#include <iostream>
#include <algorithm>
#include <vector>
#include <limits>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include "common.hpp"
#include "feature.hpp"
using namespace std;

inline real sgn (real x) {
  const real eps = 1e-3;
  return (x >= -eps) ? ((x <= eps) ? 0. : 1.) : -1.;
}

real Step(const int best_k, const int best_j, const vector<vector<real> > & w, 
	  const vector<Feature*> & phi, const real Lambda, const Dataset & S,
	  const Density & pw, const real beta_k) {
  const int k = best_k;
  const int j = best_j;
  const real wkj = w[k][j];
  const real EphiPWkj = pw.EphiPW(k)[j];
  const real EphiSkj = phi[k]->EphiS(S)[j];
  const real pbtp = EphiPWkj + Lambda;
  const real pbtm = EphiPWkj - Lambda;
  const real pbp = EphiSkj + Lambda;
  const real pbm = EphiSkj - Lambda;
  const real e2wL = exp(- 2. * wkj * Lambda);
  //cout << "Lambda=" << Lambda << " pbtp=" << pbtp << " pbm=" << pbm << " pbp="
  //<< pbp << " pbtm=" << pbtm << " e2wL=" << e2wL << endl;
  // TODO!!! this beta is different from the DeepMaxent beta (cf. paper)
  const real beta = (pbtp * pbm * e2wL - pbp * pbtm) / (pbtp * e2wL - pbtm);
  //cout << (pbtp * pbm * e2wL - pbp * pbtm) << " " << (pbtp * e2wL - pbtm) << endl;
  //cout << "beta = " << beta << " " << endl;
  if (abs(beta) <= beta_k) {
    //cout << "case 1" << endl;
    return -wkj;
  } else if (beta > beta_k) {
    //cout << "case 2: " << (pbtm * (beta_k - pbp)) << " " << (pbtp * (beta_k - pbm)) << endl;
    return 0.5 / Lambda * log((pbtm * (beta_k - pbp)) / (pbtp * (beta_k - pbm)));
  } else {
    //cout << "case 3: " << (pbtm * (beta_k + pbp)) << " " <<  (pbtp * (beta_k + pbm)) << endl;
    return 0.5 / Lambda * log((pbtm * (beta_k + pbp)) / (pbtp * (beta_k + pbm)));
  }
}

// S : dataset
// T : number of iterations
// N : array. N[k] is the dimension of the output of tree k
// w : (w_t in paper) : vector of vectors of weights
// d (best_d) : cf paper (d_kj in paper)
// phi : trees
// tolerance : epsilon for comparison with 0
// beta_k : 
// eps : (epsilon_t-1,k,j in paper)
// epsK : array of eps for fixed t and k
// phibar, phibarT : 
// Lambda : total range of the phi (TODO)
// eta : size of the update of w
// inputDim : dimension of the input
Density DeepMaxent(const Dataset & S, int T, int SpSize) {
  // initialization
  int p = 0;
  vector<vector<real> > w;
  vector<Feature*> phi;
  const real tolerance = 1e-3;
  real beta = 0.001; // TODO why ???
  int inputDim = S[0].size();
  // TODO: check that S is consistent (all samples have same dimension)

  // adding all possible features
#if 1
  //  adding raw features
  for (int i = 0; i < inputDim; ++i) {
    Feature* newFeature = new FeatureRaw(i);
    phi.push_back(newFeature);
  }
#endif

#if 1
  //  adding monomial2 features	      
  for (int i = 0; i < inputDim; ++i)
    for (int j = 0; j <= i; ++j) {
      Feature* newFeature = new FeatureMonomial2(i, j);
      phi.push_back(newFeature);
    }
#endif

  //  adding categorical features
#if 0
  for (int i = 0; i < inputDim; ++i) {
    //set category //TODO!!
    S[i]
    Feature* newFeature = new FeatureCategory(i, category);
    phi.push_back(newFeature);
  }
#endif
#if 0
  //  adding threshold features
  {
    vector<real> sortedInput;
    for (int i = 0; i < inputDim; ++i) {
      sortedInput.clear();
      for (int j = 0; j < S.size(); ++j)
	sortedInput.push_back(S[j][i]);
      sort(sortedInput.begin(), sortedInput.end());
      //for (int j = 1; j < sortedInput.size(); ++j) {
      for (int j = 1; j < sortedInput.size(); j += 10) { // NOT ALL THRESHOLDS
	real threshold = 0.5 * (sortedInput[j-1] + sortedInput[j]);
	Feature* newFeature = new FeatureThreshold(i, threshold);
	// TODO (or not): this is never destroyed
	phi.push_back(newFeature);
      }
    }
  }
#endif
  // debug feature
  {
    phi.push_back(new FeatureConstant());
    //phi.push_back(new FeatureThreshold(0, 0));
    //phi.push_back(new FeatureThreshold(1, 0));
  }
  /*
  //  adding hinge features
  {
    real b = 2;
    vector<real> sortedInput;
    for (int i = 0; i < inputDim; ++i) {
      sortedInput.clear();
      for (int j = 0; j < S.size(); ++j)
        sortedInput.push_back(S[j][i]);
      sort(sortedInput.begin(), sortedInput.end());
      for (int j = 1; j < sortedInput.size(); ++j) {
        real threshold = 0.5 * (sortedInput[j-1] + sortedInput[j]);
        Feature* newFeature = new FeatureHinge(i, threshold, b);
        phi.push_back(newFeature);
      }
    }
    }*/

  // compute lambda FOR NOW beta = beta_k \forall k
  real Lambda = 0;
  for (int i = 0; i < S.size(); ++i)
    for (int j = 0; j < S[i].size(); ++j)
      Lambda = max(Lambda, S[i][j]*S[i][j]);
  Lambda += beta + 0.1;

#ifdef USE_DENSITY_GRID
  Dataset worldGrid;
  vector<int> iterator(inputDim, 0);
  int maxK = SpSize+1;
  while (true) {
    Datapoint p(inputDim);
    for (int i = 0; i < inputDim; ++i)
      p[i] = Lambda * (2. * (real)iterator[i] / (maxK-1) - 1.);
    worldGrid.push_back(p);
    int i = 0;
    while((i < inputDim) && (iterator[i] == maxK-1))
      ++i;
    if (i == inputDim)
      break;
    for (int j = 0; j < i; ++j)
      iterator[j] = 0;
    iterator[i] += 1;
  }
#endif

  // initial distribution (all w's = 1)
  for (int i = 0; i < phi.size(); ++i)
    w.push_back(vector<real>(phi[i]->size(), 0.));
#ifdef USE_DENSITY_SAMPLE
  Density pw(w, phi, S, SpSize);
#endif
#ifdef USE_DENSITY_GRID
  Density pw(w, phi, S, worldGrid);
#endif

  cout << "size of phi=" << phi.size() << endl;
  
  // main loop (changing the w's)
  for (int t = 0; t < T; ++t) {
    real best_abs_d = -1.;
    real best_d = 0.;
    real best_beta_k = 0.;
    int best_k, best_j = 0;
    for (int k = 0; k < phi.size(); ++k) {
      real beta_k = 2. * phi[k]->RademacherComplexity() + beta;
      const vector<real> EphiPW = pw.EphiPW(k);
      const vector<real> EphiS = phi[k]->EphiS(S);
      for (int j = 0; j < phi[k]->size(); ++j) {
	real d;
	const real wkj = w[k][j];
	const real eps = EphiPW[j] - EphiS[j];
	if (abs(wkj) > tolerance)
	  d = beta_k * sgn(wkj) + eps;
	else if (abs(eps) <= beta_k)
	  d = 0.;
	else
	  d = - beta_k * sgn(eps) + eps;
	//d = eps + beta_k * sgn(wkj);
	
	//cout << "d=" << d << " eps=" << eps << " beta_k=" << beta_k
	// << " EphiPW=" << EphiPW[j] << " EphiS=" << EphiS[j] << endl;
	if (abs(d) > best_abs_d) {
	//if (t % phi.size() == k) { // DEBUG
	  best_abs_d = abs(d);
	  best_d = d;
	  best_k = k;
	  best_j = j;
	  best_beta_k = beta_k;
	}
      }
    }
    real eta = Step(best_k, best_j, w, phi, Lambda, S, pw, best_beta_k);

    //cout << "k=" << best_k << " j=" << best_j << " eta=" << eta
    //<< " d=" << best_d << endl;

    /*
    { // debug: compute gradient with finite differences
      real eps = 1e-3;
      w[best_k][best_j] += eps;
      Density pwPlus = Density(w, phi, S, SpSize);
      pwPlus.Sp = pw.Sp;
      pwPlus.precomputeFactors(pwPlus.Sp, pwPlus.factorsSp, pwPlus.expFactorsSp,
			       pwPlus.normalizerSp, pwPlus.lognormalizerSp);
      real fPlus = pwPlus.lossS(beta);
      w[best_k][best_j] -= 2.*eps;
      Density pwMinus = Density(w, phi, S, SpSize);
      pwMinus.Sp = pw.Sp;
      pwMinus.precomputeFactors(pwMinus.Sp, pwMinus.factorsSp, pwMinus.expFactorsSp,
			       pwMinus.normalizerSp, pwMinus.lognormalizerSp);
      real fMinus = pwMinus.lossS(beta);
      w[best_k][best_j] += eps;
      cout << "real gradient = " << (fPlus - fMinus) / (2*eps) << endl;
    }
    */
    
    w[best_k][best_j] += eta;

#if 0
    for (int i = 0; i < pw.w.size(); ++i)
      for (int j = 0; j < pw.w[i].size(); ++j)
	cout << pw.w[i][j] << " ";
    cout << endl;
#endif

#ifdef USE_DENSITY_SAMPLE
    pw = Density(w, phi, S, SpSize);
#endif
#ifdef USE_DENSITY_GRID
    pw = Density(w, phi, S, worldGrid);
#endif

    cout << "lossS=" << pw.loss(beta) << endl;

  }

  return pw;
}
