#include <opencv/cv.h>
#include <opencv/highgui.h>
#include "common.hpp"
#include <cmath>
#include <iostream>
#include "plot.hpp"
using namespace std;
using namespace cv;
typedef ::real Real;
typedef cv::Mat_<Real> matr;
typedef cv::Mat_<Vec3f> image;

int roundI(Real x) {
  return (int)(floor(x + (Real)0.5));
}

void display(const matr & im, int id) {
  Real minx, maxx;
  minMaxLoc(im, &minx, &maxx);
  matr im2 = (im - minx) / (maxx - minx);
  namedWindow("display" + id);
  imshow("display" + id, im2);
  cvWaitKey(0);
}

void gridPlot(Real(*f)(Real, Real), Real wmin, Real wmax, Real hmin, Real hmax, int id) {
  int w = 400, h = 400; // size of the window in pixel
  //Real wmin = -2, wmax = 2; // range of w
  //Real hmin = -2, hmax = 2; // range of h
  int wres = 50, hres = 50; // number of points in each dim
  Real wrange = wmax - wmin, hrange = hmax - hmin;
  Real wdelta = wrange / (Real)wres, hdelta = hrange / (Real)hres;
  Real wimdelta = (Real)w / (Real)wres, himdelta = (Real)h / (Real)hres;
  matr im(h, w);
  matr imSmall(hres, wres);
  for (int x = 0; x < wres; ++x) {
    for (int y = 0; y < hres; ++y) {
      Real v = f(wmin + ((Real)x+(Real)0.5)*wdelta,
		 hmin + ((Real)y+(Real)0.5)*hdelta);
      im(Range(roundI(y*himdelta), min(h, roundI((y+1)*himdelta))),
	 Range(roundI(x*wimdelta), min(w, roundI((x+1)*wimdelta)))).setTo(v);
      imSmall(y, x) = v;
    }
  }

  //cout << imSmall << endl;
  
  display(im, id);
}

/*
Real testF(Real x, Real y) {
  return x + y;
}

int main() {
  gridPlot(&testF);
  return 0;
}
*/
