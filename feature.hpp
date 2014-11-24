#ifndef __FEATURE_HPP__
#define __FEATURE_HPP__

#include "common.hpp"
#include "density.hpp"

class Feature {
public:
  Feature() {};
  virtual std::vector<real> eval(const Datapoint & X) const = 0;
  virtual int size() const = 0;
  virtual real RademacherComplexity() const = 0;
  std::vector<real> EphiS(const Dataset & S) const;
};

class FeatureConstant : public Feature {
private:
public:
  inline FeatureConstant()
    :Feature() {}
  virtual std::vector<real> eval(const Datapoint & X) const;
  virtual int size() const;
  virtual real RademacherComplexity() const;
};

class FeatureRaw : public Feature {
private:
  const int i;
public:
  inline FeatureRaw(int i)
    :Feature(), i(i) {}
  virtual std::vector<real> eval(const Datapoint & X) const;
  virtual int size() const;
  virtual real RademacherComplexity() const;
};

#if 0
class FeatureCategory : public Feature {
private:
  const int i;
  const string category;
public:
  inline FeatureCategory(int i, string category)
    :Feature(), i(i), category(category) {}
  virtual std::vector<real> eval(const Datapoint & X) const;
  virtual int size() const;
  virtual real RademacherComplexity() const;
}
#endif

class FeatureMonomial2 : public Feature {
private:
  const int i, j;
public:
  inline FeatureMonomial2(int i, int j)
    :Feature(), i(i), j(j) {}
  virtual std::vector<real> eval(const Datapoint & X) const;
  virtual int size() const;
  virtual real RademacherComplexity() const;
};

class FeatureThreshold : public Feature {
private:
  const int i;
  const real threshold;
public:
  inline FeatureThreshold(int i, real threshold)
    :Feature(), i(i), threshold(threshold) {}
  virtual std::vector<real> eval(const Datapoint & X) const;
  virtual int size() const;
  virtual real RademacherComplexity() const;
};
  
class FeatureHinge : public Feature {
private:
  const int i;
  const real threshold, b;
public:
  inline FeatureHinge(int i, real threshold, real b)
    :Feature(), i(i), threshold(threshold), b(b) {}
  virtual std::vector<real> eval(const Datapoint & X) const;
  virtual int size() const;
  virtual real RademacherComplexity() const;
};
#endif
