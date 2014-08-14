//
// softthreshold.h
// 
// Copyright 2014 Vincent Q. Vu. All rights reserved
//

#ifndef __SOFTTHRESHOLD_H
#define __SOFTTHRESHOLD_H

#include <cmath>

// Proximal operator for \lambda |x|_1
struct SoftThresholdOp
{
  SoftThresholdOp(const double z) : z(z) {}
  inline double operator() (const double x) const {
    return std::copysign(std::fdim(std::fabs(x), z), x);
  }

private:
    const double z;
};

// Proximal operator for \lambda (\alpha |x|_1 + 0.5 (1-\alpha) |x|_2^2)
struct ElasticSoftThresholdOp
{
  ElasticSoftThresholdOp(const double z, const double alpha) :
    z1(z * alpha), z2(1.0 / (1.0 + 0.5 * (1.0 - alpha) * z)) {}
  inline double operator() (const double x) const {
    return std::copysign(z2 * std::fdim(std::fabs(x), z1), x);
  }

private:
  const double z1, z2;
};

struct EntrywiseSoftThreshold
{
  EntrywiseSoftThreshold(const double lambda) : lambda(lambda) {}
  template <typename T>
  inline void operator()(T& x, const double z) const {
    x.transform( SoftThresholdOp(z * lambda) );
  }
private:
  const double lambda;
};

struct ColumnSoftThreshold
{
  ColumnSoftThreshold(const double lambda) : lambda(lambda) {}
  inline void operator()(arma::mat& x, const double z) const {
    double c = z * lambda;
    for (arma::uword j = 0; j < x.n_cols; j++) {
      double y = norm(x.col(j));
      if (y <= c) { continue; }
      else { x.col(j) = x.col(j) * (c / y); }
    }
  }

private:
  const double lambda;
};

#endif
