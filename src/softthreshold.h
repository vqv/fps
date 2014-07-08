#ifndef __SOFTTHRESHOLD_H
#define __SOFTTHRESHOLD_H

#include <algorithm>

// Proximal operator for \lambda |x|_1
class SoftThresholdOp
{
public:
  SoftThresholdOp(const double& z) : z(z) {}
  const double operator() (const double& x) const {
    return ((x > 0) - (x < 0)) * std::max(0.0, std::abs(x) - z);
  }
private:
    const double z;
};

// Proximal operator for \lambda (\alpha |x|_1 + 0.5 (1-\alpha) |x|_2^2)
class ElasticSoftThresholdOp
{
public:
  ElasticSoftThresholdOp(const double& z, const double& alpha) :
    z1(z * alpha), z2(1.0 / (1.0 + 0.5 * (1.0 - alpha) * z)) {}
  const double operator() (const double& x) const {
    return ((x > 0) - (x < 0)) * z2 * std::max(0.0, std::abs(x) - z1);
  }
private:
  const double z1, z2;
};

#endif
