//
// simplex.h
// 
// Copyright 2014 Vincent Q. Vu. All rights reserved
//

#ifndef __SIMPLEX_H
#define __SIMPLEX_H

#include <RcppArmadillo.h>
#include "block/block.h"
#include <vector>
#include <algorithm>

inline double simplex_sum(const arma::vec& x, const double theta) {

  double y = 0.0;

  for (double z : x) {
    z -= theta;
    if (z < 0.0) { continue; }
    else if (z < 1.0) { y += z; }
    else { y += 1.0; }
  }

  return y;
}

inline arma::uword simplex_transform(arma::vec& x, arma::uword& rank, const double theta) {
  rank = 0;
  x.transform( [&](double z) {
    z -= theta;
    if (z < 0.0) { z = 0.0; }
    else if (z < 1.0) { ++rank; }
    else { z = 1.0; ++rank; }
    return z;
  } );
  return rank;  
}

inline std::vector<double> simplex_knots(const arma::vec& x) {
  return std::vector<double>(x.begin(), x.end());
}

///////////////////////////////////////////////////////////////////////////////
/// Block versions
///////////////////////////////////////////////////////////////////////////////

inline double simplex_sum(const block::vec& x, const double theta) {

  double y = 0.0;
  for (auto& xi : x) { 
    for (auto z : xi) {
      z -= theta;
      if (z < 0.0) { continue; }
      else if (z < 1.0) { y += z; }
      else { y += 1.0; }
    }
  }

  return y;
}

inline void simplex_transform(block::vec& x, arma::uvec& rank, const double theta) {
  rank.zeros((arma::uword) x.size());
  auto ri = rank.begin();
  for (auto& xi : x) { simplex_transform(xi, *ri, theta); ++ri; }
}

inline std::vector<double> simplex_knots(const block::vec& x) {
  std::vector<double> knots;
  for (const auto& xi : x) { knots.insert(knots.cend(), xi.begin(), xi.end()); }
  return knots;  
}

/**
 * Projects x onto the simplex of z satisfying 0 <= z <= 1, <z,1> = d
 * @param  x            Vector to project
 * @param  rank 
 * @param  d            Target sum (expected: 0 < d < x.size())
 * @param  theta_lower  Lower bound for theta in simplex_transform(x, rank, theta)
 * @return              Number of nonzeros in the projection
 */
template <typename T1, typename T2>
inline 
void simplex(T1& x, T2& rank, double d, 
             double theta_lower = -std::numeric_limits<double>::infinity()) {
  // Let x(j) denote the jth largest x, i.e. x(1) >= x(2) >= ...
  // The root of the piecewise linear function must lie in the interval
  // [x(d) - 1, x(d+1)]
  // The knots of the function consist of the values of x and x - 1
  std::vector<double> knots = simplex_knots(x);

  // Round d up
  std::size_t d0 = std::ceil(d);

  // Partial sort of the d largest
  std::partial_sort(knots.begin(), knots.begin() + d0, knots.end(), 
                    std::greater<double>());

  // Replace x(1), ..., x(d) by x(1) - 1, ..., x(d) - 1
  std::for_each(knots.begin(), knots.begin() + d0, [](double& a) { a -= 1; });

  // Lower bound is x(d) - 1
  theta_lower = std::max(theta_lower, knots[d0 - 1]);

  // Sort knots in increasing order
  std::sort(knots.begin(), knots.end());

  // Eliminate any knots that are:
  //  (1) smaller than x(d) - 1
  //  (2) duplicate
  auto newbegin = std::lower_bound(knots.begin(), knots.end(), theta_lower);
  auto newend = std::unique(newbegin, knots.end());

  // First knot (from smallest to largest) for which simplex_sum(x, t) < d
  // This knot is the right endpoint of the interval containing 
  // the solution of the piecewise linear equation
  auto t = std::lower_bound(newbegin, newend, d, 
            [&](const double& t, const double& z) {
              return simplex_sum(x, t) >= z;
            });

  // Interpolate
  double a = t[-1], b = t[0], fa, fb, theta;

  fa = simplex_sum(x, a);
  fb = simplex_sum(x, b);
  theta = a + (b-a) * (d-fa) / (fb-fa);
  simplex_transform(x, rank, theta);
  return;
}

/**
 * Projects x onto the set of z satisfying 0 <= z <= 1, <z,1> <= d
 * @param  x    Vector to project
 * @param  rank 
 * @param  d    Target sum (expected: 0 < d < x.size())
 * @return      Number of nonzeros in the projection
 */
template <typename T1, typename T2>
inline 
void simplex_interior(T1& x, T2& rank, double d) {

  // Interior of L1 and LInfinity balls
  if (simplex_sum(x, 0) <= d) {
    simplex_transform(x, rank, 0);
    return;
  }
  simplex(x, rank, d, 0.0);
  return;
}

#endif
