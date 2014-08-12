//
// simplex.h
// 
// Copyright 2014 Vincent Q. Vu. All rights reserved
//

#ifndef __SIMPLEX_H
#define __SIMPLEX_H

#include <RcppArmadillo.h>
#include "block/block"
#include <algorithm>

inline double simplex_sum(const arma::vec& x, const double theta) {

  double y = 0.0;

  for (double z : x) {
    z -= theta;
    if (z > 1.0) { y += 1.0; } 
    else if (z > 0.0) { y += z; }
  }

  return y;
}

inline arma::uword simplex_transform(arma::vec& x, arma::uword& rank, const double theta) {
  rank = 0;
  x.transform( [&](double d) {
    d -= theta;
    if (d > 1.0) { d = 1.0; } 
    else if (d < 0.0) { d = 0.0; return d; }
    ++rank;
    return d;
  } );
  return rank;  
}

inline std::set<double> simplex_knots(const arma::vec& x) {
  return std::set<double>(x.begin(), x.end());
}

///////////////////////////////////////////////////////////////////////////////
/// Block versions
///////////////////////////////////////////////////////////////////////////////

inline double simplex_sum(const block::vec& x, const double theta) {

  double y = 0.0;
  for (auto& xi : x) { y += simplex_sum(xi, theta); }

  return y;
}

inline void simplex_transform(block::vec& x, arma::uvec& rank, const double theta) {
  rank.zeros(x.size());
  auto ri = rank.begin();
  for (auto& xi : x) { simplex_transform(xi, *ri, theta); ++ri; }
}

inline std::set<double> simplex_knots(const block::vec& x) {
  std::set<double> knots;
  for (const auto& xi : x) { knots.insert(xi.begin(), xi.end()); }
  return knots;  
}

/**
 * Projects x onto the simplex of z satisfying 0 <= z <= 1, <z,1> = d
 * @param  x            Vector to project
 * @param  d            Target sum
 * @param  interior     Include interior of simplex (default false)
 * @return              Number of nonzeros in the projection
 */
template <typename T1, typename T2>
void simplex(T1& x, T2& rank, double d, bool interior) {

  // Interior of L1 and LInfinity balls
  if(interior && simplex_sum(x, 0) <= d) {
    simplex_transform(x, rank, 0);
  }

  // Let x(j) denote the jth largest x.
  // The root of the piecewise linear function must lie in the interval
  // [x(d) - 1, x(d-1)]
  // The knots of the function consist of the values of x and x - 1
  std::set<double> knots = simplex_knots(x);

  // Construct knots corresponding to x-1 -- sorted in descending order
  std::vector<double> tmp(std::min(knots.size(), (size_t) std::ceil(d)));
  auto i = knots.crbegin();
  for (auto& t : tmp) { t = *i++ - 1; }

  // Eliminate any knots that are smaller than x(d) - 1
  knots.erase(knots.begin(), knots.lower_bound(*tmp.crbegin()));

  // Insert x-1 knots
  for (const auto& t : tmp) { knots.insert(knots.begin(), t); }

  // Find the left-most knot whose function value is < d 
  // This knot is the right endpoint of the interval containing 
  // the solution of the piecewise linear equation
  auto t = std::lower_bound(knots.begin(), knots.end(), d, 
            [&](const double& t, const double& z) {
              return simplex_sum(x, t) >= z;
            });

  // Interpolate
  double a , b, fa, fb, theta;
  b = *t;
  a = *(--t);
  fa = simplex_sum(x, a);
  fb = simplex_sum(x, b);
  theta = a + (b-a) * (d-fa) / (fb-fa);

  simplex_transform(x, rank, theta);
}

#endif
