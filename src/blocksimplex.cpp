//
// blocksimplex.cpp
// 
// Created by Vincent Q. Vu on 2014-08-07
// Copyright 2014 Vincent Q. Vu. All rights reserved
// 
#include "simplex.h"
#include <algorithm>

using namespace arma;
using namespace std;

double simplex_sum(const list<vec>& x, const double& theta) {

  double y = 0.0;
  for (const auto& xi : x) {
    for (double z : xi) {
      z -= theta;
      if (z > 1.0) {
        y += 1.0;
      } else if (z > 0.0) {
        y += z;
      }
    }
  }

  return y;
}

uvec simplex(list<vec>& x, double d, bool interior) {

  uvec rank(x.size(), fill::zeros);
  auto ri = rank.begin();

  // Interior of L1 and LInfinity balls
  if (interior && simplex_sum(x, 0.0) <= d) {
    for (auto& xi : x) {
      xi.transform( [&](double d) {
        if (d > 1.0) {
          d = 1.0;
        } else if (d < 0.0) {
          d = 0.0;
          return d;
        }
        ++(*ri);
        return d;
      } );
      ++ri;
    }
    return rank;
  }

  // Construct set of knots, sorted in ascending order
  std::set<double> knots;
  for (const auto& xi : x) { 
    knots.insert(xi.begin(), xi.end());
    for (double z : xi) { knots.insert(z - 1.0); }
  }

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

  // Perform the projection in-place
  for (auto& xi : x) {
    xi.transform( [&](double d) {
      d -= theta;
      if (d > 1.0) {
        d = 1.0;
      } else if (d < 0.0) {
        d = 0.0;
        return d;
      }
      ++(*ri);
      return d;
    } );
    ++ri;
  }
  return rank;

}
