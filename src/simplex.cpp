//
// simplex.cpp
// 
// Created by Vincent Q. Vu on 2014-03-14
// Copyright 2014 Vincent Q. Vu. All rights reserved
// 
#include "simplex.h"
#include <algorithm>

using namespace arma;

double simplex_sum(const vec& x, const double& theta) {
  double y = 0;

  for(vec::const_iterator i = x.begin(); i != x.end(); i++) {
    double z = *i - theta;
    if(z > 1.0) {
      y += 1.0;
    } else if(z > 0.0) {
      y += z;
    }
  }

  return y;
}

class simplex_compare {
  const vec& x;
public:
  explicit simplex_compare(const vec& x) : x(x) {}
  bool operator() (const double& t, const double& z) const {
    return simplex_sum(x, t) >= z;
  }
};

/**
 * Projects x onto the simplex of z satisfying 0 <= z <= 1, <z,1> = d
 * @param  x            Vector to project
 * @param  d            Target sum
 * @param  interior     Include interior of simplex (default false)
 * @return              Number of nonzeros in the projection
 */
int simplex(vec& x, double d, bool interior) {

  int rank = 0;

  // Interior of L1 and LInfinity balls
  if(interior && simplex_sum(x, 0.0) <= d) {
    for(vec::iterator i = x.begin(); i != x.end(); i++) {
      if(*i > 1.0) {
        *i = 1.0;
      } else if(*i < 0.0) {
        *i = 0.0;
        continue;
      }
      rank++;
    }
    return rank;
  }

  // Construct vector of knots
  vec knots = unique(join_vert(x - 1.0, x));

  // Find the left-most knot whose function value is < d 
  // This knot is the right endpoint of the interval containing 
  // the solution of the piecewise linear equation
  const double* t;
  t = std::lower_bound(knots.begin(), knots.end(), d, simplex_compare(x));

  // Interpolate
  double a , b, fa, fb, theta;
  a = t[-1];
  b = t[0];
  fa = simplex_sum(x, a);
  fb = simplex_sum(x, b);
  theta = a + (b-a) * (d-fa) / (fb-fa);

  // Perform the projection in-place
  for(vec::iterator i = x.begin(); i != x.end(); i++) {
    *i -= theta;
    if(*i > 1.0) {
      *i = 1.0;
    } else if(*i < 0.0) {
      *i = 0.0;
      continue;
    }
    rank++;
  }

  return rank;
}
