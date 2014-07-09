//
// utility.cpp
// 
// Created by Vincent Q. Vu on 2014-07-07
// Copyright 2014 Vincent Q. Vu. All rights reserved
// 

#include "utility.h"
#include <RcppArmadillo.h>
#include <cmath>

using namespace arma;

void compute_maxoffdiag(vec& maxoffdiag, const mat& x) {
  maxoffdiag.set_size(x.n_cols);
  for(uword j = 0; j < x.n_cols; j++) {
    maxoffdiag(j) = 0;
    for(uword i = 0; i < x.n_rows; i++) {
      if(i != j && maxoffdiag(j) < std::abs(x(i, j))) {
        maxoffdiag(j) = std::abs(x(i, j));
      }
    }
  }
}

// Generate log linear sequence
void loglinearseq(vec& out, double min, double max, uword n) {
  out.set_size(n);
  for(uword i = 0; i < n; i++) {
    double lx = (double) (n - i - 1) / (double) (n - 1);
    out(i) = (1.0 - lx) * min + lx * max;
  }
}
