//
//  main.cpp
//  fps_profiler
//
//  Copyright (c) 2014 Vincent Vu. All rights reserved.
//

#include <RcppArmadillo.h>
#include <RInside.h>
#include "fps.h"
#include "svps.h"

int main(int argc, char *argv[]) {
    
  RInside R(argc, argv);              // create an embedded R instance
  
  
  for (int i = 0; i < 5; ++i) {
    std::string cmd = "var(matrix(rnorm(1000 * 1000), ncol = 1000))";
    Rcpp::NumericMatrix S = R.parseEval(cmd);
    fps(S, 2, 50, 200, -1, -1, Rcpp::NumericVector::create(), 100, 1e-3, 1);
  }
  
  for (int i = 0; i < 5; ++i) {
    std::string cmd = "matrix(rnorm(1000 * 1000), ncol = 1000)";
    Rcpp::NumericMatrix x = R.parseEval(cmd);
    svps(x, 2, 50, 100, -1, -1, Rcpp::NumericVector::create(), 100, 1e-3, 1);
  }

  return 0;
}
