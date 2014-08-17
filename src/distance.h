#ifndef __DISTANCE_H
#define __DISTANCE_H

#include "RcppArmadillo.h"
#include "block/block.h"
#include <cmath>

struct FrobeniusDistance {
  template <typename M> 
  inline double operator()(const M& a, const M& b) const {
    return std::sqrt(arma::accu(arma::square(a - b)));
  }

  template <typename bT> 
  inline double operator()(const block::BlockMat<bT>& a, 
                           const block::BlockMat<bT>& b) const {
    double x = 0;
    auto bi = b.cbegin();
    for (const auto& ai : a) { 
      x += arma::accu(arma::square(ai - *bi));
      ++bi;
    }
    return std::sqrt(x);
  }
};

#endif
