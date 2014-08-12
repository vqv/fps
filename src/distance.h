#ifndef __DISTANCE_H
#define __DISTANCE_H

#include "RcppArmadillo.h"
#include "block/block"
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
    for (auto ai = a.cbegin(), bi = b.cbegin(); ai != a.cend(); ++ai, ++bi) {
      x += arma::accu(arma::square(*ai - *bi));
    }
    return std::sqrt(x);
  }
};

#endif
