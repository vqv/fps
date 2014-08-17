//
// projection.h
// 
// Copyright 2014 Vincent Q. Vu. All rights reserved
//

#ifndef __PROJECTION_H
#define __PROJECTION_H

#include <RcppArmadillo.h>
#include "block/block.h"

struct FantopeProjection {

  FantopeProjection(double d) : d(d) {}
  void operator()(arma::mat& x) const;
  void operator()(block::mat& x) const;

private:
  double d;
};

struct SingularValueProjection {

  SingularValueProjection(double d) : d(d) {}
  void operator()(arma::mat& x) const;
  void operator()(block::mat& x) const;

private:
  double d;
};

#endif
