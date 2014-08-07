//
// projection.cpp
// 
// Created by Vincent Q. Vu on 2014-07-03
// Copyright 2014 Vincent Q. Vu. All rights reserved
// 
#include "projection.h"
#include "simplex.h"

using namespace arma;

void FantopeProjection::operator()(mat& x) const {

  uword rank;
  vec eigval;
  mat eigvec;

  eig_sym(eigval, eigvec, x);
  rank = simplex(eigval, d);

  // Reconstruct
  x = (
    eigvec.cols(eigvec.n_cols - rank, eigvec.n_cols - 1) *
    diagmat(eigval.subvec(eigval.n_elem - rank, eigval.n_elem - 1)) *
    eigvec.cols(eigvec.n_cols - rank, eigvec.n_cols - 1).t()
  );

  return;
}

void SingularValueProjection::operator()(mat& x) const {  

  uword rank;
  vec s;
  mat u, v;

  svd(u, s, v, x);
  rank = simplex(s, d, true);

  // Reconstruct
  x = (
    u.cols(0, rank - 1) * 
    diagmat(s.subvec(0, rank - 1)) *
    v.cols(0, rank - 1).t()
  );

  return;
}
