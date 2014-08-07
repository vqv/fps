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

void FantopeProjection::operator()(BlockMat& x) const {  

  uvec rank;
  BlockVec eigval;
  BlockMat eigvec;

  eig_sym(eigval, eigvec, x);
  rank = simplex(eigval, d);

  // Reconstruct
  auto ri = rank.begin();
  auto di = eigval.cbegin();
  auto vi = eigvec.cbegin();
  for (auto& xi : x) {
    if(*ri < 1) {
      xi.zeros();
    } else {
      xi = (
        vi->cols(vi->n_cols - *ri, vi->n_cols - 1) * 
        diagmat(di->subvec(di->n_elem - *ri, di->n_elem - 1)) *
        vi->cols(vi->n_cols - *ri, vi->n_cols - 1).t()
      );
    }
    ++ri; ++di; ++vi;
  }

  return;
}

void SingularValueProjection::operator()(mat& x) const {  

  uword rank;
  vec s;
  mat u, v;

  svd(u, s, v, x);
  rank = simplex(s, d, true);

  // Reconstruct
  if (rank < 1) {
    x.zeros();
  } else {
    x = (
      u.cols(0, rank - 1) * 
      diagmat(s.subvec(0, rank - 1)) *
      v.cols(0, rank - 1).t()
    );
  }

  return;
}

void SingularValueProjection::operator()(BlockMat& x) const {  

  uvec rank;
  BlockVec s;
  BlockMat u, v;

  svd(u, s, v, x);
  rank = simplex(s, d, true);

  // Reconstruct
  auto ri = rank.begin();
  auto si = s.cbegin();
  auto ui = u.cbegin();
  auto vi = v.cbegin();
  for (auto& xi : x) {
    if(*ri < 1) {
      xi.zeros();
    } else {
      // ui->shed_cols(*ri, ui->n_cols);
      // si->shed_rows(*ri, si->n_elem);
      // vi->shed_cols(*ri, ui->n_cols);
      // xi = ui * diagmat(si) * vi.t();
      xi = (
        ui->cols(0, *ri - 1) * 
        diagmat(si->subvec(0, *ri - 1)) *
        vi->cols(0, *ri - 1).t()
      );
    }
    ++ri; ++si; ++ui; ++vi;
  }

  return;
}
