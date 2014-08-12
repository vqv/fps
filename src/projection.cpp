//
// projection.cpp
// 
// Created by Vincent Q. Vu on 2014-07-03
// Copyright 2014 Vincent Q. Vu. All rights reserved
// 
#include "projection.h"
#include "simplex.h"

void FantopeProjection::operator()(arma::mat& x) const {

  arma::uword rank;
  arma::vec eigval;
  arma::mat eigvec;

  arma::eig_sym(eigval, eigvec, x);
  simplex(eigval, rank, d, false);

  // Reconstruct
  x = (
    eigvec.cols(eigvec.n_cols - rank, eigvec.n_cols - 1) *
    arma::diagmat(eigval.subvec(eigval.n_elem - rank, eigval.n_elem - 1)) *
    eigvec.cols(eigvec.n_cols - rank, eigvec.n_cols - 1).t()
  );

  return;
}

void FantopeProjection::operator()(block::mat& x) const {  

  arma::uvec rank;
  block::vec eigval;
  block::mat eigvec;

  block::eig_sym(eigval, eigvec, x);
  simplex(eigval, rank, d, false);

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
        arma::diagmat(di->subvec(di->n_elem - *ri, di->n_elem - 1)) *
        vi->cols(vi->n_cols - *ri, vi->n_cols - 1).t()
      );
    }
    ++ri; ++di; ++vi;
  }

  return;
}

void SingularValueProjection::operator()(arma::mat& x) const {  

  arma::uword rank;
  arma::vec s;
  arma::mat u, v;

  svd(u, s, v, x);
  simplex(s, rank, d, true);

  // Reconstruct
  if (rank < 1) {
    x.zeros();
  } else {
    x = (
      u.cols(0, rank - 1) * 
      arma::diagmat(s.subvec(0, rank - 1)) *
      v.cols(0, rank - 1).t()
    );
  }

  return;
}

void SingularValueProjection::operator()(block::mat& x) const {  

  arma::uvec rank;
  block::vec s;
  block::mat u, v;

  svd(u, s, v, x);
  simplex(s, rank, d, true);

  // Reconstruct
  auto ri = rank.begin();
  auto si = s.cbegin();
  auto ui = u.cbegin();
  auto vi = v.cbegin();
  for (auto& xi : x) {
    if(*ri < 1) {
      xi.zeros();
    } else {
      xi = (
        ui->cols(0, *ri - 1) * 
        arma::diagmat(si->subvec(0, *ri - 1)) *
        vi->cols(0, *ri - 1).t()
      );
    }
    ++ri; ++si; ++ui; ++vi;
  }

  return;
}
