//
// blockmat.cpp
// 
// Created by Vincent Q. Vu on 2014-08-07
// Copyright 2014 Vincent Q. Vu. All rights reserved
// 
#include "blockmat.h"

double BlockMat::sumabs() const {
  double s = 0;
  for( const auto& b : blocks) {
    s += norm(vectorise(b), 1);
  }
  return s;
}

void BlockMat::svd(std::list<arma::mat>& u, std::list<arma::vec>& s, 
                   std::list<arma::mat>& v) const {
  u.clear(); s.clear(); v.clear();
  for (const auto& x : blocks) {
    u.push_back(arma::mat());
    s.push_back(arma::vec());
    v.push_back(arma::mat());
    arma::svd(u.back(), s.back(), v.back(), x);
  }
  return;
}

void BlockMat::eig_sym(std::list<arma::vec>& eigval, 
                       std::list<arma::mat>& eigvec) const {
  eigval.clear(); eigvec.clear();
  for (const auto& x : blocks) {
    eigval.push_back(arma::vec());
    eigvec.push_back(arma::mat());
    arma::eig_sym(eigval.back(), eigvec.back(), x);
  }
}

double dot(const BlockMat& a, const BlockMat& b) {
  double x = 0;
  for (auto ai = a.cbegin(), bi = b.cbegin(); ai != a.cend(); ++ai, ++bi) {
    x += dot(*ai, *bi);
  }
  return x;
}
