//
// functions.h
// 
// Created by Vincent Q. Vu on 2014-08-07
// Copyright 2014 Vincent Q. Vu. All rights reserved
// 
#include <cmath>

// Trace inner product
template <typename T1, typename T2>
inline 
double dot(const BlockBase<T1>& x, const BlockBase<T2>& y) {
  double z = 0;
  for (auto xi = x.cbegin(), yi = y.cbegin(); xi != x.cend(); ++xi, ++yi) {
    z += dot(*xi, *yi);
  }
  return z;
}

// dot(x'y, x'y) = trace(x'y y' x)
template <typename T1, typename T2>
inline 
double dotsquare(const BlockBase<T1>& x, const BlockBase<T2>& y) {
  double z = 0;
  for (auto xi = x.cbegin(), yi = y.cbegin(); xi != x.cend(); ++xi, ++yi) {
    z += arma::accu(arma::square((*xi).t() * (*yi)));
  }
  return z;
}

// dot(x y', x y') = trace(x'x y' y)
template <typename T1, typename T2>
inline 
double tdotsquare(const BlockBase<T1>& x, const BlockBase<T2>& y) {
  double z = 0;
  for (auto xi = x.cbegin(), yi = y.cbegin(); xi != x.cend(); ++xi, ++yi) {
    z += arma::accu(arma::square((*xi) * (*yi).t()));
  }
  return z;
}

template <typename T>
inline 
double sumabs(const BlockBase<T>& x) {
  double z = 0;
  for (const auto& xi : x) { z += norm(vectorise(xi), 1); }
  return z;
}

template <typename bT, typename vT, typename T>
void svd(BlockMat<bT>& u, BlockMat<vT>& s, BlockMat<bT>& v, 
         const BlockBase<T>& x) {
  u.resize(x.size()); s.resize(x.size()); v.resize(x.size());
  auto ui = u.begin(); auto si = s.begin(); auto vi = v.begin();
  for (const auto& xi : x) {
    arma::svd(*ui, *si, *vi, xi);
    ++ui; ++si; ++vi;
  }
}

template <typename bT, typename vT, typename T>
void eig_sym(BlockMat<vT>& eigval, BlockMat<bT>& eigvec, 
             const BlockBase<T>& x) {
  eigval.resize(x.size()); eigvec.resize(x.size());
  auto di = eigval.begin(); auto vi = eigvec.begin();
  for (const auto& xi : x) {
    arma::eig_sym(*di, *vi, xi);
    ++di; ++vi;
  }
}
