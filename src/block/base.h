//
// base.h
// 
// Copyright 2014 Vincent Q. Vu. All rights reserved
//

#ifndef __BLOCK_BASE_H
#define __BLOCK_BASE_H

#include <RcppArmadillo.h>
#include <list>

template <typename bT>
class BlockBase {
public:
  typedef typename std::list<bT>::iterator iterator;
  typedef typename std::list<bT>::const_iterator const_iterator;
  typedef typename std::list<bT>::size_type size_type;

  size_type size() const { return blocks.size(); }

  const_iterator cbegin() const { return blocks.cbegin(); }
  const_iterator cend() const { return blocks.cend(); }
  const_iterator begin() const { return blocks.begin(); }
  const_iterator end() const { return blocks.end(); }
  iterator begin() { return blocks.begin(); }
  iterator end() { return blocks.end(); }

  std::list<bT> blocks;
};

#endif
