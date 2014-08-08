//
// base.h
// 
// Copyright 2014 Vincent Q. Vu. All rights reserved
//

#ifndef __BLOCK_BASE_H
#define __BLOCK_BASE_H

#include <RcppArmadillo.h>
#include <deque>

template <typename bT, typename Container = std::deque<bT> >
class BlockBase {
public:
  typedef typename Container::iterator iterator;
  typedef typename Container::const_iterator const_iterator;
  typedef typename Container::size_type size_type;
  typedef typename Container::value_type value_type;
  typedef typename Container::reference reference;
  typedef typename Container::const_reference const_reference;

  BlockBase() : blocks() {}
  explicit BlockBase(size_type count) : blocks(count) {}


  size_type size() const { return blocks.size(); }

  const_iterator cbegin() const { return blocks.cbegin(); }
  const_iterator cend() const { return blocks.cend(); }
  const_iterator begin() const { return blocks.begin(); }
  const_iterator end() const { return blocks.end(); }
  iterator begin() { return blocks.begin(); }
  iterator end() { return blocks.end(); }

  reference back() { return blocks.back(); }
  const_reference back() const { return blocks.back(); }

  void clear() { return blocks.clear(); }
  void push_back(const bT& value) { blocks.push_back(value); }
  void push_back(bT&& value) { blocks.push_back(value); }
  void resize(size_type count) { blocks.resize(count); }

protected:
  Container blocks;
};

#endif
