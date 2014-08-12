//
// mat.h
// 
// Copyright 2014 Vincent Q. Vu. All rights reserved
//

#include <vector>

template <typename bT, typename Container = std::vector<bT> >
class BlockMat : public BlockBase< BlockMat<bT, Container> > {
public:
  typedef bT block_type;
  typedef typename block_type::elem_type elem_type;
  typedef typename Container::size_type size_type;
  typedef typename Container::const_iterator const_iterator;
  typedef typename Container::iterator iterator;
  typedef typename Container::value_type value_type;
  typedef typename Container::reference reference;
  typedef typename Container::const_reference const_reference;

  // Constructors
  BlockMat() {}
  explicit BlockMat(size_type count) : blocks(count) {}

  // Container accessors
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
  void push_back(const block_type& value) { blocks.push_back(value); }
  void push_back(block_type&& value) { blocks.push_back(value); }
  void resize(size_type count) { blocks.resize(count); }

  // Operators
  inline const BlockMat& operator*=(const double& rhs) {
    for (auto& x : blocks) { x *= rhs; }
    return *this;
  }

  inline const BlockMat& operator/=(const double& rhs) {
    for (auto& x : blocks) { x /= rhs; }
    return *this;
  }

  template <typename functor>
  inline const BlockMat& transform(functor F) {
    for (auto& x : blocks) { x.transform(F); }
    return *this;
  }

  // Expression operators
  template <typename Derived>
  inline const BlockMat& operator=(const BlockBase<Derived>& rhs) {
    blocks.resize(rhs.size());
    auto i = rhs.cbegin();
    for (auto& b : blocks) { 
      b = *i;
      ++i;
    }
    return *this;
  }

  template <typename Derived>
  inline const BlockMat& operator+=(const BlockBase<Derived>& rhs) {
    auto i = rhs.cbegin();
    for (auto& b : blocks) { 
      b += *i;
      ++i;
    }
    return *this;
  }

  template <typename Derived>
  inline const BlockMat& operator-=(const BlockBase<Derived>& rhs) {
    auto i = rhs.cbegin();
    for (auto& b : blocks) { 
      b -= *i;
      ++i;
    }
    return *this;
  }

protected:
  Container blocks;
};

template <typename bT, typename Container >
struct BlockBase_traits< BlockMat<bT, Container> > { 
  typedef bT block_type; 
  typedef typename block_type::elem_type elem_type;
  typedef typename Container::size_type size_type;
  typedef typename Container::const_iterator const_iterator;
};
