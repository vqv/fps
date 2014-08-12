//
// base.h
// 
// Copyright 2014 Vincent Q. Vu. All rights reserved
//

template <typename Derived>
struct BlockBase_traits;

template <typename Derived>
class BlockBase {
public:
  typedef typename BlockBase_traits<Derived>::block_type block_type;
  typedef typename BlockBase_traits<Derived>::elem_type elem_type;
  typedef typename BlockBase_traits<Derived>::size_type size_type;
  typedef typename BlockBase_traits<Derived>::const_iterator const_iterator;

  inline size_type size() const { 
    return static_cast<const Derived&>(*this).size(); 
  }
  inline const_iterator cbegin() const { 
    return static_cast<const Derived&>(*this).cbegin(); 
  }
  inline const_iterator cend() const { 
    return static_cast<const Derived&>(*this).cend();     
  }
  inline const_iterator begin() const { return cbegin(); }
  inline const_iterator end() const { return cend(); }

};
