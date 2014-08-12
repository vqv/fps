//
// binary-expr.h
// 
// Copyright 2014 Vincent Q. Vu. All rights reserved
//

///////////////////////////////////////////////////////////////////////////////
// Block, Block binary operator
///////////////////////////////////////////////////////////////////////////////
template <typename T1, typename T2, typename OpType>
struct BlockBinaryOpIterator {
  BlockBinaryOpIterator(const typename T1::const_iterator& a, 
                        const typename T2::const_iterator& b) 
    : _a(a), _b(b) {}
  inline void operator++() { ++_a; ++_b; }
  inline typename T1::block_type operator*() const { 
    return OpType::apply(*_a, *_b); 
  }
private:
  typename T1::const_iterator _a;
  typename T2::const_iterator _b;
};

template <typename T1, typename T2, typename OpType>
class BlockBinaryOp : public BlockBase< BlockBinaryOp<T1, T2, OpType> > {
public:
  typedef typename T1::block_type block_type;
  typedef typename T1::elem_type elem_type;
  typedef typename T1::size_type size_type;
  typedef BlockBinaryOpIterator<T1, T2, OpType> const_iterator;

  inline BlockBinaryOp(const BlockBase<T1>& a, const BlockBase<T2>& b) 
   : _a(a), _b(b) {}

  inline size_type size() const { 
    return _a.size(); 
  }
  inline const_iterator cbegin() const { 
    return const_iterator(_a.cbegin(), _b.cbegin()); 
  }
  inline const_iterator cend() const { 
    return const_iterator(_a.cend(), _b.cend()); 
  }

private:
  const BlockBase<T1>& _a;
  const BlockBase<T2>& _b;
};

template <typename T1, typename T2, typename OpType>
struct BlockBase_traits< BlockBinaryOp<T1, T2, OpType> > { 
  typedef typename T1::block_type block_type; 
  typedef typename T1::elem_type elem_type;
  typedef typename T1::size_type size_type;
  typedef BlockBinaryOpIterator<T1, T2, OpType> const_iterator;
};


///////////////////////////////////////////////////////////////////////////////
// Block, scalar binary operator
///////////////////////////////////////////////////////////////////////////////
template <typename T1, typename OpType>
struct BlockBinaryOpIterator<T1, typename T1::elem_type, OpType> {
  BlockBinaryOpIterator(const typename T1::const_iterator& a, 
                        const typename T1::elem_type b) 
    : _a(a), _b(b) {}
  inline void operator++() { ++_a; }
  inline typename BlockBase<T1>::block_type operator*() const { 
    return OpType::apply(*_a, _b); 
  }
private:
  typename T1::const_iterator _a;
  typename T1::elem_type _b;
};

template <typename T1, typename OpType>
class BlockBinaryOp<T1, typename T1::elem_type, OpType> 
  : public BlockBase< BlockBinaryOp<T1, typename T1::elem_type, OpType> > {
public:
  typedef typename T1::block_type block_type;
  typedef typename T1::elem_type elem_type;
  typedef typename T1::size_type size_type;
  typedef BlockBinaryOpIterator<T1, elem_type, OpType> const_iterator;

  inline BlockBinaryOp(const BlockBase<T1>& a, const elem_type b) 
   : _a(a), _b(b) {}

  inline size_type size() const { 
    return _a.size();
  }
  inline const_iterator cbegin() const { 
    return const_iterator(_a.cbegin(), _b);
  }
  inline const_iterator cend() const { 
    return const_iterator(_a.cend(), _a);
  }

private:
  const BlockBase<T1>& _a;
  const elem_type _b;
};

template <typename T1, typename OpType>
struct BlockBase_traits< BlockBinaryOp<T1, typename BlockBase<T1>::elem_type, OpType> > { 
  typedef typename T1::block_type block_type; 
  typedef typename T1::elem_type elem_type;
  typedef typename T1::size_type size_type;
  typedef BlockBinaryOpIterator<T1, elem_type, OpType> const_iterator;
};
