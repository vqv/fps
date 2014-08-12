//
// binary-ops.h
// 
// Copyright 2014 Vincent Q. Vu. All rights reserved
//

struct OpPlus {
  template <typename T1, typename T2>
  static inline const T1 apply(const T1 a, const T2 b) { return a + b; }
};

struct OpMinus {
  template <typename T1, typename T2>
  static inline const T1 apply(const T1 a, const T2 b) { return a - b; }
};

struct OpTimes {
  template <typename T1, typename T2>
  static inline const T1 apply(const T1 a, const T2 b) { return a * b; }
};

struct OpDiv {
  template <typename T1, typename T2>
  static inline const T1 apply(const T1 a, const T2 b) { return a / b; }
};

template <typename T1, typename T2>
inline BlockBinaryOp<T1, T2, OpPlus>
operator+(const BlockBase<T1>& a, const BlockBase<T2>& b)
{
  return BlockBinaryOp<T1, T2, OpPlus>(a, b);  
}

template <typename T1, typename T2>
inline BlockBinaryOp<T1, T2, OpMinus>
operator-(const BlockBase<T1>& a, const BlockBase<T2>& b)
{
  return BlockBinaryOp<T1, T2, OpMinus>(a, b);  
}

template <typename T1>
inline BlockBinaryOp<T1, typename T1::block_type::elem_type, OpTimes>
operator*(const BlockBase<T1>& a, const typename T1::elem_type b)
{
  return BlockBinaryOp<T1, typename T1::elem_type, OpTimes>(a, b);
}

template <typename T1>
inline BlockBinaryOp<typename T1::block_type::elem_type, T1, OpTimes>
operator*(const typename T1::elem_type a, const BlockBase<T1>& b)
{
  return BlockBinaryOp<T1, typename T1::elem_type, OpTimes>(b, a);
}

template <typename T1>
inline BlockBinaryOp<T1, typename T1::block_type::elem_type, OpDiv>
operator/(const BlockBase<T1>& a, const typename T1::elem_type b)
{
  return BlockBinaryOp<T1, typename T1::elem_type, OpDiv>(a, b);
}
