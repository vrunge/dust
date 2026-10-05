#ifndef DUST_HW_DUSTIB_H
#define DUST_HW_DUSTIB_H

#include "1D_DUSTib.h"
#include "DUST_1D_HW_Models.h"

namespace dustib {
template<class M> constexpr int vec_id =
  std::is_same<M,VecGaussPolicy>::value ? 0 : std::is_same<M,VecPoissonPolicy>::value ? 1 :
  std::is_same<M,VecExpPolicy>::value ? 2 : std::is_same<M,VecGeomPolicy>::value ? 3 :
  std::is_same<M,VecBernPolicy>::value ? 4 : std::is_same<M,VecBinomPolicy>::value ? 5 :
  std::is_same<M,VecNegbinPolicy>::value ? 6 : 7;

template<int K, class D> struct VectorMath {
  D d;
  using V = hn::Vec<D>;
  V set(double x) const { return hn::Set(d,x); }
  V log(V x) const { return hn::CallLog(d,x); }
  V log1p(V x) const { return hn::CallLog1p(d,x); }
  V exp(V x) const { return hn::CallExp(d,x); }
  V expm1(V x) const { return hn::CallExpm1(d,x); }
  auto valid(V a) const {
    auto f=hn::Lt(hn::Abs(a),set(INFINITY));
    if constexpr (K==0) return f;
    if constexpr (K==4 || K==5) return hn::And(f,hn::And(hn::Ge(a,set(0)),hn::Le(a,set(1))));
    if constexpr (K==3) return hn::And(f,hn::Ge(a,set(1)));
    if constexpr (K==2 || K==7) return hn::And(f,hn::Gt(a,set(0)));
    return hn::And(f,hn::Ge(a,set(0)));
  }
  auto boundary(V a) const {
    if constexpr (K==0 || K==2 || K==7) return hn::MaskFalse(d);
    if constexpr (K==4 || K==5) return hn::Or(hn::Eq(a,set(0)),hn::Eq(a,set(1)));
    return hn::Eq(a,set(K==3 ? 1 : 0));
  }
  V theta(V a) const {
    if constexpr (K==0) return a;
    if constexpr (K==1) return log(a);
    if constexpr (K==2 || K==7) return hn::Div(set(K==2 ? -1 : -.5),a);
    if constexpr (K==4 || K==5) return hn::Sub(log(a),log1p(hn::Neg(a)));
    if constexpr (K==3) return log1p(hn::Div(set(-1),a));
    if constexpr (K==6) return hn::IfThenElse(hn::Lt(a,set(1)),hn::Sub(log(a),log1p(a)),hn::Neg(log1p(hn::Div(set(1),a))));
  }
  V conjugate(V a) const {
    if constexpr (K==0) return hn::Mul(set(.5),hn::Mul(a,a));
    if constexpr (K==1) return hn::Mul(a,hn::Sub(log(a),set(1)));
    if constexpr (K==2 || K==7) return hn::Mul(set(K==2 ? -1 : -.5),hn::Add(log(a),set(1)));
    if constexpr (K==4 || K==5) return hn::Add(hn::Mul(a,log(a)),hn::Mul(hn::Sub(set(1),a),log1p(hn::Neg(a))));
    if constexpr (K==3) return hn::Sub(hn::Mul(hn::Sub(a,set(1)),theta(a)),log(a));
    if constexpr (K==6) return hn::Sub(hn::Mul(a,theta(a)),log1p(a));
  }
  V partition(V r) const {
    if constexpr (K==0) return hn::Mul(set(.5),hn::Mul(r,r));
    if constexpr (K==1) return exp(r);
    if constexpr (K==2) return hn::Neg(log(hn::Neg(r)));
    if constexpr (K==7) return hn::Mul(set(-.5),hn::Add(log(hn::Neg(r)),set(std::log(2.0))));
    if constexpr (K==4 || K==5) return hn::Add(hn::Max(r,set(0)),log1p(exp(hn::Neg(hn::Abs(r)))));
    if constexpr (K==3 || K==6) {
      // Both arguments remain legal even in lanes using the other expression.
      auto l=hn::IfThenElse(hn::Lt(r,set(-std::log(2.0))),log1p(hn::Neg(exp(r))),log(hn::Neg(expm1(r))));
      return hn::Sub(K==3 ? r : set(0),l);
    }
  }
  V mean(V r) const {
    if constexpr (K==0) return r;
    if constexpr (K==1) return exp(r);
    if constexpr (K==2 || K==7) return hn::Div(set(K==2 ? -1 : -.5),r);
    if constexpr (K==4 || K==5) {
      auto z=exp(hn::Neg(hn::Abs(r)));
      return hn::Div(hn::IfThenElse(hn::Ge(r,set(0)),set(1),z),hn::Add(set(1),z));
    }
    if constexpr (K==3) return hn::Div(set(-1),expm1(r));
    if constexpr (K==6) return hn::Div(exp(r),hn::Neg(expm1(r)));
  }
  auto positive(V val,V scale) const {
    return hn::And(hn::Lt(hn::Abs(val),set(INFINITY)),
      hn::Gt(val,hn::Mul(set(guard),hn::Add(set(1),scale))));
  }
};

template<int K,int Version,class D>
auto vector_test(D d, hn::Vec<D> a,hn::Vec<D> b,hn::Vec<D> c,hn::Vec<D> q) {
  VectorMath<K,D> m{d};
  const auto zero=m.set(0),one=m.set(1);
  auto valid=hn::And(hn::And(m.valid(a),m.valid(b)),hn::And(hn::Lt(hn::Abs(c),m.set(INFINITY)),hn::Lt(hn::Abs(q),m.set(INFINITY))));
  auto boundary=m.boundary(a),equal=hn::Eq(a,b);
  auto interior=hn::AndNot(boundary,valid);
  auto safea=hn::IfThenElse(interior,a,m.set(K==3 ? 2 : (K==4 || K==5 ? .5 : 1)));
  auto fa=hn::IfThenElse(boundary,zero,m.conjugate(safea));
  auto pelt=m.positive(hn::Sub(hn::Neg(fa),c),hn::Add(hn::Abs(fa),hn::Abs(c)));
  auto linear=m.positive(hn::Sub(q,c),hn::Add(hn::Abs(c),hn::Abs(q)));
  auto delta=hn::Sub(a,b),e=hn::Sub(c,q);
  auto u=hn::Mul(delta,m.theta(safea));
  auto growing=m.positive(hn::Sub(hn::Neg(u),e),hn::Add(hn::Abs(u),hn::Add(hn::Abs(c),hn::Abs(q))));
  auto regular=hn::AndNot(equal,hn::And(interior,growing));
  auto infinite=hn::MaskFalse(d);
  if constexpr (Math<K>::negative) infinite=hn::And(interior,hn::And(hn::Gt(delta,zero),hn::Le(e,zero)));
  auto denom=hn::IfThenElse(regular,delta,one);
  auto r=hn::Div(hn::Neg(e),denom);
  auto finite=hn::Lt(hn::Abs(r),m.set(INFINITY));
  if constexpr (Math<K>::negative) finite=hn::And(finite,hn::Lt(r,zero));
  auto use=hn::AndNot(infinite,hn::And(regular,finite));
  r=hn::IfThenElse(use,r,m.set(Math<K>::negative ? -1 : 0));
  auto ar=m.partition(r),ra=hn::Mul(r,safea);
  auto j=hn::Sub(hn::Sub(ar,ra),c);
  auto result=m.positive(j,hn::Add(hn::Abs(ar),hn::Add(hn::Abs(ra),hn::Abs(c))));
  if constexpr (Version != 2) {
    auto z=m.mean(r);
    auto good=hn::AndNot(m.boundary(z),m.valid(z));
    z=hn::IfThenElse(good,z,safea);
    auto x=hn::Div(hn::Sub(z,a),denom);
    good=hn::And(good,hn::And(hn::Gt(x,zero),hn::Lt(x,m.set(INFINITY))));
    if constexpr (K != 0) {
      const auto left=hn::Div(hn::Sub(m.set(K==3 ? 1 : 0),a),denom);
      good=hn::And(good,hn::Or(hn::Ge(delta,zero),hn::Lt(x,left)));
      if constexpr (K==4 || K==5) {
        const auto right=hn::Div(hn::Sub(one,a),denom);
        good=hn::And(good,hn::Or(hn::Le(delta,zero),hn::Lt(x,right)));
      }
    }
    x=hn::IfThenElse(good,x,one);
    auto fz=m.conjugate(z);
    if constexpr (Version==1) {
      auto xe=hn::Mul(x,e);
      auto h=hn::Sub(hn::Sub(hn::Neg(fz),c),xe);
      result=hn::Or(hn::And(good,m.positive(h,hn::Add(hn::Abs(fz),hn::Add(hn::Abs(c),hn::Abs(xe))))),hn::AndNot(good,result));
    } else {
      auto w=hn::Div(one,hn::Add(one,x)),mu=hn::Mul(x,w);
      auto wf=hn::Mul(w,fz),mq=hn::Mul(mu,q);
      good=hn::And(good,hn::Lt(mu,one));
      auto g=hn::Add(hn::Sub(hn::Neg(wf),c),mq);
      result=hn::Or(hn::And(good,m.positive(g,hn::Add(hn::Abs(wf),hn::Add(hn::Abs(c),hn::Abs(mq))))),hn::AndNot(good,result));
    }
  }
  return hn::And(valid,hn::Or(pelt,hn::Or(hn::And(equal,linear),hn::Or(infinite,hn::And(use,result)))));
}
} // namespace dustib
#endif
