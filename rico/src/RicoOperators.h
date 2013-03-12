#ifndef __RICO_OPERATORS_H__
#define __RICO_OPERATORS_H__

#include "Rico.h"
#include "RicoNumber.h"
#include "RicoFunction.h"
#include "RicoBasicRV.h"

////////////////////////////////////////////////////////////////////////////////////

static ostream& operator<<(ostream& os, RicoPtr x) {
  os << x->getExpression();
  return os;
}

static ostream& operator<<(ostream& os, const Rico& x) {
  os << x.getExpression();
  return os;
}

static bool operator==(RicoPtr A, RicoPtr B) {
  if(A->isNumber())   return (const RicoNumber&)  (*A) == B;
  if(A->isFunction()) return (const RicoFunction&)(*A) == B;
  if(A->isBasicRV())  return (const RicoBasicRV&) (*A) == B;
  cout << "ERROR: Unknown type " << A->type() << " for operator==()" << endl;
  return false;
}

static bool operator!=(RicoPtr A, RicoPtr B) {return !operator==(A,B);}

static RicoPtr operator+(RicoPtr A, RicoPtr B) {
  if(A->isNumber()) {
    const RicoNumber& a = (const RicoNumber&) *A;
    if(a.isZero()) return B;
    if(B->isNumber()) {
      const RicoNumber& b = (const RicoNumber&) *B;
      if(b.isZero()) return A;
      return a + b;
    }
  }
  return Rico::Add(A,B);
}

static RicoPtr operator*(RicoPtr A, RicoPtr B) {
  if(A->isNumber()) {
    const RicoNumber& a = (const RicoNumber&) *A;
    if(a.isOne()) return B;
    if(B->isNumber()) {
      const RicoNumber& b = (const RicoNumber&) *B;
      return a * b;
    }
    if(a.isZero()) return Rico::Zero;
  }
  return Rico::Multiply(A,B);
}

static RicoPtr operator/(RicoPtr A, RicoPtr B) {
  if(B->isNumber()) {
    const RicoNumber& b = (const RicoNumber&) *B;
    if(b.isOne()) return A;
    if(A->isNumber()) {
      const RicoNumber& a = (const RicoNumber&) *A;
      return a / b;
    }
  }
  return Rico::Divide(A,B);
}

////////////////////////////////////////////////////////////////////////////////////
#endif // __RICO_OPERATORS_H__
