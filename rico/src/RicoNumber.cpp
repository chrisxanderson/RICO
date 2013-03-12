#include <iostream>
#include <sstream>
#include <string>

#include "RicoNumber.h"
#include "RicoInteger.h"
#include "RicoFraction.h"
#include "RicoDouble.h"
#include "RicoOperators.h"

#include <math.h>

#include "UnitTest++/UnitTest++.h"

using namespace std;

// RicoNumber ///////////////////////////////////////////////////////////////

RicoPtr RicoNumber::getNeg() const {
  switch(ntype()) {
  case INTEGER:  return ((const RicoInteger &)(*this)).getNeg();
  case FRACTION: return ((const RicoFraction&)(*this)).getNeg();
  case DOUBLE:   return ((const RicoDouble  &)(*this)).getNeg();
  }
  cout << "ERROR: Unknown type in RicoNumber::getNeg" << endl;
  return Rico::NaN;
}

// 5     => {1/5}
// {2/3} => {3/2}
// {1/5} => 5
// 5.0   => 0.2

RicoPtr RicoNumber::getInverse() const {
  switch(ntype()) {
  case INTEGER:  return ((const RicoInteger &)(*this)).getInverse();
  case FRACTION: return ((const RicoFraction&)(*this)).getInverse();
  case DOUBLE:   return ((const RicoDouble  &)(*this)).getInverse();
  }
  cout << "ERROR: Unknown type in RicoNumber::getInverse" << endl;
  return Rico::NaN;  
}

RicoPtr RicoNumber::Log(RicoPtr A) {
  if(!A->isNumber()) {
    cout << "Warning: RicoNumber::Log called with non-number operand." << endl;
    return Rico::Log(A);
  }
  const RicoNumber& a = (const RicoNumber&)(*A);
  switch(a.ntype()) {
  case INTEGER:  return ((const RicoInteger &)a).Log();
  case FRACTION: return ((const RicoFraction&)a).Log();
  case DOUBLE:   return ((const RicoDouble  &)a).Log();
  }
  cout << "ERROR: Unknown type in RicoNumber::Log" << endl;
  return Rico::NaN;
}

RicoPtr RicoNumber::Exp(RicoPtr A) {
  if(!A->isNumber()) {
    cout << "Warning: RicoNumber::Exp called with non-number operand." << endl;
    return Rico::Exp(A);
  }
  const RicoNumber& a = (const RicoNumber&)(*A);
  switch(a.ntype()) {
  case INTEGER:  return ((const RicoInteger &)a).Exp();
  case FRACTION: return ((const RicoFraction&)a).Exp();
  case DOUBLE:   return ((const RicoDouble  &)a).Exp();
  }
  cout << "ERROR: Unknown type in RicoNumber::Exp" << endl;
  return Rico::NaN;
}

// Evaluated in order. i ~ integer, f ~ fraction, d ~ double, x,y ~ non-special number, s ~ sign
// Address special cases {0, 1, +inf, -inf, NaN}
// base e > 0    1 +inf -inf  NaN   y
// vvvv r-----------------------------
//    0 | NaN    0    0  NaN  NaN {NaN, 0) ~ {-,+}
//    1 |   1    1  NaN  NaN  NaN   1
// +inf | NaN +inf  NaN  NaN  NaN {0, inf} ~ {-,+}
// -inf | NaN -inf  NaN  NaN  NaN {0,-inf} ~ {-,+}
//  NaN | NaN  NaN  NaN  NaN  NaN NaN
//    x |   1    x  s*inf  0  NaN x^y

RicoPtr RicoNumber::Power(RicoPtr A, RicoPtr B) {
  if(!A->isNumber() || !B->isNumber()) {
    cout << "WARNING: Called RicoNumber::Power with at least one non-number" << endl;
    return Rico::Power(A,B);
  }

  const RicoNumber& x = (const RicoNumber&) *A;
  const RicoNumber& y = (const RicoNumber&) *B;

  if(x.isZero()) {
    if(y.isZero())   return Rico::NaN;
    if(y.isOne())    return Rico::Zero;
    if(y.isPosInf()) return Rico::Zero;
    if(y.isNegInf()) return Rico::NaN;
    if(y.isNeg())    return Rico::NaN;
    return Rico::Zero;
  }

  if(x.isOne()) {
    if(y.isZero())   return Rico::One;
    if(y.isOne())    return Rico::One;
    if(y.isPosInf()) return Rico::NaN;
    if(y.isNegInf()) return Rico::NaN;
    return Rico::One;
  }

  if(x.isPosInf()) {
    if(y.isZero())   return Rico::NaN;
    if(y.isOne())    return Rico::Fraction(1,0);
    if(y.isPosInf()) return Rico::NaN;
    if(y.isNegInf()) return Rico::NaN;
    if(y.isPos())    return Rico::Fraction(1,0);
    return Rico::Zero;
  }

  if(x.isNegInf()) {
    if(y.isZero())   return Rico::NaN;
    if(y.isOne())    return Rico::Fraction(-1,0);
    if(y.isPosInf()) return Rico::NaN;
    if(y.isNegInf()) return Rico::NaN;
    if(y.isPos())    return Rico::Fraction(-1,0);
    return Rico::Zero;
  }

  if(x.isNaN() || y.isNaN()) return Rico::NaN;

  if(y.isZero())   return Rico::One;
  if(y.isOne())    return A;
  if(y.isPosInf()) {
    if(x.isPos()) return Rico::Fraction( 1,0);
    else          return Rico::Fraction(-1,0);
  }
  if(y.isNegInf()) return Rico::Zero;

  if(x.isDouble()   || y.isDouble())   return RicoDouble  ::Pow(x,y);
  if(x.isFraction() || y.isFraction()) return RicoFraction::Pow(x,y);
  if(x.isInteger()  || y.isInteger())  return RicoInteger ::Pow(x,y);

  cout << "WARNING: RicoNumber::Power fell through cases for " << A << " ^ " << B << endl;
  return Rico::NaN;
}

bool RicoNumber::operator==(RicoPtr That) const {
  if(That->type() != type()) return false;
  const RicoNumber& that = (const RicoNumber&) *That;
  if(isNaN() || that.isNaN()) return false;
  if(isPosInf()) return that.isPosInf();
  if(isNegInf()) return that.isNegInf();
  switch(ntype()) {
  case INTEGER: {
    const RicoInteger& a = (const RicoInteger&) *this;
    switch(that.ntype()) {
    case INTEGER:  {const RicoInteger&  b = (const RicoInteger&)  that; return(a == b);}
    case FRACTION: {const RicoFraction& b = (const RicoFraction&) that; return(b == a);}
    case DOUBLE:   {const RicoDouble&   b = (const RicoDouble&)   that; return(b == a);}
    } break; }
  case FRACTION: {
    const RicoFraction& a = (const RicoFraction&) *this;
    switch(that.ntype()) {
    case INTEGER:  {const RicoInteger&  b = (const RicoInteger&)  that; return(a == b);}
    case FRACTION: {const RicoFraction& b = (const RicoFraction&) that; return(a == b);}
    case DOUBLE:   {const RicoDouble&   b = (const RicoDouble&)   that; return(b == a);}
    } break; }
  case DOUBLE: {
    const RicoDouble& a = (const RicoDouble&) *this;
    switch(that.ntype()) {
    case INTEGER:  {const RicoInteger&  b = (const RicoInteger&)  that; return(a == b);}
    case FRACTION: {const RicoFraction& b = (const RicoFraction&) that; return(a == b);}
    case DOUBLE:   {const RicoDouble&   b = (const RicoDouble&)   that; return(a == b);}
    } break; }
  }
  cout << "RicoNumber::operator==() not yet supported." << endl;
  return false;
}

RicoPtr RicoNumber::operator+(const RicoNumber& B) const {
  switch(ntype()) {
  case INTEGER: {
    const RicoInteger& a = (const RicoInteger&) *this;
    switch(B.ntype()) {
    case INTEGER:  {const RicoInteger&  b = (const RicoInteger&)  B; return(a + b);}
    case FRACTION: {const RicoFraction& b = (const RicoFraction&) B; return(b + a);}
    case DOUBLE:   {const RicoDouble&   b = (const RicoDouble&)   B; return(b + a);}
    } break; }
  case FRACTION: {
    const RicoFraction& a = (const RicoFraction&) *this;
    switch(B.ntype()) {
    case INTEGER:  {const RicoInteger&  b = (const RicoInteger&)  B; return(a + b);}
    case FRACTION: {const RicoFraction& b = (const RicoFraction&) B; return(a + b);}
    case DOUBLE:   {const RicoDouble&   b = (const RicoDouble&)   B; return(b + a);}
    } break; }
  case DOUBLE: {
    const RicoDouble& a = (const RicoDouble&) *this;
    switch(B.ntype()) {
    case INTEGER:  {const RicoInteger&  b = (const RicoInteger&)  B; return(a + b);}
    case FRACTION: {const RicoFraction& b = (const RicoFraction&) B; return(a + b);}
    case DOUBLE:   {const RicoDouble&   b = (const RicoDouble&)   B; return(a + b);}
    } break; }
  }
  cout << "RicoNumber::operator+() not yet supported." << endl;
  return Rico::Zero;
}

RicoPtr RicoNumber::operator*(const RicoNumber& B) const {
  switch(ntype()) {
  case INTEGER: {
    const RicoInteger& a = (const RicoInteger&) *this;
    switch(B.ntype()) {
    case INTEGER:  {const RicoInteger&  b = (const RicoInteger&)  B; return(a * b);}
    case FRACTION: {const RicoFraction& b = (const RicoFraction&) B; return(b * a);}
    case DOUBLE:   {const RicoDouble&   b = (const RicoDouble&)   B; return(b * a);}
    } break; }
  case FRACTION: {
    const RicoFraction& a = (const RicoFraction&) *this;
    switch(B.ntype()) {
    case INTEGER:  {const RicoInteger&  b = (const RicoInteger&)  B; return(a * b);}
    case FRACTION: {const RicoFraction& b = (const RicoFraction&) B; return(a * b);}
    case DOUBLE:   {const RicoDouble&   b = (const RicoDouble&)   B; return(b * a);}
    } break; }
  case DOUBLE: {
    const RicoDouble& a = (const RicoDouble&) *this;
    switch(B.ntype()) {
    case INTEGER:  {const RicoInteger&  b = (const RicoInteger&)  B; return(a * b);}
    case FRACTION: {const RicoFraction& b = (const RicoFraction&) B; return(a * b);}
    case DOUBLE:   {const RicoDouble&   b = (const RicoDouble&)   B; return(a * b);}
    } break; }
  }
  cout << "RicoNumber::operator*() not yet supported." << endl;
  return Rico::Zero;
}

RicoPtr RicoNumber::operator/(const RicoNumber& B) const {
  return this->operator*((const RicoNumber&)(*(B.getInverse()))); // This is ugly.
}

// Testing ////////////////////////////////////////////////////////////////////////////////////////

TEST(TestRicoNumberInverse) {
  // 5     => {1/5}
  CHECK_EQUAL(Rico::Integer(5), Rico::Inverse(Rico::Fraction(1,5)));

  // {2/3} => {3/2}
  CHECK_EQUAL(Rico::Fraction(3,2), Rico::Inverse(Rico::Fraction(2,3)));

  // {1/5} => 5
  CHECK_EQUAL(Rico::Integer(5), Rico::Inverse(Rico::Fraction(1,5)));

  // 5.0   => 0.2
  CHECK_EQUAL(Rico::Double(0.2), Rico::Inverse(Rico::Double(5.0)));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
