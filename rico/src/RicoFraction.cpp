#include <iostream>
#include <sstream>
#include <string>

#include <cmath>

#include "RicoFraction.h"
#include "RicoInteger.h"
#include "RicoDouble.h"
#include "RicoOperators.h"

#include "UnitTest++/UnitTest++.h"

using namespace std;

// RicoFraction ///////////////////////////////////////////////////////////////

RicoFraction::RicoFraction(int x, int y) : RicoNumber(RicoNumber::FRACTION), m_x(x), m_y(y) {
  if(y == 0) return;
  int d = gcd(x,y);
  if(d != 1) {
    m_x = x/d;
    m_y = y/d;
  }
  if(m_y < 0) {m_x = -m_x; m_y = -m_y;}  // ensure sign carried by numerator.
}

RicoPtr RicoFraction::Log() const {
  if(isNaN())    return Rico::NaN;
  if(isNeg())    return Rico::NaN;
  if(isOne())    return Rico::Zero;
  if(isZero())   return Rico::Fraction(-1,0);
  if(isPosInf()) return Rico::Fraction(1,0);
  return Rico::Double(log(((double)x())/((double)y())));
}

RicoPtr RicoFraction::Exp() const {
  if(isNaN())    return Rico::NaN;
  if(isZero())   return Rico::One;
  if(isNegInf()) return Rico::Zero;
  if(isPosInf()) return Rico::Fraction(1,0);
  return Rico::Double(exp(((double)x())/((double)y())));
}

string RicoFraction::getExpression() const {
  stringstream sstr;
  sstr << "(" << x() << "/" << y() << ")";
  return sstr.str();
}

// Assume fractions are proper (not integers) i.e. upconverted.
RicoPtr RicoFraction::Pow(const RicoNumber& a, const RicoNumber& b) {
  if(a.isFraction() && b.isInteger()) {
    int xx = ((const RicoFraction&)a).x();
    int xy = ((const RicoFraction&)a).y();

    RicoPtr X = RicoInteger::Pow(RicoInteger(xx), b);
    if(X->isDouble()) return RicoDouble::Pow(RicoDouble(xx/xy), b);
    RicoPtr Y = RicoInteger::Pow(RicoInteger(xy), b);
    if(Y->isDouble()) return RicoDouble::Pow(RicoDouble(xx/xy), b);
    if(!X->isInteger() || !Y->isInteger()) {
      cout << "ERROR: RicoFraction::Pow should have reduced to int^int, but didn't" << endl;
      return Rico::NaN;
    }
    int x = ((const RicoInteger&)(*X)).x();
    int y = ((const RicoInteger&)(*Y)).x();
    return Rico::Fraction(x,y);
  }
  if(a.isFraction() && b.isFraction()) {
    if(a.isNeg()) return Rico::NaN;
    return Rico::Double(pow(((const RicoFraction&)a).d(), ((const RicoFraction&)b).d()));
  }
  if(a.isInteger() && b.isFraction()) {
    if(a.isNeg()) return Rico::NaN;
    return Rico::Double(pow((double)((const RicoInteger&)a).x(), ((const RicoFraction&)b).d()));
  }
  cout << "ERROR: RicoFraction::Pow() unhandled case." << endl;
  return Rico::NaN;
}

bool RicoFraction::operator==(const RicoInteger & b) const {return x() == b.x() && y() == 1;}
bool RicoFraction::operator==(const RicoFraction& b) const {return x() == b.x() && y() == b.y();}

RicoPtr RicoFraction::upConvert(RicoPtr F) {
  if(F->isFraction()) {
    const RicoFraction& f = (const RicoFraction&)(*F);
    if(f.y() == 1) return Rico::Integer(f.x());
  }
  return F;
}

RicoPtr RicoFraction::operator+(const RicoInteger & b) const {
  RicoPtr Z = Rico::Fraction(x() + b.x()*y(), y());
  return RicoFraction::upConvert(Z);
}
RicoPtr RicoFraction::operator+(const RicoFraction& b) const {
  RicoPtr Z = Rico::Fraction(x()*b.y() + b.x()*y(), y()*b.y());
  return RicoFraction::upConvert(Z);
}
RicoPtr RicoFraction::operator*(const RicoInteger& b) const {
  RicoPtr Z = Rico::Fraction( x() * b.x(), y());
  return RicoFraction::upConvert(Z);
}
RicoPtr RicoFraction::operator*(const RicoFraction& b) const {
  RicoPtr Z = Rico::Fraction(x()*b.x(), y()*b.y());
  return RicoFraction::upConvert(Z);
}

// Util ////////////////////////////////////////////////////////////////////////

int RicoFraction::gcd(int a, int b) {
  int t;
  while(b != 0) {
    t = b;
    b = a % b;
    a = t;
  }
  return a;
}

///////////////////////////////////////////////////////////////////////////////

TEST(TestRicoFraction) {
  RicoFraction * f;

  f = new RicoFraction(15, 36);
  CHECK(f->x() == 5 && f->y() == 12);

  f = new RicoFraction(3, 4);
  CHECK(f->x() == 3 && f->y() == 4);

  f = new RicoFraction(1071, 462);
  CHECK(f->x() == 51 && f->y() == 22);

  f = new RicoFraction(3, -2);
  CHECK(f->x() == -3 && f->y() == 2);

  f = new RicoFraction(-3, 2);
  CHECK(f->x() == -3 && f->y() == 2);

  f = new RicoFraction(-3, -2);
  CHECK(f->x() == 3 && f->y() == 2);
}

TEST(TestRicoFractionAddition) {
  RicoPtr F1 = Rico::Fraction(3,4);
  RicoPtr F2 = Rico::Fraction(0,4);
  RicoPtr F3 = Rico::Fraction(3,0);
  RicoPtr F4 = Rico::Fraction(0,0);
  RicoPtr F5 = Rico::Fraction(1,4);
  RicoPtr I1 = Rico::Integer(7);

  CHECK_EQUAL(Rico::Fraction(31, 4), F1 + I1);
  CHECK_EQUAL(Rico::One,             F1 + F5);
  CHECK_EQUAL(I1,                    F2 + I1);
  CHECK_EQUAL(F1,                    F2 + F1);
  CHECK_EQUAL(Rico::Zero,            F2);
  CHECK_EQUAL(Rico::Zero,            F2 + F2);
  CHECK(F3->isPosInf());
  CHECK(F4->isNaN());
  CHECK((F1+F4)->isNaN());
  CHECK((I1+F4)->isNaN());
  CHECK((F4+F4)->isNaN());
}

TEST(TestFractionLog) {
  // Log(NaN) => NaN
  CHECK(RicoFraction(0,0).Log()->isNaN());
  
  // Log(-1) => NaN
  CHECK(RicoFraction(-1,1).Log()->isNaN());

  // Log(0) => -Infty
  CHECK(RicoFraction( 0,1).Log()->isNegInf());

  // Log(1) => 0
  CHECK(RicoFraction( 1,1).Log()->isZero());

  // Log(+Infty) => +Infty
  CHECK(RicoFraction( 1,0).Log()->isPosInf());
}

TEST(TestFractionExp) {
  // Exp(NaN) => NaN
  CHECK(RicoFraction(0,0).Exp()->isNaN());

  // Exp(0) => 1
  CHECK(RicoFraction(0,1).Exp()->isOne());

  // Exp(-Infty) => 0
  CHECK(RicoFraction(-1,0).Exp()->isZero());

  // Exp(+Infty) => +Infty
  CHECK(RicoFraction(1,0).Exp()->isPosInf());
}

TEST(TestFractionPow) {
  // 1/2 ^ 2 = 1/4
  CHECK_EQUAL(Rico::Fraction(1,4), RicoFraction::Pow(RicoFraction(1,2), RicoInteger(2)));

  // 1/2 ^ 1/2 = 0.7071067811865475244008443621048490392848359376884740
  CHECK_EQUAL(Rico::Double(0.7071067811865475244008443621048490392848359376884740),
	      RicoFraction::Pow(RicoFraction(1,2), RicoFraction(1,2)));

  // 2 ^ 1/2 = 1.4142135623730950488016887242096980785696718753769480
  CHECK_EQUAL(Rico::Double(1.4142135623730950488016887242096980785696718753769480),
	      RicoFraction::Pow(RicoInteger(2), RicoFraction(1,2)));
}

///////////////////////////////////////////////////////////////////////////////
