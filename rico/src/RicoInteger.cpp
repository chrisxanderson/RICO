#include "RicoInteger.h"
#include "RicoDouble.h"
#include "RicoOperators.h"

#include <cmath>

#include <iostream>
#include <sstream>
#include <string>

#include "UnitTest++/UnitTest++.h"

using namespace std;

// RicoInteger ///////////////////////////////////////////////////////////////

string RicoInteger::getExpression() const {
  stringstream sstr;
  sstr << x();
  return sstr.str();
}

RicoPtr RicoInteger::Log() const {
  if(isNeg())  return Rico::NaN;
  if(isOne())  return Rico::Zero;
  if(isZero()) return Rico::Fraction(-1,0);
  return Rico::Double(log((double)x()));
}

RicoPtr RicoInteger::Exp() const {
  if(isZero()) return Rico::One;
  return Rico::Double(exp((double)x()));
}

RicoPtr RicoInteger::Pow(const RicoNumber& a, const RicoNumber& b) {
  if(!a.isInteger() && !b.isInteger()) {
    cout << "Error: RicoInteger::Pow() handles only RicoInteger operands." << endl;
    return Rico::NaN;
  }
  int x = ((const RicoInteger&)a).x();
  int y = ((const RicoInteger&)b).x();

  bool  neg = (y < 0); if(neg) y = -y;
  int    zi = 1; for(int i = 0; i < y; i++) zi *= x;
  double zd = pow((double)x, y);

  if(RicoDouble(zd) == RicoDouble(zi)) {
    if(neg) return Rico::Fraction(1,zi);
    else    return Rico::Integer(zi);
  } else {
    if(neg) return Rico::Double(1/zd);
    else    return Rico::Double(zd);
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

TEST(TestIntegerLog) {
  // Log(-1) => NaN
  CHECK(RicoInteger(-1).Log()->isNaN());

  // Log(0) => -Infty
  CHECK(RicoInteger( 0).Log()->isNegInf());

  // Log(1) => 0
  CHECK(RicoInteger( 1).Log()->isZero());
}

TEST(TestIntegerExp) {
  // Exp(0) => 1
  CHECK(RicoInteger(0).Exp()->isOne());
}

TEST(TestIntegerPow) {
  // 2^3 = 8
  CHECK_EQUAL(Rico::Integer(8), RicoInteger::Pow(RicoInteger(2), RicoInteger(3)));

  // 2^{-3} = 1/8
  CHECK_EQUAL(Rico::Fraction(1,8), RicoInteger::Pow(RicoInteger(2), RicoInteger(-3)));

  // 2^100 = 1.267650600228229401496703205376 × 10^30
  CHECK_EQUAL(Rico::Double(1.267650600228229401496703205376E30), 
	      RicoInteger::Pow(RicoInteger(2), RicoInteger(100)));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
