#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>

#include <cmath>
#include <limits>

#include "RicoDouble.h"
#include "RicoInteger.h"
#include "RicoFraction.h"
#include "RicoOperators.h"

#include "UnitTest++/UnitTest++.h"

using namespace std;

const double DOUBLE_COMPARE_TOL = 1E-15; // relative error

// RicoDouble ///////////////////////////////////////////////////////////////

RicoDouble::RicoDouble() : RicoNumber(RicoNumber::DOUBLE), m_x(0) {
  m_double_type = NaN;
}

RicoDouble::RicoDouble(double x) : RicoNumber(RicoNumber::DOUBLE), m_x(x) {
  if(m_x != m_x)                      {m_double_type = NaN;     return;}
  if(abs(m_x)   < DOUBLE_COMPARE_TOL) {m_double_type = ZERO;    return;}
  if(abs(m_x-1) < DOUBLE_COMPARE_TOL) {m_double_type = ONE;     return;}
  if(abs(m_x+1) < DOUBLE_COMPARE_TOL) {m_double_type = NEG_ONE; return;}
  if(numeric_limits<double>::has_infinity) {
    if(m_x > 0 && m_x ==  numeric_limits<double>::infinity()) {m_double_type = POS_INF; return;}
    if(m_x < 0 && m_x == -numeric_limits<double>::infinity()) {m_double_type = NEG_INF; return;}
  }
  m_double_type = REGULAR;
}

RicoPtr RicoDouble::Log() const {
  if(isNaN())    return Rico::NaN;
  if(isNeg())    return Rico::NaN;
  if(isOne())    return Rico::Zero;
  if(isZero())   return Rico::Fraction(-1,0);
  if(isPosInf()) return Rico::Fraction( 1,0);
  return Rico::Double(log(x()));
}

RicoPtr RicoDouble::Exp() const {
  if(isNaN())    return Rico::NaN;
  if(isZero())   return Rico::One;
  if(isNegInf()) return Rico::Zero;
  if(isPosInf()) return Rico::Fraction(1,0);
  return Rico::Double(exp(x()));
}

string RicoDouble::getExpression() const {
  stringstream sstr;
  sstr << x();
  return sstr.str();
}

RicoPtr RicoDouble::Pow(const RicoNumber & a, const RicoNumber & b) {
  if(a.isDouble() && b.isInteger())
    return Rico::Double(pow(((const RicoDouble&)a).x(), ((const RicoInteger&)b).x()));
  if(a.isDouble() && b.isFraction())
    return Rico::Double(pow(((const RicoDouble&)a).x(), ((const RicoFraction&)b).d()));
  if(a.isDouble() && b.isDouble())
    return Rico::Double(pow(((const RicoDouble&)a).x(), ((const RicoDouble&)b).x()));
  if(a.isInteger() && b.isDouble())
    return Rico::Double(pow((double)((const RicoInteger&)a).x(), ((const RicoDouble&)b).x()));
  if(a.isFraction() && b.isDouble())
    return Rico::Double(pow(((const RicoFraction&)a).d(), ((const RicoDouble&)b).x()));
  cout << "ERROR: RicoDouble::Pow() unhandled case." << endl;
  return Rico::NaN;
}

bool RicoDouble::operator==(const RicoDouble& b) const {
  if(isNaN() || b.isNaN())       return false;
  if(isPosInf() && b.isPosInf()) return true;
  if(isNegInf() && b.isNegInf()) return true;
  return abs(x() - b.x())/abs(x()) < DOUBLE_COMPARE_TOL;
}

bool RicoDouble::operator==(const RicoInteger& b) const {
  if(isPosInf() || isNegInf() || isNaN()) return false;
  return abs(x() - (double)b.x())/abs(x()) < DOUBLE_COMPARE_TOL;
}

bool RicoDouble::operator==(const RicoFraction& b) const {
  if(isNaN() || b.isNaN())       return false;
  if(isPosInf() && b.isPosInf()) return true;
  if(isNegInf() && b.isNegInf()) return true;
  if(isZero()   && b.isZero())   return true;
  return abs(x() - (double)b.x() / (double)b.y())/abs(x()) < DOUBLE_COMPARE_TOL;
}

RicoPtr RicoDouble::operator+(const RicoInteger & b) const {
  if(isPosInf() || isNegInf() || isNaN()) return Rico::Double(x());
  if(isZero()) return Rico::Integer(b.x());
  if(isOne())  return Rico::Integer(1 + b.x());
  return Rico::Double(x() + b.x());
}

RicoPtr RicoDouble::operator+(const RicoFraction& b) const {
  if(isNaN() || b.isNaN())       return Rico::Double();
  if(isPosInf() && b.isPosInf()) return Rico::Double(x());
  if(isPosInf() && b.isNegInf()) return Rico::Double();
  if(isNegInf() && b.isPosInf()) return Rico::Double();
  if(isNegInf() && b.isNegInf()) return Rico::Double(x());
  if(b.isPosInf())               return Rico::Fraction( 1,0);
  if(b.isNegInf())               return Rico::Fraction(-1,0);
  if(isZero())                   return Rico::Fraction(b.x(), b.y());
  if(isOne())                    return Rico::Fraction(b.x() + b.y(), b.y());
  return Rico::Double(x() + b.x()/b.y());
}

RicoPtr RicoDouble::operator+(const RicoDouble  & b) const {
  return Rico::Double(x() + b.x());
}

RicoPtr RicoDouble::operator*(const RicoInteger & b) const {
  if(isPosInf() || isNegInf() || isNaN()) return Rico::Double(x());
  if(isZero()) return Rico::Zero;
  if(isOne())  return Rico::Integer(b.x());
  return Rico::Double(x() * b.x());
}

RicoPtr RicoDouble::operator*(const RicoFraction& b) const {
  if(isNaN() || b.isNaN())       return Rico::Double();
  if(isPosInf() && b.isPosInf()) return Rico::Double(x());
  if(isPosInf() && b.isNegInf()) return Rico::Double(-x());
  if(isNegInf() && b.isPosInf()) return Rico::Double(-x());
  if(isNegInf() && b.isNegInf()) return Rico::Double(-x());
  if(b.isPosInf())               return Rico::Fraction( (x() < 0 ? -1:  1)*b.x(), 0);
  if(b.isNegInf())               return Rico::Fraction( (x() < 0 ?  1: -1)*b.x(),0);
  if(isZero())                   return Rico::Zero;
  if(isOne())                    return Rico::Fraction(b.x(), b.y());
  return Rico::Double(x() * b.x()/b.y());
}

RicoPtr RicoDouble::operator*(const RicoDouble  & b) const {
  return Rico::Double(x() * b.x());
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TEST(TestDoubleLog) {
  // Log(NaN) => NaN
  CHECK(RicoDouble(std::numeric_limits<double>::quiet_NaN()).Log()->isNaN());
  
  // Log(-1) => NaN
  CHECK(RicoDouble(-1).Log()->isNaN());

  // Log(0) => -Infty
  CHECK(RicoDouble(0).Log()->isNegInf());

  // Log(1) => 0
  CHECK(RicoDouble(1).Log()->isZero());

  // Log(+Infty) => +Infty
  CHECK(RicoDouble(numeric_limits<long double>::infinity()).Log()->isPosInf());
}

TEST(TestDoubleExp) {
  // Exp(NaN) => NaN
  CHECK(RicoDouble(std::numeric_limits<double>::quiet_NaN()).Exp()->isNaN());

  // Exp(0) => 1
  CHECK(RicoDouble(0).Exp()->isOne());

  // Exp(-Infty) => 0
  CHECK(RicoDouble(-numeric_limits<long double>::infinity()).Exp()->isZero());

  // Exp(+Infty) => +Infty
  CHECK(RicoDouble(numeric_limits<long double>::infinity()).Exp()->isPosInf());
}

////////////////////////////////////////////////////////////////////////////////////////////////////
