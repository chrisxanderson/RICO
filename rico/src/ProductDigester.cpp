#include <iostream>
#include <string>
#include <sstream>

#include "RicoOperators.h"
#include "ProductDigester.h"
#include "SumDigester.h"

#include "UnitTest++/UnitTest++.h"

// ProductDigester /////////////////////////////////////////////////////////////////////

//                  a  b c             // header of unique base components
//      a^2*b*c  -> 2  1 1             // We can assume each component us unique (simplest)
//      a*b^{-1} -> 1 -1 o             // some component may overlap with first row.

// Tested Rules:
// 0 x B          => 0
// A x 0          => 0
// 1 x 1          => 1
// 5 x {1/5}      => 1
// 3 x 7          => 21
// 1 x B          => B
// A x 1          => A
// 5 x B          => 5*B
// 5 x 3B         => 15B
// 3A x 5         => 15A
// A x 5          => 5*A
// A / 5          => {1/5}A
// X x Y          => XY
// X*Y x X*Z      => X^2YZ              // Powers captured
// (X*3) x (X*7)  => 21X^2              // Coefficients pulled from any location, put in front.
// X^2 x X*Y      => X^3Y
// X^Y x X^Z      => X^(Y+Z)
// X^2 x X^{-2}   => 1
// X^2Y x X^{-2}Z => YZ
// Y x X / Y      => X
// YX / XY        => 1
// YX / XYZ       => 1/Z
// X/Y * Z/W      => XZ / YW

void ProductDigester::addTerm(RicoPtr A) {
  if(A->isMultiply()) addTerm(((const RicoFunction&)(*A)).get_params());
  else if(A->isNumber()) m_coefficient = m_coefficient * A;
  else if(A->isPower()) {
    const RicoFunction& P = (const RicoFunction&)(*A);
    addTerm(P.x(), P.y());
  } else {
    addTerm(A, Rico::One);
  }
}

void ProductDigester::addInvTerm(RicoPtr A) {
  if(A->isMultiply()) addInvTerm(((const RicoFunction&)(*A)).get_params());
  else if(A->isNumber()) m_coefficient = m_coefficient / A;
  else if(A->isPower()) {
    const RicoFunction& P = (const RicoFunction&)(*A);
    addTerm(P.x(), Negate(P.y()));
  } else 
    addTerm(A, Rico::NegOne);
}

void ProductDigester::addTerm(RicoPtr Base, RicoPtr Exponent) {
  for(int i = 0; i < m_bases.size(); i++) {
    if(m_bases[i] == Base) {
      m_exponents[i].push_back(Exponent);
      return;
    }
  }
  m_bases.push_back(Base);
  m_exponents.push_back(vRicoPtr());
  m_exponents.back().push_back(Exponent);
}

void ProductDigester::addTerm(const vRicoPtr& V) {
  for(vRicoPtr::const_iterator it = V.begin(); it != V.end(); it++) addTerm(*it);
}

void ProductDigester::addInvTerm(const vRicoPtr& V) {
  for(vRicoPtr::const_iterator it = V.begin(); it != V.end(); it++) addInvTerm(*it);
}

RicoPtr ProductDigester::get() const {
  vRicoPtr P, Q;
  P.push_back(m_coefficient);
  for(int i = 0; i < m_bases.size(); i++) {
    RicoPtr base = m_bases[i];
    RicoPtr exp  = SumDigester::Digest(m_exponents[i]);
    if(isNegative(exp)) Q.push_back(Rico::Power(base, Negate(exp)));
    else                P.push_back(Rico::Power(base, exp));
  }
  if(Q.size() == 0) return Rico::Multiply(P);  // can't call RicoSimplify since recursive.
  return Rico::Divide(Rico::Multiply(P), Rico::Multiply(Q));
}

bool ProductDigester::isNegative(RicoPtr A) {
  if(A->isNeg()) return true;
  if(A->isMultiply() && A->x()->isNeg()) return true;
  return false;
}

RicoPtr ProductDigester::Negate(RicoPtr A) {
  if(A->isNumber()) return Rico::Negate(A);
  if(A->isMultiply() && A->x()->isNeg()) {
    const vRicoPtr& P = ((const RicoFunction&)(*A)).get_params();
    vRicoPtr Q;
    for(vRicoPtr::const_iterator it = P.begin(); it != P.end(); it++) {
      if(it == P.begin()) Q.push_back(Rico::Negate(*it));
      else                Q.push_back(*it);
    }
    return Rico::Multiply(Q);
  }
  return Rico::Multiply(Rico::NegOne, A);
}

// Testing //////////////////////////////////////////////////////////////////////////////

TEST(TestProductDigester1) {
  RicoPtr A = Rico::Integer(5);
  RicoPtr B = Rico::Integer(7);

  // 5 == 5
  CHECK_EQUAL(A, ProductDigester(A).get());
  // 5x7 = 35
  CHECK_EQUAL(Rico::Integer(35), ProductDigester(Rico::Multiply(A, B)).get());

  CHECK_EQUAL(Rico::Zero,        ProductDigester(Rico::Multiply(Rico::Zero, B)).get());
  CHECK_EQUAL(Rico::Zero,        ProductDigester(Rico::Multiply(A, Rico::Zero)).get());
  CHECK_EQUAL(A,                 ProductDigester(Rico::Multiply(A, Rico::One)).get());
  CHECK_EQUAL(B,                 ProductDigester(Rico::Multiply(Rico::One, B)).get());
  CHECK_EQUAL(Rico::One,         ProductDigester(Rico::Multiply(Rico::One, Rico::One)).get());
  CHECK_EQUAL(Rico::One,         ProductDigester(Rico::Multiply(Rico::Integer(5), 
								Rico::Fraction(1,5))).get());
  CHECK_EQUAL(Rico::Integer(35), ProductDigester(Rico::Multiply(A, B)).get());
}

TEST(TestProductDigester2) {
  RicoPtr X = Rico::Normal();
  RicoPtr Y = Rico::Cauchy();
  RicoPtr Z = Rico::Triangular();
  RicoPtr T = Rico::Integer(3);
  RicoPtr F = Rico::Integer(5);

  CHECK_EQUAL(Rico::Multiply(F, Y), ProductDigester(Rico::Multiply(F, Y)).get());
  CHECK_EQUAL(Rico::Multiply(F, X), ProductDigester(Rico::Multiply(X, F)).get());

  RicoPtr TX = Rico::Multiply(T, X);

  CHECK_EQUAL(Rico::Multiply(Rico::Integer(15), X), ProductDigester(Rico::Multiply(F, TX)).get());
  CHECK_EQUAL(Rico::Multiply(Rico::Integer(15), X), ProductDigester(Rico::Multiply(TX, F)).get());

  RicoPtr XY = Rico::Multiply(X,Y);
  RicoPtr XZ = Rico::Multiply(X,Z);

  CHECK_EQUAL(XY, ProductDigester(Rico::Multiply(X, Y)).get());
  CHECK_EQUAL(Rico::Add(Rico::Square(X), Y, Z), 
	      ProductDigester(Rico::Multiply(XY, XZ)).get());

  RicoPtr X2Y = Rico::Multiply(Rico::Square(X),Y);
  CHECK_EQUAL(Rico::Add(Rico::Power(X, Rico::Three), Y, Z), 
	      ProductDigester(Rico::Multiply(X2Y, XZ)).get());
}

TEST(TestProductDigester3) {
  RicoPtr X = Rico::Normal();

  RicoPtr X3 = Rico::Multiply(X, Rico::Integer(3));
  RicoPtr X7 = Rico::Multiply(X, Rico::Integer(7));

  CHECK_EQUAL(Rico::Multiply(Rico::Integer(21), Rico::Square(X)), 
	      ProductDigester(Rico::Multiply(X3, X7)).get());
}

TEST(TestProductDigester4) {
  RicoPtr X = Rico::Normal();
  RicoPtr Y = Rico::Cauchy();
  RicoPtr Z = Rico::Triangular();

  RicoPtr XY = Rico::Power(X,Y);
  RicoPtr XZ = Rico::Power(X,Z);

  CHECK_EQUAL(Rico::Power(X, Rico::Add(Y,Z)), ProductDigester(Rico::Multiply(XY, XZ)).get());

  RicoPtr X2 = Rico::Square(X);
  RicoPtr X_2 = Rico::Power(X, Rico::Integer(-2));

  CHECK_EQUAL(Rico::One, ProductDigester(Rico::Multiply(X2, X_2)).get());

  RicoPtr X2Y  = Rico::Multiply(X2, Y);
  RicoPtr X_2Z = Rico::Multiply(X_2, Z);
  CHECK_EQUAL(Rico::Multiply(Y,Z), ProductDigester(Rico::Multiply(X2Y, X_2Z)).get());

  // TEMP
  CHECK_EQUAL(Rico::Fraction(1,5), Rico::One / Rico::Integer(5));

  // X / 5 => {1/5}X
  CHECK_EQUAL(Rico::Multiply(Rico::Fraction(1,5),X), 
	      ProductDigester(X,Rico::Integer(5)).get());

  // Y x X / Y => X
  CHECK_EQUAL(X, ProductDigester(Rico::Multiply(X,Y), Y).get());

  // YX / XY => 1
  CHECK_EQUAL(Rico::One, ProductDigester(Rico::Multiply(X,Y), Rico::Multiply(Y,X)).get());

  // YX / XYZ => 1/Z
  CHECK_EQUAL(Rico::Divide(Rico::One,Z), 
	      ProductDigester(Rico::Multiply(X,Y), Rico::Multiply(Y,X,Z)).get());
}

/////////////////////////////////////////////////////////////////////////////////////////
