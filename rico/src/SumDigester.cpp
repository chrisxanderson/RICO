#include <iostream>
#include <string>
#include <sstream>

#include "SumDigester.h"
#include "RicoSimplify.h"
#include "RicoOperators.h"

#include "UnitTest++/UnitTest++.h"

// SumDigester /////////////////////////////////////////////////////////////////////////
// (2a + b  + c) + (a + {-1}b)
//
//        1 | a b c             // extendable header of unique base components with number slot
//        7 | 3 0 1             // cumulative coefficients of base elements as terms added.

// Tested Rules:
// 0 + B          => B
// A + 0          => A
// 0 + 0          => 0
// 5 + {-5}       => 0
// 3 + 7          => 10
// 5 + B          => 5 + B
// A + 5          => 5 + A
// 3A + 5         => 5 + 3A
// A + B          => A + B
// (A+B) + (A+C)  => 2A + B + C
// (A+3) + (A+7)  => 10 + 2A
// 2A + A + B     => 3A + B
// AB + AC        => AB + AC
// 2A + {-2}A     => 0

void SumDigester::addTerm(RicoPtr A) {
  if(A->isNumber()) {
    m_constant = m_constant + A;
    return;
  }

  if(A->isAdd()) {
    const vRicoPtr& P = ((const RicoFunction&)(*A)).get_params();
    addTerm(P);
    return;
  }

  if(A->isMultiply()) {
    RicoPtr C = Rico::One;
    const vRicoPtr& P = ((const RicoFunction&)(*A)).get_params();
    vRicoPtr        Q;
    for(vRicoPtr::const_iterator it = P.begin(); it != P.end(); it++) {
      if((*it)->isNumber()) C = C * (*it);
      else Q.push_back(*it);
    }
    addTerm(C, Rico::Multiply(Q));
    return;
  }

  addTerm(Rico::One, A);
}

void SumDigester::addTerm(const vRicoPtr& V) {
  for(vRicoPtr::const_iterator it = V.begin(); it != V.end(); it++) 
    addTerm(*it);
}

void SumDigester::addTerm(RicoPtr Coefficient, RicoPtr Base) {
  for(int i = 0; i < m_bases.size(); i++) {
    if(m_bases[i] == Base) {
      m_coefficients[i] = m_coefficients[i] + Coefficient;
      return;
    }
  }
  m_coefficients.push_back(Coefficient);
  m_bases.push_back(Base);
}

RicoPtr SumDigester::get() const {
  vRicoPtr P;
  P.push_back(m_constant);
  for(int i = 0; i < m_bases.size(); i++) 
    P.push_back(Rico::Multiply(m_coefficients[i], m_bases[i]));
  return Rico::Add(P);
}

// Testing //////////////////////////////////////////////////////////////////////////////

TEST(TestSumDigester1) {
  RicoPtr A   = Rico::Integer(5);
  RicoPtr B   = Rico::Integer(7);
  RicoPtr X   = Rico::Normal();
  RicoPtr Y   = Rico::Cauchy();
  RicoPtr Z   = Rico::Triangular();

  RicoPtr X3  = Rico::Multiply(Rico::Integer( 3), X);
  RicoPtr X4  = Rico::Multiply(Rico::Integer( 4), X);
  RicoPtr X3n = Rico::Multiply(Rico::Integer(-3), X);
  RicoPtr XY  = Rico::Add(X, Y);
  RicoPtr XZ  = Rico::Add(X, Z);
  RicoPtr Xp3 = Rico::Add(X, Rico::Integer(3));
  RicoPtr Xp7 = Rico::Add(Rico::Integer(7), X);

  // 5+7 = 32
  CHECK_EQUAL(Rico::Integer(12), SumDigester(A, B).get());
  // 0 + X          => X
  CHECK_EQUAL(X,                 SumDigester(Rico::Zero, X).get());
  // X + 0          => X
  CHECK_EQUAL(X,                 SumDigester(X, Rico::Zero).get());
  // 0 + 0          => 0
  CHECK_EQUAL(Rico::Zero,        SumDigester(Rico::Zero, Rico::Zero).get());
  // 5 + {-5}       => 0
  CHECK_EQUAL(Rico::Zero,        SumDigester(Rico::Integer(5), Rico::Integer(-5)).get());
  // 5 + X          => 5 + X
  CHECK_EQUAL(Rico::Add(A, X) ,  SumDigester(A, X).get());
  // X + 5          => 5 + X
  CHECK_EQUAL(Rico::Add(A, X) ,  SumDigester(X, A).get());
  // 3X + 5         => 5 + 3X
  CHECK_EQUAL(Rico::Add(A, X3),  SumDigester(X3, A).get());
  // X + Y          => X + Y
  CHECK_EQUAL(Rico::Add(X,Y),    SumDigester(X, Y).get());
  // (X+Y) + (X+Z)  => 2X + Y + Z
  CHECK_EQUAL(Rico::Add(Rico::Multiply(Rico::Two, X), Y, Z),  SumDigester(XY, XZ).get());
  // (X+3) + (X+7)  => 10 + 2X
  CHECK_EQUAL(Rico::Add(Rico::Integer(10), Rico::Multiply(Rico::Two, X)),    
	      SumDigester(Xp3, Xp7).get());
  // 3X + X + Y     => 4X + Y
  CHECK_EQUAL(Rico::Add(X4, Y),  SumDigester(Rico::Add(X3, X, Y)).get()); 
  // 3X + {-3}X     => 0
  CHECK_EQUAL(Rico::Zero, SumDigester(X3, X3n).get());
}

/////////////////////////////////////////////////////////////////////////////////////////
