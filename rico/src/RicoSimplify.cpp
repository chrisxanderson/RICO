#include <iostream>
#include <string>
#include <sstream>

#include "RicoSimplify.h"
#include "RicoAdd.h"
#include "RicoMultiply.h"
#include "RicoPower.h"
#include "RicoNumber.h"
#include "RicoInteger.h"
#include "RicoBasicRV.h"
#include "RicoOperators.h"

#include "SumProductMatrix.h"
#include "PowerOfSums.h"
#include "ProductDigester.h"
#include "SumDigester.h"

#include "UnitTest++/UnitTest++.h"

// RicoSimplify ////////////////////////////////////////////////////////////////////////////

RicoPtr RicoSimplify::Simplify(RicoPtr A) {
  if(!A->isFunction()) return A;

  const RicoFunction& F = (const RicoFunction&)(*A);
  const vRicoPtr&     P = F.get_params();

  vRicoPtr Q; // Depth-first
  for(vRicoPtr::const_iterator it = P.begin(); it != P.end(); it++) 
    Q.push_back(RicoSimplify::Simplify(*it));

  if(F.isNegate())   return RicoSimplify::Negate  (Q.front());
  if(F.isAdd())      return RicoSimplify::Add     (Q);
  if(F.isSubtract()) return RicoSimplify::Subtract(Q.front(), Q.back());
  if(F.isMultiply()) return RicoSimplify::Multiply(Q);
  if(F.isDivide())   return RicoSimplify::Divide  (Q.front(), Q.back());
  if(F.isPower())    return RicoSimplify::Power   (Q.front(), Q.back());
  if(F.isLog())      return RicoSimplify::Log     (Q.front());
  if(F.isExp())      return RicoSimplify::Exp     (Q.front());
  if(F.isSqrt())     return RicoSimplify::Sqrt    (Q.front());

  cout << "ERROR: RicoSimplify::Simplify does not support function type: " << F.type() << endl;
  return Rico::NaN;
}

// Simplification Rules are conditionally applied in order.

// -(NaN) => NaN
// -(n)   => (-n)
// -(n*A) => (-n)*A
//  -(A)  => -1 * A

RicoPtr RicoSimplify::Negate(RicoPtr A) {
  if(A->isNaN())    return Rico::NaN;
  if(A->isNegate()) return A->x();           // included for non-depth-first case
  if(A->isNumber()) return Rico::Negate(A);
  if(A->isMultiply() && A->x()->isNumber()) {
    const vRicoPtr& P = ((const RicoFunction&)(*A)).get_params();
    vRicoPtr Q;
    for(vRicoPtr::const_iterator it = P.begin(); it != P.end(); it++) {
      if(it == P.begin()) Q.push_back(Rico::Negate(*it));
      else                Q.push_back(*it);
    }
    return Rico::Multiply(Q);  // will naturally remove any 1's
  }
  return Rico::Multiply(Rico::NegOne, A);
}

// A + NaN   => NaN
// NaN + B   => NaN
// NaN + NaN => NaN
// 0 + 5     => 5
// 5 + 0     => 5
// 5 + 7     => 12
// 5 + (-5)  => 0
// A + 0     => A
// 0 + A     => A
// A + 5     => 5 + A                           // numbers first
// A + B     => A + B                           // respect order
// A + A     => 2A
// (A + B) + (A + C)     => 2A + B + C
// A + (B + C)           => A + B + C
// (A1+A2+A3) + (3A2+B2) => A1 + 4A2 + A3 + B2 // collect and flatten
// (5 + A) + (C + A)     => 5 + 2A + C         // TODO - numbers first.

RicoPtr RicoSimplify::Add(RicoPtr A, RicoPtr B) {
  if(A->isNaN() || B->isNaN()) return Rico::NaN;
  return SumDigester::Digest(A,B);
}

RicoPtr RicoSimplify::Add(const vRicoPtr& P) {
  for(vRicoPtr::const_iterator it = P.begin(); it != P.end(); it++)
    if((*it)->isNaN()) return Rico::NaN;
  return SumDigester::Digest(P); 
}

// A - NaN   => NaN
// NaN - B   => NaN
// NaN - NaN => NaN
// A - B     => A + SimpleNeg(B)

RicoPtr RicoSimplify::Subtract (RicoPtr A, RicoPtr B) {
  if(A->isNaN() || B->isNaN()) return Rico::NaN;
  return RicoSimplify::Add(A, RicoSimplify::Negate(B));
}

///////////////// The workhorse is ProductDigester. Just distributive law here. ////////
// Sum Handling
// A x (B + C)       => A x B + A x C           // recurse to RicoSimplify::Multiply and ::Add
// (A + B) x C       => A x C + B x C           // recurse to RicoSimplify::Multiply and ::Add

// How-to product handling: Use sparse two-row matrix. Allow for non-numeric exponents.
//                  a  b c                      // header of unique base components
//      a^2*b*c  -> 2  1 1                      // We can assume each component us unique (simplest)
//      a*b^{-1} -> 1 -1 o                      // some component may overlap with first row.
//
// __Preliminary___
// A * NaN   => NaN
// NaN * B   => NaN
// NaN * NaN => NaN
// __Product of Fractions____
// (A/B) * C     => (A*C) / B
// A * (B/C)     => (A*B) / C
// (A/B) * (C/D) => (A*C) / (B*D)
// __Product of Sums_____
// (A+B)*C   => A*C + B*C
// A*(B+C)   => A*B + A*C
// __Product of Products__
// (A*2*B)*(B*C) => 2AB^2C

RicoPtr RicoSimplify::Multiply(RicoPtr A, RicoPtr B) {
  vRicoPtr P; P.push_back(A); P.push_back(B); return RicoSimplify::Multiply(P);
}

RicoPtr RicoSimplify::Multiply(RicoPtr A, RicoPtr B, RicoPtr C) {
  vRicoPtr P; P.push_back(A); P.push_back(B); P.push_back(C); return RicoSimplify::Multiply(P);
}

RicoPtr RicoSimplify::Multiply(const vRicoPtr& P) {
  for(vRicoPtr::const_iterator it = P.begin(); it != P.end(); it++)
    if((*it)->isNaN()) return Rico::NaN;

  // Detect fractions
  vRicoPtr Num, Denom;
  for(vRicoPtr::const_iterator it = P.begin(); it != P.end(); it++) {
    if((*it)->isDivide()) {Num.push_back((*it)->x()); Denom.push_back((*it)->y());}
    else                   Num.push_back(*it);
  }
  if(Denom.size() > 0) 
    return RicoSimplify::Divide(RicoSimplify::Multiply(Num), 
				RicoSimplify::Multiply(Denom));
  
  // return on first sum.
  for(vRicoPtr::const_iterator it = P.begin(); it != P.end(); it++) {
    if((*it)->isAdd()) { 
      vRicoPtr SUM;
      const vRicoPtr& Q = ((const RicoFunction&)(**it)).get_params();
      for(vRicoPtr::const_iterator jt = Q.begin(); jt != Q.end(); jt++) {
	vRicoPtr PROD;
	for(vRicoPtr::const_iterator kt = P.begin(); kt != P.end(); kt++) {
	  if(kt == it) PROD.push_back(*jt);
	  else         PROD.push_back(*kt);
	}
	SUM.push_back(RicoSimplify::Multiply(PROD));
      }
      return RicoSimplify::Add(SUM);
    }
  }

  return ProductDigester::Digest(P); 

  //TODO: Before we return the Digester result we try to simplify 1D polynomial quotients, nD later.
}

// A   / NaN     => NaN
// NaN / B       => NaN
// NaN / NaN     => NaN
// n / m         => <evaluate>
// A / 1         => A
// A / (Bx/By)   => (A*By)/Bx
// A*B/A         => B

RicoPtr RicoSimplify::Divide(RicoPtr A, RicoPtr B) {
  if(A->isNaN() || B->isNaN()) return Rico::NaN;
  if(A->isNumber() && B->isNumber()) return A / B;
  if(B->isOne()) return A;
  if(B->isDivide()) return RicoSimplify::Divide(RicoSimplify::Multiply(A,B->y()), B->x());
  return ProductDigester::Digest(A,B); 
}

// NaN^B, A^NaN => NaN
// (A1^A2)^B    => A1^{A1*B}                // product of exponents
// exp(A1)^B    => exp(A1*B)                // product of exponents
// n1^n2        => <evaluate>               // See RicoNumber::Pow(n1,n2)
// A^0          => 1                        // A is non-numeric
// A^1          => A                        // A is non-numeric
// A^{-n}       => 1/A^n                    // avoid negative exponents
// A^{-n*C}     => 1/A^{n*C}
// (A1+A2+A3)^4 => A1^4 + ... + A3^4        // special algorithm for this complex case.
// (A1*A2*A3)^B => 2^B*A1^B*...*A2^B        // distributive law for exponents

RicoPtr RicoSimplify::Power(RicoPtr A, RicoPtr B) {
  if(A->isNaN() || B->isNaN()) return Rico::NaN;

  if(A->isPower())                   return Rico::Power(A->x(), RicoSimplify::Multiply(A->y(), B));
  if(A->isExp())                     return Rico::Exp(          RicoSimplify::Multiply(A->x(), B));
  if(A->isNumber() && B->isNumber()) return RicoNumber::Power(A,B); 
  if(B->isZero())                    return Rico::One;
  if(B->isOne())                     return A;

  if(B->isNeg()) return Rico::Divide(Rico::One, RicoSimplify::Power(A, Rico::Negate(B)));

  if(B->isMultiply() && B->x()->isNeg()) {
    const vRicoPtr& P = ((const RicoFunction&)(*B)).get_params();
    vRicoPtr Q;
    for(vRicoPtr::const_iterator it = P.begin(); it != P.end(); it++) {
      if(it == P.begin()) Q.push_back(Rico::Negate(*it));
      else                Q.push_back(*it);
    }
    return Rico::Divide(Rico::One, RicoSimplify::Power(A, Rico::Multiply(Q)));
  }

  // Special Case: (A1 + A2)^n  (n will not be positive)
  if(A->isAdd() && B->isInteger()) {
    const vRicoPtr& a = ((const RicoFunction&)(*A)).get_params();
    int             b = ((const RicoInteger &)(*B)).x();
    return PowerOfSums::Power(a, b);
  }

  if(A->isMultiply()) {
    vRicoPtr V;
    const vRicoPtr& P = ((const RicoFunction&)(*A)).get_params();
    for(vRicoPtr::const_iterator it = P.begin(); it != P.end(); it++) 
      V.push_back(RicoSimplify::Power(*it, B));
    return RicoSimplify::Multiply(V);
  }

  return Rico::Power(A, B);
}

// log(NaN)    => NaN
// log(exp(A)) => A
// log(A^B)    => B*log(A)
// log(A/B)    => log(A) + -1*log(B)
// log(A*B*C)  => log(A)+log(B)+log(C)
// log(1)      => 0                         
// log(n)      => do the math

RicoPtr RicoSimplify::Log(RicoPtr A) {
  if(A->isNaN()) return Rico::NaN;

  if(A->isExp()) return A->x();

  if(A->isPower()) return RicoSimplify::Multiply(A->y(), RicoSimplify::Log(A->x()));

  if(A->isDivide()) 
    return RicoSimplify::Add(RicoSimplify::Log(A->x()), 
			     Rico::Multiply(Rico::NegOne, RicoSimplify::Log(A->y())));

  if(A->isMultiply()) {
    const RicoFunction& F = (const RicoFunction&) *A;
    vRicoPtr P;
    for(vRicoPtr::const_iterator it = F.get_params().begin(); it != F.get_params().end(); it++) 
      P.push_back(RicoSimplify::Log(*it));
    return RicoSimplify::Add(P);
  }

  if(A->isNumber()) return RicoNumber::Log(A);

  return Rico::Log(A);
}

// exp(NaN)    => NaN
// exp(log(A)) => A
// exp(A+B+C)  => exp(A)*exp(B)*exp(C)
// exp(0)      => 1
// exp(n)      => do the math

RicoPtr RicoSimplify::Exp (RicoPtr A) {
  if(A->isNaN()) return Rico::NaN;

  if(A->isLog()) return A->x();

  if(A->isAdd()) {
    const RicoFunction& F = (const RicoFunction&) *A;
    vRicoPtr P;
    for(vRicoPtr::const_iterator it = F.get_params().begin(); it != F.get_params().end(); it++) 
      P.push_back(RicoSimplify::Exp(*it));
    return RicoSimplify::Multiply(P);
  }

  if(A->isNumber()) return RicoNumber::Exp(A);

  return Rico::Exp(A);
}

// sqrt(NaN) => NaN
// sqrt(A)   => A^{1/2}

RicoPtr RicoSimplify::Sqrt(RicoPtr A) {
  if(A->isNaN()) return Rico::NaN;
  return RicoSimplify::Power(A, Rico::Half);
}

/* Factor: 
 * A*B + A*C + D => A*(B+C) + D
 * L1*d + L1*L2*d^2 + L1*L2*L3*d^3  =>  L1*d(1+L2*d(1 + L3*d))
 *
 * 1) Load the full expression into the SumProductMatrix
 * 2) Scan columns for first most-common expression. Break ties by sum of powers.
 * 3) Bifurcate rows into A*B+AC and D based on presence or absense of common A.
 * 4) In "member" rows pull out as many A's as we can and decrement the matrix accordingly.
 * Example:
 *      Given A*B + A*C + D
 *           | A B C D |
 *	     | 1 1     |
 *	     | 1   1   |
 *	     |       1 |
 *	Since A dominates we have 2 rows and a singleton
 *           | A B C D |
 *	A *  | 0 1     | = A(B + C)
 *	     | 0   1   |
 *	+               
 *           | A B C D |
 *	     |       1 | = D
 *
 *      And the result is,
 *           A*(B+C) + D
 *
 */

RicoPtr RicoSimplify::Factor(RicoPtr A) {
  if(!A->isAdd()) return A;
  const RicoFunction& F = (const RicoFunction&) *A;

  vRicoPtr P;  // Depth-first
  for(vRicoPtr::const_iterator it = F.get_params().begin(); it != F.get_params().end(); it++) 
    P.push_back(Factor(*it));

  return SumProductMatrix(Rico::Add(P)).getFactored();
}

// Testing /////////////////////////////////////////////////////////////////////////////////

TEST(TestSimplifyNegate) {
  RicoPtr X = Rico::Normal();

  // -(NaN) -> NaN
  CHECK(RicoSimplify::Negate(Rico::NaN)->isNaN());

  // -(-X) => X
  CHECK_EQUAL(X, RicoSimplify::Negate(Rico::Negate(X)));

  // -{3/{-5}} => {3/5}
  CHECK_EQUAL(Rico::Fraction(3, 5), RicoSimplify::Negate(Rico::Fraction(3,-5)));

  // -{3/5} => {-3/5}
  CHECK_EQUAL(Rico::Fraction(-3, 5), RicoSimplify::Negate(Rico::Fraction(3,5)));

  // -{1} => {-1}
  CHECK_EQUAL(Rico::NegOne, RicoSimplify::Negate(Rico::One));

  // -(n*X)  => {-n}*X
  CHECK_EQUAL(Rico::Multiply(Rico::Two, X), 
	      RicoSimplify::Negate(Rico::Multiply(Rico::Integer(-2),X)));

  // -(-1*X)  => X
  CHECK_EQUAL(X, RicoSimplify::Negate(Rico::Multiply(Rico::NegOne,X)));

  //  -{X}  => -1 * X
  CHECK_EQUAL(Rico::Multiply(Rico::NegOne, X), RicoSimplify::Negate(X));
}

TEST(TestSimplifyLog) {
  RicoPtr X  = Rico::Normal();
  RicoPtr Y  = Rico::Cauchy();

  // Log(Exp(X)) => X
  CHECK_EQUAL(X, RicoSimplify::Log(Rico::Exp(X)));

  // Log(X^Y) => Y*log(X)
  CHECK_EQUAL(Rico::Multiply(Y, Rico::Log(X)), RicoSimplify::Log(Rico::Power(X,Y)));

  // Log(X/Y) => log(X) + -1*log(Y)
  CHECK_EQUAL(Rico::Add(Rico::Log(X), Rico::Multiply(Rico::NegOne, Rico::Log(Y))), 
	      RicoSimplify::Log(Rico::Divide(X,Y)));

  // Log(XY) => Log(X) + Log(Y)
  CHECK_EQUAL(Rico::Add(Rico::Log(X), Rico::Log(Y)), RicoSimplify::Log(Rico::Multiply(X,Y)));

  // Log(Exp(X)*3*Y) => X + Log(3) + Log(Y)
  CHECK_EQUAL(Rico::Add(X, RicoNumber::Log(Rico::Three), Rico::Log(Y)), 
	      RicoSimplify::Log(Rico::Multiply(Rico::Exp(X), Rico::Three, Y)));

  // Log(1) => 0
  CHECK_EQUAL(Rico::Zero, RicoSimplify::Log(Rico::One));

  // Log(7.2) => 1.9740810260220096270241953835352169059203800300974916
  CHECK_EQUAL(Rico::Double(1.9740810260220096270241953835352169059203800300974916), 
	      RicoSimplify::Log(Rico::Double(7.2)));
}

TEST(TestSimplifyExp) {
  RicoPtr X  = Rico::Normal();
  RicoPtr Y  = Rico::Cauchy();

  // Exp(Log(X)) => X
  CHECK_EQUAL(X, RicoSimplify::Exp(Rico::Log(X)));

  // Exp(X+Y) => Exp(X)*Exp(Y)
  CHECK_EQUAL(Rico::Multiply(Rico::Exp(X), Rico::Exp(Y)), RicoSimplify::Exp(Rico::Add(X,Y)));

  // Exp(Log(X) + 3 + Y) => X * Exp(3) * Exp(Y)
  CHECK_EQUAL(Rico::Multiply(X, RicoNumber::Exp(Rico::Three), Rico::Exp(Y)), 
	      RicoSimplify::Exp(Rico::Add(Rico::Log(X), Rico::Three, Y)));

  // Exp(0) => 1
  CHECK_EQUAL(Rico::One, RicoSimplify::Exp(Rico::Zero));

  // Exp(7.2) => 1339.4307643944178296873515152987188838689536926177469
  CHECK_EQUAL(Rico::Double(1339.4307643944178296873515152987188838689536926177469), 
	      RicoSimplify::Exp(Rico::Double(7.2)));
}


TEST(TestSimplifyPower) {
  RicoPtr X = Rico::Normal();
  RicoPtr Y = Rico::Cauchy();
  RicoPtr Z = Rico::Triangular();

  // NaN ^ NaN => NaN
  CHECK(RicoSimplify::Power(Rico::NaN, Rico::NaN)->isNaN());

  // X ^ NaN => NaN
  CHECK(RicoSimplify::Power(X, Rico::NaN)->isNaN());

  // NaN ^ Y => NaN
  CHECK(RicoSimplify::Power(Rico::NaN, Y)->isNaN());

  // (X+Y+Z)^2 => X^2 + ... + Z^2  // special algorithm for this complex case.
  RicoPtr XYZ = Rico::Add(X, Y, Z);
  RicoPtr X2  = Rico::Power(X, Rico::Two);
  RicoPtr Y2  = Rico::Power(Y, Rico::Two);
  RicoPtr Z2  = Rico::Power(Z, Rico::Two);
  RicoPtr XY  = Rico::Multiply(Rico::Two, X, Y);
  RicoPtr XZ  = Rico::Multiply(Rico::Two, X, Z);
  RicoPtr YZ  = Rico::Multiply(Rico::Two, Y, Z);
  vRicoPtr P; P.push_back(X2); P.push_back(XY); P.push_back(XZ);
  P.push_back(Y2); P.push_back(YZ); P.push_back(Z2);
  CHECK_EQUAL(Rico::Add(P), RicoSimplify::Power(XYZ, Rico::Two));

  // (XY)^Z => X^Z * Y^Z
  CHECK_EQUAL(Rico::Multiply(Rico::Power(X,Z), Rico::Power(Y,Z)), 
	      RicoSimplify::Power(Rico::Multiply(X,Y),Z)); 

  // (X^Y)^Z => X^(YZ)
  CHECK_EQUAL(Rico::Power(X, Rico::Multiply(Y,Z)), RicoSimplify::Power(Rico::Power(X,Y),Z));

  // exp(X)^Y => exp(XY)          // product of exponents
  CHECK_EQUAL(Rico::Exp(Rico::Multiply(X,Y)), RicoSimplify::Power(Rico::Exp(X),Y));

  // n1^n2 => <evaluate>         // See RicoNumber::Pow(n1,n2)
  CHECK_EQUAL(Rico::Integer(9), RicoSimplify::Power(Rico::Three, Rico::Two));

  // X^0 => 1
  CHECK_EQUAL(Rico::One, RicoSimplify::Power(X, Rico::Zero));

  // X^1 => X
  CHECK_EQUAL(X, RicoSimplify::Power(X, Rico::One));

  // X^{-1} => 1/X
  CHECK_EQUAL(Rico::Divide(Rico::One,X), RicoSimplify::Power(X, Rico::NegOne));

  // X^{-2} => 1/X^2
  CHECK_EQUAL(Rico::Divide(Rico::One,Rico::Square(X)), RicoSimplify::Power(X, Rico::Integer(-2)));

  // X^{-n*Y}     => 1/X^{n*Y}
  CHECK_EQUAL(Rico::Divide(Rico::One, Rico::Power(X, Rico::Multiply(Rico::Two, Y))),
	      RicoSimplify::Power(X, Rico::Multiply(Rico::Integer(-2), Y)));
}

TEST(TestSimplifyAdd) {
  RicoPtr X = Rico::Normal();
  RicoPtr Y = Rico::Cauchy();
  RicoPtr Z = Rico::Triangular();

  RicoPtr F = Rico::Integer( 5);
  RicoPtr M = Rico::Integer(-5);
  RicoPtr S = Rico::Integer( 7);
  RicoPtr O = Rico::Zero;

  CHECK_EQUAL(F,                 RicoSimplify::Add(F, O));
  CHECK_EQUAL(F,                 RicoSimplify::Add(O, F));
  CHECK_EQUAL(Rico::Integer(12), RicoSimplify::Add(F, S));
  CHECK_EQUAL(Rico::Zero,        RicoSimplify::Add(F, M));
  CHECK_EQUAL(Rico::Integer(5),  RicoSimplify::Add(X, F)->x());
  CHECK_EQUAL(X,                 RicoSimplify::Add(X, Y)->x());
  CHECK_EQUAL(X,                 RicoSimplify::Add(X, O));

  RicoPtr XY = Rico::Add(X,Y);
  RicoPtr XZ = Rico::Add(X,Z);
  RicoPtr YZ = Rico::Add(Y,Z);
  CHECK_EQUAL(XY,                                            RicoSimplify::Add(X, Y));
  CHECK_EQUAL(Rico::Multiply(Rico::Two, X),                  RicoSimplify::Add(X, X));
  CHECK_EQUAL(Rico::Add(Rico::Multiply(Rico::Two, X), Y, Z), RicoSimplify::Add(XY,XZ));
  CHECK_EQUAL(Rico::Add(X,Y,Z),                              RicoSimplify::Add(XY, Z));
  CHECK_EQUAL(Rico::Add(X,Y,Z),                              RicoSimplify::Add(X, YZ));

  RicoPtr Y3Z  = Rico::Add(Rico::Multiply(Rico::Three, Y), Z);
  RicoPtr XY4Z = Rico::Add(X, Rico::Multiply(Rico::Integer(4), Y), Z);
  CHECK_EQUAL(XY4Z, RicoSimplify::Add(XY, Y3Z));

  // (5+X)+(Y+X) => 5 + 2X + Y         // numbers first.
  CHECK_EQUAL(Rico::Add(F, Rico::Multiply(Rico::Two,X),Y),
	      RicoSimplify::Add(Rico::Add(F,X),Rico::Add(Y,X)));
}

TEST(TestSimplifyDivide) {
  RicoPtr X  = Rico::Normal();
  RicoPtr Y  = Rico::Cauchy();
  RicoPtr Z  = Rico::Triangular();

  // X   / NaN    => NaN
  CHECK(RicoSimplify::Divide(X,Rico::NaN)->isNaN());

  // NaN / Y      => NaN
  CHECK(RicoSimplify::Divide(Rico::NaN,Y)->isNaN());
  
  // NaN / NaN    => NaN
  CHECK(RicoSimplify::Divide(Rico::NaN,Rico::NaN)->isNaN());

  // X / 1        => X
  CHECK_EQUAL(X, RicoSimplify::Divide(X,Rico::One));

  // X / Y        => X / Y
  CHECK_EQUAL(Rico::Divide(X,Y), RicoSimplify::Divide(X,Y));

  // X / (Y / Z)  => (XZ) / Y
  CHECK_EQUAL(Rico::Divide(Rico::Multiply(X,Z),Y), RicoSimplify::Divide(X,Rico::Divide(Y,Z)));
}

TEST(TestSimplifyMultiply) {
  RicoPtr X  = Rico::Normal();
  RicoPtr Y  = Rico::Cauchy();
  RicoPtr Z  = Rico::Triangular();
  RicoPtr XY = Rico::Add(X,Y);

  RicoPtr F = Rico::Integer(5);
  RicoPtr M = Rico::Fraction(1,5);
  RicoPtr S = Rico::Integer(7);

  CHECK_EQUAL(Rico::Zero,        RicoSimplify::Multiply(Rico::Zero, F));
  CHECK_EQUAL(F,                 RicoSimplify::Multiply(Rico::One, F));
  CHECK_EQUAL(Rico::Zero,        RicoSimplify::Multiply(F, Rico::Zero));
  CHECK_EQUAL(Rico::Integer(35), RicoSimplify::Multiply(F, S));
  CHECK_EQUAL(Rico::One,         RicoSimplify::Multiply(F, M));
  CHECK_EQUAL(F,                 RicoSimplify::Multiply(F, Rico::One));
  CHECK_EQUAL(X,                 RicoSimplify::Multiply(X, Rico::One));
  CHECK_EQUAL(Rico::Zero,        RicoSimplify::Multiply(X, Rico::Zero));
  CHECK_EQUAL(F,                 RicoSimplify::Multiply(X, F)->x()); // Number first
  CHECK_EQUAL(X,                 RicoSimplify::Multiply(X, Y)->x()); // Preserve order

  // X*(Y+1) = XY + X
  CHECK_EQUAL(Rico::Add(Rico::Multiply(X,Y), X), RicoSimplify::Multiply(X, Rico::Add(Y,Rico::One)));

  // (X + Y)*(X + Y) = X^2 + 2XY + Y^2
  RicoPtr XY2 = Rico::Add(Rico::Square(X), Rico::Multiply(Rico::Two, X, Y), Rico::Square(Y));
  CHECK_EQUAL(XY2, RicoSimplify::Multiply(XY,XY));

  // X*(Y+1)*Z => XYZ + XZ
  CHECK_EQUAL(Rico::Add(Rico::Multiply(X,Y,Z), Rico::Multiply(X,Z)),
	      RicoSimplify::Multiply(X, Rico::Add(Y,Rico::One), Z));

  // (X/Y) * Z => (XZ) / Y
  CHECK_EQUAL(Rico::Divide(Rico::Multiply(X,Z),Y), RicoSimplify::Multiply(Rico::Divide(X,Y),Z));

  // X * (Y/Z) => (XY) / Z
  CHECK_EQUAL(Rico::Divide(Rico::Multiply(X,Y),Z), RicoSimplify::Multiply(X,Rico::Divide(Y,Z)));

  // (X/Y) * (Z/F) => {1/F}XZ / Y
  CHECK_EQUAL(Rico::Divide(Rico::Multiply(Rico::Inverse(F),X,Z), Y),
	      RicoSimplify::Multiply(Rico::Divide(X,Y), Rico::Divide(Z,F)));

  // (X/Y) * (Y/Z) => X / Z
  CHECK_EQUAL(Rico::Divide(X,Z), RicoSimplify::Multiply(Rico::Divide(X,Y), Rico::Divide(Y,Z)));

  // (X+Y)*(X+Y)*(X+Y) => (X+Y)^3
  CHECK_EQUAL(RicoSimplify::Power(XY,Rico::Three), RicoSimplify::Multiply(XY,XY,XY));
}

TEST(TestSimplify) {
  RicoPtr X = Rico::Normal();
  RicoPtr Y = Rico::Cauchy();
  RicoPtr Z = Rico::Triangular();

  // (X+1)^2 / (X+1) => X+1   // won't work, yet. Need polynomial reduction.
  //RicoPtr X1 = Rico::Add(X,Rico::One);
  //CHECK_EQUAL(X1, RicoSimplify::Simplify(Rico::Divide(Rico::Square(X1),X1)));

  // ((X+Y)Z)^2   => X^2Z^2 +2XYZ^2 + Y^2Z^2
  CHECK_EQUAL(Rico::Add(Rico::Multiply(Rico::Square(X),Rico::Square(Z)), 
			Rico::Multiply(Rico::Two,Rico::Multiply(X,Y,Rico::Square(Z))),
			Rico::Multiply(Rico::Square(Y),Rico::Square(Z))),
	      RicoSimplify::Simplify(Rico::Square(Rico::Multiply(Rico::Add(X,Y),Z))));
}

////////////////////////////////////////////////////////////////////////////////////////////
