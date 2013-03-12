#include <iostream>
#include <sstream>
#include <string>

#include "Rico.h"
#include "Util.h"

#include "RicoRegistry.h"

#include "RicoInteger.h"
#include "RicoDouble.h"
#include "RicoFraction.h"

#include "RicoNegate.h"
#include "RicoSqrt.h"
#include "RicoLog.h"
#include "RicoExp.h"

#include "RicoAdd.h"
#include "RicoSubtract.h"
#include "RicoMultiply.h"
#include "RicoDivide.h"
#include "RicoPower.h"

#include "RicoTriangular.h"
#include "RicoNormal.h"
#include "RicoLogNormal.h"
#include "RicoCauchy.h"
#include "RicoChiSquared.h"
#include "RicoStudents_t.h"
#include "RicoFisher.h"
#include "RicoBeta.h"
#include "RicoExponential.h"
#include "RicoLogistic.h"
#include "RicoWeibull.h"
#include "RicoUniform.h"

#include "RicoOperators.h"

#include "UnitTest++/UnitTest++.h"

// Rico ////////////////////////////////////////////////////////////////////////////
int Rico::registerRicoPtr(RicoPtr x) {return RicoRegistry::registerRicoPtr(x);}

RicoPtr Rico::NaN    = Rico::Fraction(0,0);
RicoPtr Rico::NegOne = Rico::Integer(-1);
RicoPtr Rico::Zero   = Rico::Integer(0);
RicoPtr Rico::Half   = Rico::Fraction(1,2);
RicoPtr Rico::One    = Rico::Integer(1);
RicoPtr Rico::Two    = Rico::Integer(2);
RicoPtr Rico::Three  = Rico::Integer(3);

RicoPtr Rico::Integer   (int x)                         {return RicoPtr(new RicoInteger(x));}
RicoPtr Rico::Double    ()                              {return RicoPtr(new RicoDouble());}
RicoPtr Rico::Double    (double x)                      {return RicoPtr(new RicoDouble(x));}
RicoPtr Rico::Fraction  (int x, int y)                  {return RicoPtr(new RicoFraction(x,y));}

// Neg(NaN)    => NaN
// Neg(Neg(x)) => A
// Neg(0)      => 0
// Neg(-n)     => n
// Neg(n)      => -n

RicoPtr Rico::Negate(RicoPtr x) {
  if(x->isNaN())    return Rico::NaN;
  if(x->isNegate()) return x->x();
  if(x->isZero())   return Rico::Zero;
  if(x->isNumber()) return ((const RicoNumber&)(*x)).getNeg();
  return RicoPtr(new RicoNegate(x));
}

RicoPtr Rico::Sqrt            (RicoPtr x)               {return RicoPtr(new RicoSqrt  (x));}
RicoPtr Rico::Log             (RicoPtr x)               {return RicoPtr(new RicoLog   (x));}
RicoPtr Rico::Exp             (RicoPtr x)               {return RicoPtr(new RicoExp   (x));}
RicoPtr Rico::Square          (RicoPtr x)               {return Power(x, Two);}

// 5 => {1/5}
// X => 1 / X

RicoPtr Rico::Inverse(RicoPtr x) {
  if(x->isNumber()) return ((const RicoNumber&)(*x)).getInverse();
  return Rico::Divide(Rico::One, x);
}

// (X+Y)+Z           => (X+Y+Z)
// X+(Y+Z)           => (X+Y+Z)
// (X+Y)+(Z+K)       => (X+Y+Z+K)
// (X+0)+Y+(Z+K+0)   => (X+Y+Z+K)  // filter out zero
// (X+(Y+(Z+0)+K+0)) => (X,Y,Z,K)  // Cascaded additions are flattened
// X+0               => X
// 0+Y               => Y
// X+NaN             => NaN
// NaN+Y             => NaN

RicoPtr Rico::Add(RicoPtr x, RicoPtr y) {
  vRicoPtr P; Add(P, x); Add(P, y); return Add(P);
}
RicoPtr Rico::Add(RicoPtr x, RicoPtr y, RicoPtr z) {
  vRicoPtr P; Add(P, x); Add(P, y); Add(P, z); return Add(P);
}
RicoPtr Rico::Add(RicoPtr x, RicoPtr y, RicoPtr z, RicoPtr k) {
  vRicoPtr P; Add(P, x); Add(P, y); Add(P, z); Add(P, k); return Add(P);
}
void Rico::Add(vRicoPtr& V, RicoPtr A) {
  if(A->isAdd()) {
    const vRicoPtr& P = ((const RicoFunction&)(*A)).get_params();
    for(vRicoPtr::const_iterator it = P.begin(); it != P.end(); it++) Add(V, *it);
  } else
    V.push_back(A);
}
// Final stage. Filters 0,NaN.
RicoPtr Rico::Add(const vRicoPtr& P) {
  vRicoPtr Q;
  for(vRicoPtr::const_iterator it = P.begin(); it != P.end(); it++) {
    if((*it)->isNaN()) return Rico::NaN;
    if((*it)->isZero()) continue;
    Q.push_back(*it);
  }
  if(Q.size() == 0) return Rico::Zero;
  if(Q.size() == 1) return Q.front();
  return RicoPtr(new RicoAdd(Q));
}

RicoPtr Rico::Subtract(RicoPtr x, RicoPtr y)            {return RicoPtr(new RicoSubtract(x,y));}
RicoPtr Rico::Subtract(const vRicoPtr& v)               {return RicoPtr(new RicoSubtract(v));}

// (X*Y)*Z           => (X*Y*Z)
// X*(Y*Z)           => (X*Y*Z)
// (X*Y)*(Z*K)       => (X*Y*Z*K)
// (X*1)*Y*(Z*K*1)   => (X*Y*Z*K)  // filter out one
// (X*(Y*(Z*1)*K*1)) => (X,Y,Z,K)  // Cascaded additions are flattened
// (X*(Y*(Z*0)*K*1)) => 0          // zero kills the expression

RicoPtr Rico::Multiply(RicoPtr x, RicoPtr y) {
  vRicoPtr P; Multiply(P, x); Multiply(P, y); return Multiply(P);
}
RicoPtr Rico::Multiply(RicoPtr x, RicoPtr y, RicoPtr z) {
  vRicoPtr P; Multiply(P, x); Multiply(P, y); Multiply(P, z); return Multiply(P);
}
RicoPtr Rico::Multiply(RicoPtr x, RicoPtr y, RicoPtr z, RicoPtr k) {
  vRicoPtr P; Multiply(P, x); Multiply(P, y); Multiply(P, z); Multiply(P, k); return Multiply(P);
}
void Rico::Multiply(vRicoPtr& V, RicoPtr A) {
  if(A->isMultiply()) {
    const vRicoPtr& P = ((const RicoFunction&)(*A)).get_params();
    for(vRicoPtr::const_iterator it = P.begin(); it != P.end(); it++) Multiply(V, *it);
  } else 
    V.push_back(A);
}
// Final stage. Filters 0,1,NaN.
RicoPtr Rico::Multiply(const vRicoPtr& P) {
  for(vRicoPtr::const_iterator it = P.begin(); it != P.end(); it++)
    if((*it)->isNaN())  return Rico::NaN; // must be in own loop.
  vRicoPtr Q;
  for(vRicoPtr::const_iterator it = P.begin(); it != P.end(); it++) {
    if((*it)->isZero()) return Rico::Zero;
    if((*it)->isOne()) continue;
    Q.push_back(*it);
  }
  if(Q.size() == 0) return Rico::One;
  if(Q.size() == 1) return Q.front();
  return RicoPtr(new RicoMultiply(Q));
}

// NaN / NaN => NaN
// NaN / y   => NaN
// x   / NaN => NaN
// x   / 1   => x
// n   / m   => <evaluate>
// 0   / y   => 0

RicoPtr Rico::Divide(RicoPtr x, RicoPtr y) {
  if(x->isNaN() || y->isNaN()) return Rico::NaN;
  if(y->isOne()) return x;
  if(x->isNumber() && y-> isNumber()) return x / y;
  if(x->isZero()) return Rico::Zero;
  return RicoPtr(new RicoDivide(x,y));
}

// NaN^A => NaN
// A^NaN => NaN
// A^0   => 1
// A^1   => A

RicoPtr Rico::Power(RicoPtr x, RicoPtr y) {
  if(x->isNaN() || y->isNaN()) return Rico::NaN;
  if(y->isOne())  return x;
  if(y->isZero()) return Rico::One;
  return RicoPtr(new RicoPower(x,y));
}

RicoPtr Rico::Function(const RicoFunction& f, const vRicoPtr& v) {
  return RicoPtr(RicoFunction::Function(f, v));
}

RicoPtr Rico::Triangular()                              {return RicoPtr(new RicoTriangular());}
RicoPtr Rico::Triangular(double M)                      {return RicoPtr(new RicoTriangular(M));}
RicoPtr Rico::Triangular(double L, double R)            {return RicoPtr(new RicoTriangular(L,R));}
RicoPtr Rico::Triangular(double L, double M, double R)  {return RicoPtr(new RicoTriangular(L,M,R));}

RicoPtr Rico::Normal()                                  {return RicoPtr(new RicoNormal());}
RicoPtr Rico::Normal(double mu)                         {return RicoPtr(new RicoNormal(mu));}
RicoPtr Rico::Normal(double mu, double sigma)           {return RicoPtr(new RicoNormal(mu, sigma));}

RicoPtr Rico::LogNormal()                               {return RicoPtr(new RicoLogNormal());}
RicoPtr Rico::LogNormal(double u)                       {return RicoPtr(new RicoLogNormal(u));}
RicoPtr Rico::LogNormal(double u, double s)             {return RicoPtr(new RicoLogNormal(u, s));}

RicoPtr Rico::Cauchy()                                  {return RicoPtr(new RicoCauchy());}
RicoPtr Rico::Cauchy(double u)                          {return RicoPtr(new RicoCauchy(u));}
RicoPtr Rico::Cauchy(double u, double s)                {return RicoPtr(new RicoCauchy(u, s));}

RicoPtr Rico::ChiSquared(int df)                        {return RicoPtr(new RicoChiSquared(df));}
RicoPtr Rico::Students_t(int df)                        {return RicoPtr(new RicoStudents_t(df));}
RicoPtr Rico::Fisher(int df1, int df2)                  {return RicoPtr(new RicoFisher(df1, df2));}
RicoPtr Rico::Beta(double a, double b)                  {return RicoPtr(new RicoBeta(a, b));}
RicoPtr Rico::Exponential(double L)                     {return RicoPtr(new RicoExponential(L));}
RicoPtr Rico::Logistic(double u, double s)              {return RicoPtr(new RicoLogistic(u, s));}
RicoPtr Rico::Weibull(double sh)                        {return RicoPtr(new RicoWeibull(sh));}
RicoPtr Rico::Weibull(double sh, double sc)             {return RicoPtr(new RicoWeibull(sh, sc));}
RicoPtr Rico::Uniform()                                 {return RicoPtr(new RicoUniform());}
RicoPtr Rico::Uniform(double L, double R)               {return RicoPtr(new RicoUniform(L, R));}

//////////////////////////////////////////////////////////////////////////////////

vector<double> Rico::CDF(vector<double> X) const {
  vector<double> P;
  for(int i = 0; i < X.size(); i++) P.push_back(CDF(X[i]));
  return P;
}

vector<double> Rico::Pr(vector<double> x) const {
  int N = x.size();
  vector<double> p;
  vector<double> P = CDF(x);
  vector<double> Q = Util::diff(P);
  p.push_back(P[0]);
  for(int i = 0; i < N-1; i++) p.push_back(Q[i]);
  p.push_back(1 - P[N-1]);
  return p;
}

vector<double> Rico::Density(vector<double> x) const {
  vector<double> z  = PointsToPartition(x);
  vector<double> dx = Util::diff(z);
  vector<double> p  = Pr(z);
  vector<double> q  = Util::Slice(p, 1, -1);
  vector<double> y  = Util::divide(q,dx);
  return y;
}

// Access ////////////////////////////////////////////////////////////////

bool Rico::isInteger () const {return isNumber  () && ((const RicoNumber&)(*this)).isInteger   ();}
bool Rico::isFraction() const {return isNumber  () && ((const RicoNumber&)(*this)).isFraction  ();}
bool Rico::isDouble  () const {return isNumber  () && ((const RicoNumber&)(*this)).isDouble    ();}
bool Rico::isZero    () const {return isNumber  () && ((const RicoNumber&)(*this)).isZero      ();}
bool Rico::isPos     () const {return isNumber  () && ((const RicoNumber&)(*this)).isPos       ();}
bool Rico::isNeg     () const {return isNumber  () && ((const RicoNumber&)(*this)).isNeg       ();}
bool Rico::isOne     () const {return isNumber  () && ((const RicoNumber&)(*this)).isOne       ();}
bool Rico::isNaN     () const {return isNumber  () && ((const RicoNumber&)(*this)).isNaN       ();}
bool Rico::isPosInf  () const {return isNumber  () && ((const RicoNumber&)(*this)).isPosInf    ();}
bool Rico::isNegInf  () const {return isNumber  () && ((const RicoNumber&)(*this)).isNegInf    ();}

bool Rico::isNegate  () const {return isFunction() && ((const RicoFunction&)(*this)).isNegate  ();}
bool Rico::isAdd     () const {return isFunction() && ((const RicoFunction&)(*this)).isAdd     ();}
bool Rico::isSubtract() const {return isFunction() && ((const RicoFunction&)(*this)).isSubtract();}
bool Rico::isMultiply() const {return isFunction() && ((const RicoFunction&)(*this)).isMultiply();}
bool Rico::isDivide  () const {return isFunction() && ((const RicoFunction&)(*this)).isDivide  ();}
bool Rico::isPower   () const {return isFunction() && ((const RicoFunction&)(*this)).isPower   ();}
bool Rico::isLog     () const {return isFunction() && ((const RicoFunction&)(*this)).isLog     ();}
bool Rico::isExp     () const {return isFunction() && ((const RicoFunction&)(*this)).isExp     ();}
bool Rico::isSqrt    () const {return isFunction() && ((const RicoFunction&)(*this)).isSqrt    ();}

RicoPtr Rico::x() const {
  if(isFunction()) return ((const RicoFunction&)(*this)).x();
  cout << "ERROR: Rico::x() invalid for non-function object" << endl;
  return Rico::NaN;
}

RicoPtr Rico::y() const {
  if(isFunction()) return ((const RicoFunction&)(*this)).y();
  cout << "ERROR: Rico::y() invalid for non-function object" << endl;
  return Rico::NaN;
}

// R/C interface //////////////////////////////////////////////////////////

string Rico::getDefinedExpression() const {
  stringstream sstr;
  sstr << getExpression() << endl;

  // Show definitions
  set<RicoID> rv_set;
  collectBasicRVs(rv_set);

  for(set<RicoID>::iterator ii = rv_set.begin(); ii != rv_set.end(); ++ii) {
    const RicoBasicRV& rv = RicoBasicRV::getBasicRV(*ii);
    sstr << "X" << *ii << " = " << rv.getDefinition() << endl; 
  }

  return sstr.str();
}

const vRicoPtr Rico::differentiate(RicoPtr This) const {
  set<RicoID> rv_set;
  collectBasicRVs(rv_set);

  vRicoPtr drv;
  for(set<RicoID>::iterator ii = rv_set.begin(); ii != rv_set.end(); ++ii) {
    drv.push_back(differentiate(This, *ii));
  }
  return drv;
}

RicoPtr Rico::differentiate(RicoPtr This, RicoPtr Wrt) const {
  // Translate generic RicoPtr to BasicRV id, if possible.
  if(!Wrt->isBasicRV()) {
    cout << "ERROR: Cannot differentiate with respect to non-BasicRV" << endl;
    return This;
  }
  return differentiate(This, id());
}

RicoPtr Rico::differentiate_single(RicoPtr This) const {
  set<RicoID> rv_set;
  collectBasicRVs(rv_set);

  if(rv_set.size() != 1) 
    cout << "WARNING: Rico::differentiate_single found " << rv_set.size() << " variables" << endl;

  set<RicoID>::iterator ii = rv_set.begin();  // choose first one only! Ignore others.
  return differentiate(This, *ii);
}

// Utility ///////////////////////////////////////////////////////////////

vector<double> Rico::PointsToPartition(vector<double> x) {
  int N = x.size();
  vector<double> y;
  if(N == 0) return y;
  if(N == 1) {
    y.push_back(x[0] - 1);
    y.push_back(x[0] + 1);
  } else {
    y.push_back(x[0] - (x[1] - x[0])/2);
    for(int i = 0; i < N-1; i++)
      y.push_back((x[i+1] + x[i])/2);
    y.push_back(x[N-1] + (x[N-1] - x[N-2])/2);
  }
  return y;
}

// Testing ///////////////////////////////////////////////////////////////

TEST(CheckWorksWithPointers)
{
    void* p = (void *)0x100;
    CHECK(p);
    CHECK(p != 0);
}

TEST(TestRicoNegate) {
  RicoPtr X = Rico::Normal();

  // Neg(NaN)    => NaN
  CHECK(Rico::Negate(Rico::NaN)->isNaN());

  // Neg(Neg(x)) => A
  CHECK_EQUAL(X, Rico::Negate(Rico::Negate(X)));

  // Neg(0)      => 0
  CHECK_EQUAL(Rico::Zero, Rico::Negate(Rico::Zero));

  // Neg(-n)     => n
  CHECK_EQUAL(Rico::Two, Rico::Negate(Rico::Integer(-2)));

  // Neg(n)     => -n
  CHECK_EQUAL(Rico::Integer(-2), Rico::Negate(Rico::Integer(2)));
}

TEST(TestRicoAdd) {
  RicoPtr X = Rico::Normal();
  RicoPtr Y = Rico::Cauchy();
  RicoPtr Z = Rico::Triangular();
  RicoPtr K = Rico::Integer(5);
  RicoPtr W = Rico::Integer(7);
  RicoPtr N = Rico::Zero;

  // (X+Y)+Z => (X+Y+Z)
  CHECK_EQUAL(Rico::Add(X,Y,Z), Rico::Add(Rico::Add(X,Y),Z));
 
 // X+(Y+Z) => (X+Y+Z)
  CHECK_EQUAL(Rico::Add(X,Y,Z), Rico::Add(X,Rico::Add(Y,Z)));

  // (X+Y)+(Z+K) => (X+Y+Z+K)
  CHECK_EQUAL(Rico::Add(X,Y,Z,K), Rico::Add(Rico::Add(X,Y),Rico::Add(Z,K)));

  // (X+0)+Y+(Z+K+0) => (X+Y+Z+K)  // filter out zero
  CHECK_EQUAL(Rico::Add(X,Y,Z,K), Rico::Add(Rico::Add(X,N), Y, Rico::Add(Z,K,N)));

  // (X+(Y+(Z+0)+K+0)) => (X,Y,Z,K)  // Cascaded additions are flattened
  CHECK_EQUAL(Rico::Add(X,Y,Z,K), Rico::Add(X,Rico::Add(Y,Rico::Add(Z,N),K,N)));

  // X+0               => X
  CHECK_EQUAL(X, Rico::Add(X, Rico::Zero));

  // 0+Y               => Y
  CHECK_EQUAL(Y, Rico::Add(Rico::Zero, Y));

  // X+NaN             => NaN
  CHECK(Rico::Add(X, Rico::NaN)->isNaN());

  // NaN+Y             => NaN
  CHECK(Rico::Add(Rico::NaN, Y)->isNaN());
}

TEST(TestRicoMultiply) {
  RicoPtr X = Rico::Normal();
  RicoPtr Y = Rico::Cauchy();
  RicoPtr Z = Rico::Triangular();
  RicoPtr K = Rico::Integer(5);
  RicoPtr W = Rico::Integer(7);
  RicoPtr N = Rico::One;
  RicoPtr M = Rico::Zero;

  // NaN*Y             &\rightarrow NaN
  CHECK(Rico::Multiply(Rico::NaN, Y)->isNaN());
  CHECK(Rico::Multiply(Rico::NaN, Rico::Zero)->isNaN());

  // X*NaN             &\rightarrow NaN
  CHECK(Rico::Multiply(X, Rico::NaN)->isNaN());
  CHECK(Rico::Multiply(Rico::Zero, Rico::NaN)->isNaN());

  // X*0               &\rightarrow 0
  CHECK_EQUAL(Rico::Zero, Rico::Multiply(X, Rico::Zero));

  // 0*Y               &\rightarrow 0
  CHECK_EQUAL(Rico::Zero, Rico::Multiply(Rico::Zero, Y));

  // X*1               &\rightarrow X
  CHECK_EQUAL(X, Rico::Multiply(X, Rico::One));

  // 1*Y               &\rightarrow Y
  CHECK_EQUAL(Y, Rico::Multiply(Rico::One, Y));

  // (X*Y)*Z => (X*Y*Z)
  CHECK_EQUAL(Rico::Multiply(X,Y,Z), Rico::Multiply(Rico::Multiply(X,Y),Z));

  // X*(Y*Z) => (X*Y*Z)
  CHECK_EQUAL(Rico::Multiply(X,Y,Z), Rico::Multiply(X,Rico::Multiply(Y,Z)));

  // (X*Y)*(Z*K) => (X*Y*Z*K)
  CHECK_EQUAL(Rico::Multiply(X,Y,Z,K), Rico::Multiply(Rico::Multiply(X,Y),Rico::Multiply(Z,K)));

  // (X*1)*Y*(Z*K*1) => (X*Y*Z*K)  // filter out Ones
  CHECK_EQUAL(Rico::Multiply(X,Y,Z,K), Rico::Multiply(Rico::Multiply(X,N),Y,Rico::Multiply(Z,K,N)));

  // (X*(Y*(Z*1)*K*1)) => (X,Y,Z,K)  // Cascades are flattened
  CHECK_EQUAL(Rico::Multiply(X,Y,Z,K), Rico::Multiply(X,Rico::Multiply(Y,Rico::Multiply(Z,N),K,N)));

  // (X*(Y*(Z*0)*K*1)) => 0          // zero kills the expression
  CHECK_EQUAL(Rico::Zero, Rico::Multiply(X,Rico::Multiply(Y,Rico::Multiply(Z,M),K,N)));
}

TEST(TestRicoDivide) {
  RicoPtr X = Rico::Normal();
  RicoPtr Y = Rico::Cauchy();

  // NaN / NaN => NaN
  CHECK(Rico::Divide(Rico::NaN, Rico::NaN)->isNaN());

  // NaN / y   => NaN
  CHECK(Rico::Divide(Rico::NaN, Y)->isNaN());

  // x   / NaN => NaN
  CHECK(Rico::Divide(X, Rico::NaN)->isNaN());

  // x   / 1   => x
  CHECK_EQUAL(X, Rico::Divide(X, Rico::One));

  // n   / m   => <evaluate>
  CHECK_EQUAL(Rico::Fraction(2,3), Rico::Divide(Rico::Two, Rico::Three));

  // 0   / y   => 0
  CHECK_EQUAL(Rico::Zero, Rico::Divide(Rico::Zero, Y));
}

TEST(TestRicoPower) {
  RicoPtr X = Rico::Normal();

  // X^NaN => NaN
  CHECK(Rico::Power(X, Rico::NaN)->isNaN());

  // NaN^X => NaN
  CHECK(Rico::Power(Rico::NaN, X)->isNaN());

  // X^0 => 1
  CHECK_EQUAL(Rico::One, Rico::Power(X, Rico::Zero));

  // X^1 => X
  CHECK_EQUAL(X, Rico::Power(X, Rico::One));
}

////////////////////////////////////////////////////////////////////////////////////////////////////
