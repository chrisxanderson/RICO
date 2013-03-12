#include "RicoMultiply.h"

#include <iostream>
#include <sstream>
#include <string>

// RicoMultiply ///////////////////////////////////////////////////////////////

RicoPtr RicoMultiply::Simple(RicoPtr x, RicoPtr y) {
  vRicoPtr v;
  v.push_back(x);
  v.push_back(y);
  return Simple(v);
}

RicoPtr RicoMultiply::Simple(const vRicoPtr& v) {
  int N = 0;
  vRicoPtr x;
  for(vRicoPtr::const_iterator it = v.begin(); it != v.end(); it++) {
    if((*it)->isZero()) return Rico::Zero;
    if((*it)->isOne()) continue;
    if((*it)->isNaN()) return *it;
    x.push_back(*it);
    N++;
  }
  if(N == 0) return Rico::One;
  if(N == 1) return x.front();
  return Rico::Multiply(x);
}

string RicoMultiply::getExpression() const {
  stringstream sstr;

  sstr << "(";

  for(vRicoPtr::const_iterator it = get_params().begin(); it != get_params().end(); it++) {
    vRicoPtr::const_iterator jt = it;
    if(++jt != get_params().end()) {
      sstr << (*it)->getExpression();
      sstr << " * ";
    } else 
      sstr << (*it)->getExpression();
  }

  sstr << ")";

  return sstr.str();
}

RicoPtr RicoMultiply::differentiate(RicoPtr This, RicoID Wrt) const {
  // D(xy)  = xy' + x'y
  // D(xyz) = x(yz)' + x'(yz) = x(yz' + y'z) + x'yz = xyz' + xy'z + x'yz
  const vRicoPtr& A = get_params();
        vRicoPtr  B = differentiate_params(Wrt);
        vRicoPtr  C;
	vRicoPtr  T = A;

 vRicoPtr::iterator       it;
 vRicoPtr::const_iterator jt;
 for(it = T.begin(), jt = B.begin(); it != T.end(); it++, jt++) {
   RicoPtr tmp = *it;
   *it = *jt;
   C.push_back(Rico::Multiply(T));
   *it = tmp;
  }
 
 return Rico::Add(C);
}

bool RicoMultiply::operator==(RicoPtr That) const {
  if(That->type() != type()) return false;
  const RicoFunction& a = (const RicoFunction&) *this;
  const RicoFunction& b = (const RicoFunction&) *That;
  return a == b;
}

///////////////////////////////////////////////////////////////////////////////
