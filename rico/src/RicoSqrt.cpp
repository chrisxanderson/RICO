#include <iostream>
#include <sstream>
#include <string>

#include "RicoSqrt.h"

// RicoSqrt ///////////////////////////////////////////////////////////////

string RicoSqrt::getExpression() const {
  stringstream sstr;
  sstr << "sqrt(" << x()->getExpression() << ")";
  return sstr.str();
}

RicoPtr RicoSqrt::differentiate(RicoPtr This, RicoID Wrt) const {
  // D(sqrt(X))  == 1/2*D(X)/sqrt(X)

  const vRicoPtr& A  = get_params();
        vRicoPtr  B  = differentiate_params(Wrt);
  const RicoPtr&  x  = A.front();
  const RicoPtr&  xp = B.front();

  RicoPtr  a  = Rico::Fraction(1,2);
  RicoPtr  b  = Rico::Divide(xp, This);
  return Rico::Multiply(a,b);
}

///////////////////////////////////////////////////////////////////////////////
