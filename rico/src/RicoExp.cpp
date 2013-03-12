#include "RicoExp.h"

#include <iostream>
#include <sstream>
#include <string>

// RicoExp ///////////////////////////////////////////////////////////////

string RicoExp::getExpression() const {
  stringstream sstr;
  sstr << "exp(" << x()->getExpression() << ")";
  return sstr.str();
}

RicoPtr RicoExp::differentiate(RicoPtr This, RicoID Wrt) const {
  // D(exp(X))  == exp(x)*x'
  const vRicoPtr& A = get_params();
        vRicoPtr  B = differentiate_params(Wrt);
  const RicoPtr&  x = A.front();
  const RicoPtr& xp = B.front();
  return Rico::Multiply(This, xp);
}

///////////////////////////////////////////////////////////////////////////////
