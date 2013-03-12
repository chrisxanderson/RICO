#include "RicoLog.h"

#include <iostream>
#include <sstream>
#include <string>

// RicoLog ///////////////////////////////////////////////////////////////

string RicoLog::getExpression() const {
  stringstream sstr;
  sstr << "log(" << x()->getExpression() << ")";
  return sstr.str();
}

RicoPtr RicoLog::differentiate(RicoPtr This, RicoID Wrt) const {
  // D(log(X))  == x'/x
  const vRicoPtr& A = get_params();
        vRicoPtr  B = differentiate_params(Wrt);
  const RicoPtr&  x = A.front();
  const RicoPtr& xp = B.front();
  return Rico::Divide(xp,x);
}

///////////////////////////////////////////////////////////////////////////////
