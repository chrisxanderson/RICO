#include "RicoNegate.h"

#include <iostream>
#include <sstream>
#include <string>

// RicoNegate ///////////////////////////////////////////////////////////////

string RicoNegate::getExpression() const {
  stringstream sstr;
  sstr << "(";
  sstr << "-";
  sstr << x()->getExpression();
  sstr << ")";
  return sstr.str();
}

RicoPtr RicoNegate::differentiate(RicoPtr This, RicoID Wrt) const {
  // D(-x)  == -x'
  vRicoPtr A = differentiate_params(Wrt);
  return Rico::Negate(A.front());
}
