#include "RicoSubtract.h"

#include <iostream>
#include <sstream>
#include <string>

// RicoSubtract ///////////////////////////////////////////////////////////////

string RicoSubtract::getExpression() const {
  stringstream sstr;
  sstr << "(";
  sstr << x()->getExpression();
  sstr << " - ";
  sstr << y()->getExpression();
  sstr << ")";
  return sstr.str();
}

RicoPtr RicoSubtract::differentiate(RicoPtr This, RicoID Wrt) const {
  return Rico::Subtract(differentiate_params(Wrt));
}

///////////////////////////////////////////////////////////////////////////////
