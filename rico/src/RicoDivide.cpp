#include "RicoDivide.h"

#include <iostream>
#include <sstream>
#include <string>

// RicoDivide ///////////////////////////////////////////////////////////////

RicoDivide::~RicoDivide() {}

string RicoDivide::getExpression() const {
  stringstream sstr;
  sstr << "(";
  sstr << x()->getExpression();
  sstr << " / ";
  sstr << y()->getExpression();
  sstr << ")";
  return sstr.str();
}

RicoPtr RicoDivide::differentiate(RicoPtr This, RicoID Wrt) const {
  // D(x/y)  = x'/y-xy'/y^2
  const vRicoPtr& A  = get_params();
        vRicoPtr  B  = differentiate_params(Wrt);
  const RicoPtr&  x  = A.front();
  const RicoPtr&  y  = A.back();
  const RicoPtr&  xp = B.front();
  const RicoPtr&  yp = B.back();
  RicoPtr         c  = Rico::Divide(xp,y);
  RicoPtr         d  = Rico::Divide( Rico::Multiply(x,yp), Rico::Square(y) );
  return Rico::Subtract(c,d);
}

///////////////////////////////////////////////////////////////////////////////
