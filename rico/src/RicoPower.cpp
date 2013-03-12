#include <iostream>
#include <sstream>
#include <string>

#include "RicoPower.h"
#include "RicoNumber.h"

// RicoPower ///////////////////////////////////////////////////////////////

string RicoPower::getExpression() const {
  stringstream sstr;
  sstr << x()->getExpression();
  sstr << " ^ ";
  sstr << y()->getExpression();
  return sstr.str();
}

RicoPtr RicoPower::differentiate(RicoPtr This, RicoID Wrt) const {
  // D(X^Y)  == D( exp(Y*ln(X)) ) = Y*X^(Y-1)*X' + X^Y*ln(X)*Y'

  const vRicoPtr& A  = get_params();
        vRicoPtr  B  = differentiate_params(Wrt);
  const RicoPtr&  x  = A.front();
  const RicoPtr&  y  = A.back();
  const RicoPtr&  xp = B.front();
  const RicoPtr&  yp = B.back();

  RicoPtr y1  = Rico::Subtract(y, Rico::One);
  RicoPtr xy1 = Rico::Power(x, y1);
  RicoPtr a   = Rico::Multiply(y, xy1, xp);
  RicoPtr lnx = Rico::Log(x);
  RicoPtr b   = Rico::Multiply(This, lnx, yp);
  return Rico::Add(a,b);
}

///////////////////////////////////////////////////////////////////////////////
