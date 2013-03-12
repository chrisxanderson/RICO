#include "RicoAdd.h"

#include <iostream>
#include <sstream>
#include <string>

// RicoAdd ///////////////////////////////////////////////////////////////

string RicoAdd::getExpression() const {
  stringstream sstr;

  sstr << "(";

  for(vRicoPtr::const_iterator it = get_params().begin(); it != get_params().end(); it++) {
    vRicoPtr::const_iterator jt = it;
    if(++jt != get_params().end()) {
      sstr << (*it)->getExpression();
      sstr << " + ";
    } else 
      sstr << (*it)->getExpression();
  }

  sstr << ")";

  return sstr.str();
}

RicoPtr RicoAdd::differentiate(RicoPtr This, RicoID Wrt) const {
  return Rico::Add(differentiate_params(Wrt));
}

bool RicoAdd::operator==(RicoPtr That) const {
  if(That->type() != type()) return false;
  const RicoFunction& a = (const RicoFunction&) *this;
  const RicoFunction& b = (const RicoFunction&) *That;
  return a == b;
}


///////////////////////////////////////////////////////////////////////////////
