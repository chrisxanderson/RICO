#include "RicoRegistry.h"
#include "RicoBasicRV.h"
#include "RicoInteger.h"

#include <iostream>
#include <sstream>
#include <string>

// #include <utility> // why is this here?

using namespace std;

// RicoBasicRV ///////////////////////////////////////////////////////////////

RicoID RicoBasicRV::rico_id = 1;

map< RicoID, const RicoBasicRV* > * RicoBasicRV::basic_rv_map = 
  new map< RicoID, const RicoBasicRV* >();

// Automatically set the id and register self locally.
RicoBasicRV::RicoBasicRV() : Rico(Rico::BASIC_RV), m_id(rico_id++) {
  (*basic_rv_map)[m_id] = this;
}

string RicoBasicRV::getExpression() const {
  stringstream sstr;
  sstr << "X" << id();  // All BasicRV's issued and id on create.
  return sstr.str();
}

RicoPtr RicoBasicRV::differentiate(RicoPtr This, RicoID Wrt) const {
  if(id() == Wrt) return Rico::One;
  else            return Rico::Zero;
}

bool RicoBasicRV::operator==(RicoPtr Y) const {
  if(type() != Y->type()) return false;
  return id() == Y->id();
}

///////////////////////////////////////////////////////////////////////////////
