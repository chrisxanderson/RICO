#include "RicoFunction.h"
#include "RicoBasicRV.h"

#include "RicoNegate.h"
#include "RicoAdd.h"
#include "RicoSubtract.h"
#include "RicoMultiply.h"
#include "RicoDivide.h"
#include "RicoPower.h"
#include "RicoLog.h"
#include "RicoExp.h"
#include "RicoSqrt.h"
#include "RicoOperators.h"

#include <iostream>
#include <sstream>
#include <string>

// RicoFunction ///////////////////////////////////////////////////////////////

void RicoFunction::params(vRicoPtr& P) const {
  for(vRicoPtr::const_iterator it = get_params().begin(); it != get_params().end(); it++) 
    P.push_back(*it);
}

void RicoFunction::params(vRicoPtr& P, int i) const {
  int I = 0;
  for(vRicoPtr::const_iterator it = get_params().begin(); it != get_params().end(); it++, I++) 
    if(I != i) P.push_back(*it);
}

RicoFunction* RicoFunction::Function(const RicoFunction& f, const vRicoPtr& v) {
  switch(f.ftype()) {
  case NEG : return new RicoNegate  (v.front());
  case ADD : return new RicoAdd     (v);
  case SUB : return new RicoSubtract(v.front(), v.back());
  case MUL : return new RicoMultiply(v);
  case DIV : return new RicoDivide  (v.front(), v.back());
  case POW : return new RicoPower   (v.front(), v.back());
  case LOG : return new RicoLog     (v.front());
  case EXP : return new RicoExp     (v.front());
  case SQRT: return new RicoSqrt    (v.front());
  }
  cout << "ERROR: unsupported function type: " << f.ftype() << endl;
  return new RicoAdd(v); // just something to do.
}

// Utility ////////////////////////////////////////////////////////////////////

void RicoFunction::collectBasicRVs(set<RicoID>& rv_set) const {
  for(vRicoPtr::const_iterator it = get_params().begin(); it != get_params().end(); it++) 
    (*it)->collectBasicRVs(rv_set);
}

vRicoPtr RicoFunction::differentiate_params(RicoID Wrt) const  {
  vRicoPtr B;
  for(vRicoPtr::const_iterator it = get_params().begin(); it != get_params().end(); it++) {
    RicoPtr d = (*it)->differentiate(*it, Wrt);
    B.push_back(d);
  }
  return B;
}

bool RicoFunction::operator==(RicoPtr That) const {
  // Assumes if N = 1 or 2 then order matters. If N > 2 then error to call here.
  // this is overridden by Add and Multiply and comes back to specific operator==
  if(That->type() != type()) return false;
  const RicoFunction& that = (const RicoFunction&) *That;

  if(size() == 1) return x() == that.x();
  if(size() == 2) return x() == that.x() && y() == that.y();

  cout << "RicoFunction::operator==() called with " << size() << " > 2 operands." << endl;

  return false;
}

bool RicoFunction::operator==(const RicoFunction& that) const {
  // assumes order doesn't matter.
  
  const vRicoPtr& A = get_params();

  if(size() != that.size()) return false;

  vRicoPtr B = that.get_params(); // consumable list

  for(vRicoPtr::const_iterator it = A.begin(); it != A.end(); it++)
    if(!RicoFunction::test_membership(*it, B)) return false;

  return true;
}

bool RicoFunction::test_membership(RicoPtr ai, vRicoPtr& B) {
  for(vRicoPtr::iterator it = B.begin(); it != B.end(); it++) {
    if(ai == *it) {   // be sure to access operator==(RicoPtr, RicoPtr)
      B.erase(it);
      return true;
    }
  }
  return false;
}

///////////////////////////////////////////////////////////////////////////////
