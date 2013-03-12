#ifndef __RICO_EXP_H__
#define __RICO_EXP_H__

#include "RicoFunction.h"

/////////////////////////////////////////////////////////////////////////

class RicoExp : public RicoFunction {

 public:
  RicoExp(RicoPtr x) : RicoFunction(RicoFunction::EXP, x) {} 
  ~RicoExp() {}

 public: // Interface
  virtual string  getExpression()                         const;
  virtual RicoPtr differentiate(RicoPtr This, RicoID Wrt) const;

}; // RicoExp

/////////////////////////////////////////////////////////////////////////
#endif // __RICO_EXP_H__

