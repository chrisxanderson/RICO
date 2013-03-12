#ifndef __RICO_MULTIPLY_H__
#define __RICO_MULTIPLY_H__

#include "RicoFunction.h"

/////////////////////////////////////////////////////////////////////////

class RicoMultiply : public RicoFunction {
 public:
 RicoMultiply(RicoPtr x, RicoPtr y)            : RicoFunction(RicoFunction::MUL,x,y)   {}
 RicoMultiply(RicoPtr x, RicoPtr y, RicoPtr z) : RicoFunction(RicoFunction::MUL,x,y,z) {}
 RicoMultiply(const vRicoPtr& v)               : RicoFunction(RicoFunction::MUL,v)     {}
  ~RicoMultiply() {}

 public: // Interface
  virtual string  getExpression()                         const;
  virtual RicoPtr differentiate(RicoPtr This, RicoID Wrt) const;

  static RicoPtr Simple(RicoPtr x, RicoPtr y);
  static RicoPtr Simple(const vRicoPtr& v);

 public: // Operators
  virtual bool operator==(RicoPtr That) const;

}; // RicoMultiply

/////////////////////////////////////////////////////////////////////////
#endif //__RICO_MULTIPLY_H__

