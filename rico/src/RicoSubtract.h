#ifndef __RICO_SUBTRACT_H__
#define __RICO_SUBTRACT_H__

#include "RicoFunction.h"

/////////////////////////////////////////////////////////////////////////

class RicoSubtract : public RicoFunction {

 public:
 RicoSubtract(RicoPtr x, RicoPtr y) : RicoFunction(RicoFunction::SUB, x, y) {} 
 RicoSubtract(const vRicoPtr& v)    : RicoFunction(RicoFunction::SUB, v)    {}
  ~RicoSubtract() {}

 public: // Interface
  virtual string  getExpression()                         const;
  virtual RicoPtr differentiate(RicoPtr This, RicoID Wrt) const;

}; // RicoSubtract

/////////////////////////////////////////////////////////////////////////
#endif // __RICO_SUBTRACT_H__

