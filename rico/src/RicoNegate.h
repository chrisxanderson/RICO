#ifndef __RICO_NEGATE_H__
#define __RICO_NEGATE_H__

#include "RicoFunction.h"

/////////////////////////////////////////////////////////////////////////

class RicoNegate : public RicoFunction {
 public:
  RicoNegate(RicoPtr x) : RicoFunction(RicoFunction::NEG, x) {} 
  ~RicoNegate() {}

 public: // Interface
  virtual string  getExpression()                        const;
  virtual RicoPtr differentiate(RicoPtr This, RicoID Wrt) const;

}; // RicoNegate

/////////////////////////////////////////////////////////////////////////
#endif // __RICO_NEGATE_H__

