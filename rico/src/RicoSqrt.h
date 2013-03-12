#ifndef __RICO_SQRT_H__
#define __RICO_SQRT_H__

#include "RicoFunction.h"

/////////////////////////////////////////////////////////////////////////

class RicoSqrt : public RicoFunction {

 public:
  RicoSqrt(RicoPtr x) : RicoFunction(RicoFunction::SQRT, x) {} 
  ~RicoSqrt() {}

 public: // Interface
  virtual string  getExpression()                         const;
  virtual RicoPtr differentiate(RicoPtr This, RicoID Wrt) const;

}; // RicoSqrt

/////////////////////////////////////////////////////////////////////////
#endif // __RICO_SQRT_H__

