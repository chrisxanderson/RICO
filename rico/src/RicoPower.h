#ifndef __RICO_POWER_H__
#define __RICO_POWER_H__

#include "RicoFunction.h"

/////////////////////////////////////////////////////////////////////////

class RicoPower : public RicoFunction {
 public:
  RicoPower(RicoPtr x, RicoPtr y) : RicoFunction(RicoFunction::POW, x, y) {}
  RicoPower(const vRicoPtr& v)    : RicoFunction(RicoFunction::POW, v)    {}
  ~RicoPower() {}

 public: // Interface
  virtual string  getExpression()                         const;
  virtual RicoPtr differentiate(RicoPtr This, RicoID Wrt) const;

}; // RicoPower

/////////////////////////////////////////////////////////////////////////
#endif //__RICO_POWER_H__

