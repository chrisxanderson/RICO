#ifndef __RICO_DIVIDE_H__
#define __RICO_DIVIDE_H__

#include "RicoFunction.h"

/////////////////////////////////////////////////////////////////////////

class RicoDivide : public RicoFunction {
 public:
  RicoDivide(RicoPtr x, RicoPtr y) : RicoFunction(RicoFunction::DIV, x, y) {}
  ~RicoDivide();

 public: // Interface
  virtual string  getExpression()                        const;
  virtual RicoPtr differentiate(RicoPtr This, RicoID Wrt) const;

}; // RicoDivide

/////////////////////////////////////////////////////////////////////////
#endif //__RICO_DIVIDE_H__

