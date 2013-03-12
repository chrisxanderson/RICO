#ifndef __RICO_ADD_H__
#define __RICO_ADD_H__

#include "RicoFunction.h"

/////////////////////////////////////////////////////////////////////////

class RicoAdd : public RicoFunction {
 public:
  RicoAdd(RicoPtr x, RicoPtr y)            : RicoFunction(RicoFunction::ADD, x, y)    {} 
  RicoAdd(RicoPtr x, RicoPtr y, RicoPtr z) : RicoFunction(RicoFunction::ADD, x, y, z) {} 
  RicoAdd(const vRicoPtr& v)               : RicoFunction(RicoFunction::ADD, v)       {}
  ~RicoAdd() {}

 public: // Interface
  virtual string  getExpression()                         const;
  virtual RicoPtr differentiate(RicoPtr This, RicoID Wrt) const;

 public: // Operators
  virtual bool operator==(RicoPtr That) const;

}; // RicoAdd

/////////////////////////////////////////////////////////////////////////
#endif // __RICO_ADD_H__

