#ifndef __RICO_LOG_H__
#define __RICO_LOG_H__

#include "RicoFunction.h"

/////////////////////////////////////////////////////////////////////////

class RicoLog : public RicoFunction {

 public:
  RicoLog(RicoPtr x) : RicoFunction(RicoFunction::LOG, x) {} 
  ~RicoLog() {}

 public: // Interface
  virtual string  getExpression()                         const;
  virtual RicoPtr differentiate(RicoPtr This, RicoID Wrt) const;

}; // RicoLog

/////////////////////////////////////////////////////////////////////////
#endif // __RICO_LOG_H__

