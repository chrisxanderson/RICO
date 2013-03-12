#ifndef __RICO_INTEGER_H__
#define __RICO_INTEGER_H__

#include "RicoNumber.h"

/////////////////////////////////////////////////////////////////////////

class RicoInteger : public RicoNumber {

 private:
  int m_x;

 public:
  RicoInteger(int x) : RicoNumber(RicoNumber::INTEGER), m_x(x) {}
  ~RicoInteger() {}

 public: // Interface
  int x() const {return m_x;}

  virtual string getExpression() const;

  virtual bool isZero()   const {return x() == 0;}
  virtual bool isPos()    const {return x() > 0;}
  virtual bool isNeg()    const {return x() < 0;}
  virtual bool isOne ()   const {return x() == 1;}
  virtual bool isNegOne() const {return x() == -1;}
  virtual bool isPosInf() const {return false;}
  virtual bool isNegInf() const {return false;}
  virtual bool isNaN()    const {return false;}

  RicoPtr getNeg    () const {return Rico::Integer(-x());}
  RicoPtr getInverse() const {return Rico::Fraction(1, x());}
  RicoPtr Log       () const;
  RicoPtr Exp       () const;

  // Assume all special cases are already handled
  static RicoPtr Pow(const RicoNumber& a, const RicoNumber& b);

  bool    operator==(const RicoInteger& b) const {return x() == b.x();}
  RicoPtr operator+ (const RicoInteger& b) const {return Rico::Integer(x() + b.x());}
  RicoPtr operator* (const RicoInteger& b) const {return Rico::Integer(x() * b.x());}

}; // RicoInteger

#endif  // __RICO_INTEGER_H__
