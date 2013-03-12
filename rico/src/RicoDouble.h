#ifndef __RICO_DOUBLE_H__
#define __RICO_DOUBLE_H__

#include "RicoNumber.h"

class RicoInteger;
class RicoFraction;

/////////////////////////////////////////////////////////////////////////

class RicoDouble : public RicoNumber {
 private:
  enum DOUBLE_TYPE {REGULAR, ZERO, ONE, NEG_ONE, POS_INF, NEG_INF, NaN};

  double      m_x;
  DOUBLE_TYPE m_double_type;

 public:
  RicoDouble(); // NaN
  RicoDouble(double x);
  ~RicoDouble() {}

 public: // Interface
  double x()        const {return m_x;}
  bool   isZero()   const {return m_double_type == ZERO;}
  bool   isPos()    const {return x() > 0;}
  bool   isNeg()    const {return x() < 0;}
  bool   isOne ()   const {return m_double_type == ONE; }
  bool   isNegOne() const {return m_double_type == NEG_ONE;}
  bool   isPosInf() const {return m_double_type == POS_INF;}
  bool   isNegInf() const {return m_double_type == NEG_INF;}
  bool   isNaN()    const {return m_double_type == NaN;}

  RicoPtr getNeg    () const {return Rico::Double(-x());}
  RicoPtr getInverse() const {return Rico::Double(1/x());}
  RicoPtr Log       () const;
  RicoPtr Exp       () const;

  virtual string getExpression() const;

  // Assume all special cases are already handled
  static RicoPtr Pow(const RicoNumber  & a, const RicoNumber & b);

  bool    operator==(const RicoDouble  & b) const;
  bool    operator==(const RicoInteger & b) const;
  bool    operator==(const RicoFraction& b) const;

  RicoPtr operator+ (const RicoInteger & b) const;
  RicoPtr operator+ (const RicoFraction& b) const;
  RicoPtr operator+ (const RicoDouble  & b) const;

  RicoPtr operator* (const RicoInteger & b) const;
  RicoPtr operator* (const RicoFraction& b) const;
  RicoPtr operator* (const RicoDouble  & b) const;

}; // RicoDouble

///////////////////////////////////////////////////////////////////////////
#endif  // __RICO_DOUBLE_H__
