#ifndef __RICO_FRACTION_H__
#define __RICO_FRACTION_H__

#include "RicoNumber.h"

class RicoInteger;

/////////////////////////////////////////////////////////////////////////

class RicoFraction : public RicoNumber {
  // Record in lowest form (DONE).
  // Numerator carries sign.
 private:
  int m_x;
  int m_y;

 public:
  RicoFraction(int x, int y);
  ~RicoFraction() {}

 public: // Interface
  int    x() const {return m_x;}
  int    y() const {return m_y;}
  double d() const {return ((double)x())/((double)y());}  // assumes this is valid operation

  bool isZero()   const {return x() == 0 && y() != 0;}
  bool isPos()    const {return x() > 0;}
  bool isNeg()    const {return x() < 0;}
  bool isOne ()   const {return x() ==  1 && y() == 1;}
  bool isNegOne() const {return x() == -1 && y() == 1;}
  bool isPosInf() const {return x() > 0   && y() == 0;}
  bool isNegInf() const {return x() < 0   && y() == 0;}
  bool isNaN()    const {return x() == 0  && y() == 0;}

  RicoPtr getNeg    () const {return Rico::Fraction(-x(), y());}
  RicoPtr getInverse() const {return Rico::Fraction(y(), x());}
  RicoPtr Log       () const;
  RicoPtr Exp       () const;

  virtual string getExpression() const;

  // Assume all special cases are already handled.
  // Only accepts RicoFraction and RicoInteger
  static RicoPtr Pow(const RicoNumber& a, const RicoNumber& b);

  bool operator==(const RicoFraction& b) const;
  bool operator==(const RicoInteger&  b) const;

  RicoPtr operator+(const RicoInteger & b) const;
  RicoPtr operator+(const RicoFraction& b) const;
  RicoPtr operator*(const RicoInteger & b) const;
  RicoPtr operator*(const RicoFraction& b) const;

 private: // Util
  static int     gcd(int a, int b);
  static RicoPtr upConvert(RicoPtr F);

}; // RicoFraction

////////////////////////////////////////////////////////////////////////////////////
#endif  // __RICO_FRACTION_H__
