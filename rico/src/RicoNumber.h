#ifndef __RICO_NUMBER_H__
#define __RICO_NUMBER_H__

#include "Rico.h"

/////////////////////////////////////////////////////////////////////////

class RicoNumber : public Rico {
 public:
  enum NUMBER_TYPE {INTEGER, FRACTION, DOUBLE};

 protected:
  NUMBER_TYPE m_number_type;

 public: // Constructor
  RicoNumber(NUMBER_TYPE number_type) : m_number_type(number_type), Rico(Rico::NUMBER) {}
  ~RicoNumber() {}

  static RicoPtr Log  (RicoPtr A);
  static RicoPtr Exp  (RicoPtr A);
  static RicoPtr Power(RicoPtr A, RicoPtr B);

 public: // Access
  NUMBER_TYPE ntype() const {return m_number_type;}
  bool   isInteger () const {return m_number_type == INTEGER;}
  bool   isFraction() const {return m_number_type == FRACTION;}
  bool   isDouble  () const {return m_number_type == DOUBLE;}

  virtual bool isZero()   const = 0;
  virtual bool isPos()    const = 0;
  virtual bool isNeg()    const = 0;
  virtual bool isOne ()   const = 0;
  virtual bool isNegOne() const = 0;
  virtual bool isPosInf() const = 0;
  virtual bool isNegInf() const = 0;
  virtual bool isNaN()    const = 0;

  virtual ppPair getPair() const {
    cout << "ERROR: RicoNumber::getPair() undefined." << endl;
    return ppPair(Rico::One, Rico::One);
  }

  RicoPtr getNeg    () const;
  RicoPtr getInverse() const;
  
 public: // Operations
  virtual RicoPtr differentiate(RicoPtr This, RicoID Wrt) const {return Rico::Zero;}

  virtual bool operator==(RicoPtr That)  const;
  
  RicoPtr operator+(const RicoNumber& B) const;
  RicoPtr operator*(const RicoNumber& B) const;
  RicoPtr operator/(const RicoNumber& B) const;

}; // RicoNumber

//////////////////////////////////////////////////////////////////////////////
#endif  // __RICO_NUMBER_H__
