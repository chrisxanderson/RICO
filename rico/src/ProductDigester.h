#ifndef __PRODUCT_PRODUCT_MATRIX_H__
#define __PRODUCT_PRODUCT_MATRIX_H__

#include <vector>
#include <list>
#include <map>

#include "Rico.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////

//                  a  b c             // header of unique base components
//      a^2*b*c  -> 2  1 1             // We can assume each component us unique (simplest)
//      a*b^{-1} -> 1 -1 o             // some component may overlap with first row.

class ProductDigester {
 private:
  RicoPtr          m_coefficient;      // must be created cumulatively
  vector< RicoPtr> m_bases;
  vector<vRicoPtr> m_exponents;

 public:
  ProductDigester(RicoPtr A)            : m_coefficient(Rico::One) {addTerm(A);}
  ProductDigester(RicoPtr A, RicoPtr B) : m_coefficient(Rico::One) {addTerm(A); addInvTerm(B);}
  ProductDigester(const vRicoPtr& V)    : m_coefficient(Rico::One) {addTerm(V);}
  ~ProductDigester() {}

  static RicoPtr Digest(RicoPtr A)            {return ProductDigester(A)  .get();}
  static RicoPtr Digest(RicoPtr A, RicoPtr B) {return ProductDigester(A,B).get();}
  static RicoPtr Digest(const vRicoPtr& V)    {return ProductDigester(V)  .get();}

 public: // Access
  RicoPtr get() const;

 private:
  void addTerm   (RicoPtr A);
  void addInvTerm(RicoPtr A);
  void addTerm   (const vRicoPtr& V);
  void addInvTerm(const vRicoPtr& V);
  void addTerm   (RicoPtr Base, RicoPtr Exponent);

  static bool isNegative(RicoPtr A);
  static RicoPtr Negate (RicoPtr A);
};

////////////////////////////////////////////////////////////////////////////////////
#endif // __PRODUCT_PRODUCT_MATRIX_H__
