#ifndef __PRODUCT_DIGESTER_H__
#define __PRODUCT_DIGESTER_H__

#include <vector>

#include "Rico.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////
// (2a + b  + c) + (a + {-1}b)
//
//          a b c             // header of unique base components
//          3 0 1

class SumDigester {
 private:
  RicoPtr         m_constant;
  vector<RicoPtr> m_coefficients;
  vector<RicoPtr> m_bases;

 public:
  SumDigester(RicoPtr A)            : m_constant(Rico::Zero) {addTerm(A);}
  SumDigester(RicoPtr A, RicoPtr B) : m_constant(Rico::Zero) {addTerm(A); addTerm(B);}
  SumDigester(const vRicoPtr& V)    : m_constant(Rico::Zero) {addTerm(V);}
  ~SumDigester() {}

  static RicoPtr Digest(RicoPtr A, RicoPtr B) {return SumDigester(A,B).get();}
  static RicoPtr Digest(const vRicoPtr& V)    {return SumDigester(V)  .get();}

 public: // Access
  RicoPtr get() const;

 private:
  void addTerm(RicoPtr A);
  void addTerm(RicoPtr Coefficient, RicoPtr Base);
  void addTerm(const vRicoPtr& V);
};

////////////////////////////////////////////////////////////////////////////////////
#endif // __PRODUCT_DIGESTER_H__
