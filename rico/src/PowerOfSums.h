#ifndef __POWER_OF_SUMS_H__
#define __POWER_OF_SUMS_H__

#include "Rico.h"
#include <list>
#include <vector>

using namespace std;


// Power of Sums
/*
  The goal is to produce a sequence of all M-codes where each value represents the power of the
  summand in the product. Mcodes are powers of products all summed together.

Example:  (a+b+c)^4 = (a+b+c)*(a+b+c)*(a+b+c)*(a+b+c) = a^4 + a^3b + ... + c^4
          M = 3, N = 4
  products actual Mcode Ncode K
  aaaa     a^4    400   0000  4!/(4!*0!*0!) = 1     aaaa
  aaab     a^3b   310   0001  4!/(3!*1!*0!) = 4     aaab + aaba + abaa + baaa
  aaac	   a^3c   301   0002  4!/(3!*0!*1!) = 4     aaac + aaca + acaa + caaa
  aabb     a^2b^2 220   0011  4!/(2!*2!*0!) = 6     aabb + abab + abba + baab + baba + bbaa
  aabc 	   a^2bc  211   0012  4!/(2!*1!*1!) = 12    aabc + aacb + abac + acab + abca + acba + ..
  aacc	   a^2c^2 202   0022  4!/(2!*0!*2!) = 6
  abbb	   ab^3   130   0111  4!/(1!*3!*0!) = 4     abbb + babb + bbab + bbba
  abbc	   ab^2c  121   0112  4!/(1!*2!*1!) = 12    abbc + abcb + acbb + ...
  abcc	   abc^2  112   0122  4!/(1!*1!*2!) = 12
  accc	   ac^3   103   0222  4!/(1!*0!*3!) = 4
  bbbb     b^4    040   1111  4!/(0!*4!*0!) = 1     bbbb
  bbbc	   b^3c   031   1112  4!/(0!*3!*1!) = 4 
  bbcc	   b^2c^2 022   1122  4!/(0!*2!*2!) = 6
  bccc	   bc^3   013   1222  4!/(0!*1!*3!) = 4
  cccc	   c^4    004   2222  4!/(0!*0!*4!) = 1
  	   	                      TOTAL = 81 = 3*3*3*3
  Notice: The Ncode is easiest to produce since it's just a counting up over the possible m-tuples
          The Ncode has one caveat: the numbers are non-decreasing.
	  The Mcode is found from the Ncode by counting number of instances of each m-tuple.
	  a = 0, b = 1, c = 2, etc.
	  The Mcode is easiest to use for the expression and the associated coefficient.

  case (2,3)
  Ncode  Mcode
  000    30
  001    21
  011    12
  111    03

  case (3,1)
  Ncode  Mcode
  0      100
  1      010
  2      001

  case (1,3)
  Ncode Mcode
  000   3
*/

class PowerOfSums {
  vector<int> m_v;
  int m_M;  // number discrete values possible
  int m_N;  // number of M-tuples "digits"

public:
  PowerOfSums(int M, int N) : m_M(M), m_N(N) { 
    for(int i = 0; i < N; i++) m_v.push_back(0); 
  }

  static void    getCode(list<int>& Mcode, int M, int N);
  static RicoPtr Power(const vRicoPtr& P, int N);

private:
  static RicoPtr getTerm(const vRicoPtr& P, 
			 const list<int>& Mcode, 
			 list<int>::const_iterator& it, 
			 int M, int N);

  void convertCodeType(const list<int>& Ncode, list<int>& Mcode);
  bool getNextCode(list<int>& code);
  bool increment();
  static long factorial(int N);
};

////////////////////////////////////////////////////////////////////////////////////
#endif // __POWER_OF_SUMS_H__
