#ifndef __RICO_SIMPLIFY_H__
#define __RICO_SIMPLIFY_H__

#include "Rico.h"
#include <list>

using namespace std;

////////////////////////////////////////////////////////////////////////////////////

typedef list<ppPair> lppPair;

class RicoSimplify {
 private: 
 RicoSimplify() {} // Don't instantiate

 public: // Stage One simplificiation
 static RicoPtr Simplify(RicoPtr A);

 public: // Internal and Testing
 static RicoPtr Negate  (RicoPtr A);
 static RicoPtr Add     (RicoPtr A, RicoPtr B);
 static RicoPtr Add     (const vRicoPtr& P);
 static RicoPtr Subtract(RicoPtr A, RicoPtr B);
 static RicoPtr Multiply(RicoPtr A, RicoPtr B);
 static RicoPtr Multiply(RicoPtr A, RicoPtr B, RicoPtr C);
 static RicoPtr Multiply(const vRicoPtr& P);
 static RicoPtr Divide  (RicoPtr A, RicoPtr B);
 static RicoPtr Power   (RicoPtr A, RicoPtr B);
 static RicoPtr Log     (RicoPtr A);
 static RicoPtr Exp     (RicoPtr A);
 static RicoPtr Sqrt    (RicoPtr A);
 static RicoPtr Factor  (RicoPtr A);

 public: // Util
 static void    PowerOfSums(list<int>& code, int M, int N, int m = 0);
};

////////////////////////////////////////////////////////////////////////////////////
#endif // __RICO_SIMPLIFY_H__
