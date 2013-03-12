#include <iostream>
#include <string>
#include <sstream>

#include "PowerOfSums.h"
#include "RicoOperators.h"
#include "RicoSimplify.h"

#include "UnitTest++/UnitTest++.h"

///////////////////////////////////////////////////////////

void PowerOfSums::getCode(list<int>& Mcode, int M, int N) { 
  PowerOfSums pos = PowerOfSums(M, N);
  list<int> Ncode; // this code is modulo N
  while(pos.getNextCode(Ncode)) {} 
  pos.convertCodeType(Ncode, Mcode);
}

RicoPtr PowerOfSums::Power(const vRicoPtr& P, int N) {
  int M = P.size();
  list<int> Mcode;
  PowerOfSums::getCode(Mcode, M, N);
  vRicoPtr Q;
  list<int>::const_iterator it = Mcode.begin();
  while(it != Mcode.end()) Q.push_back(getTerm(P, Mcode, it, M, N));
  if(Q.size() == 1) return Q.front();
  return RicoSimplify::Add(Q);
}

RicoPtr PowerOfSums::getTerm(const vRicoPtr& P, 
			     const list<int>& Mcode, 
			     list<int>::const_iterator& it, 
			     int M, int N) {
  vRicoPtr Q;
  list<int>::const_iterator jt = it;
  vRicoPtr::const_iterator  pt = P.begin();
  long NF = factorial(N);
  for(int i = 0; i < M; i++, jt++) NF /= factorial(*jt);
  if(NF > 1) Q.push_back(Rico::Integer(NF));
  for(int i = 0; i < M; i++, it++, pt++) {
    if(*it == 0) continue;
    Q.push_back(RicoSimplify::Power(*pt, Rico::Integer(*it)));
  }
  if(Q.size() == 1) return Q.front();
  // If P has a number it will be first in list
  if(NF > 1 && Mcode.front() > 0 && P.front()->isNumber()) {
    RicoPtr N1   = Q.front(); Q.pop_front();
    RicoPtr N2   = Q.front(); Q.pop_front();
    RicoPtr N1N2 = RicoSimplify::Multiply(N1, N2);
    Q.push_front(N1N2);
  }
  return RicoSimplify::Multiply(Q);
}

void PowerOfSums::convertCodeType(const list<int>& Ncode, list<int>& Mcode) {
  // Ex[3,4] convert 0011 to 220 since 2 zeros, 2 ones and 0 two's
  list<int>::const_iterator it = Ncode.begin();
  while(it != Ncode.end()) {
    vector<int> mcode; 
    for(int j = 0; j < m_M; j++) mcode.push_back(0);
    for(int i = 0; i < m_N; i++, it++) mcode[(*it)]++;
    for(int j = 0; j < m_M; j++) Mcode.push_back(mcode[j]);
  }
}

bool PowerOfSums::getNextCode(list<int>& code) {
  for(int i = 0; i < m_N; i++) code.push_back(m_v[i]);
  return increment();
}

bool PowerOfSums::increment() {
  if(m_v[m_N-1] < m_M - 1) {m_v[m_N - 1]++; return true;}
  for(int i = m_N-2; i >= 0; i--)  // find the next digit to increment
    if(m_v[i] < m_M - 1) {
      int d = m_v[i] + 1;
      for(int j = i; j < m_N; j++) m_v[j] = d; // reset to new "ground state" (non-decreasing)
      return true;
    }
  return false;
}

long PowerOfSums::factorial(int N) {
  if(N == 0) return 1;
  long F = 1;
  for(int i = 1; i <= N; i++) F *= i;
  return F;
}

// Testing /////////////////////////////////////////////////////////////////////////////////

TEST(TestPowerOfSums1) {
  // Remember that the codes all run together. They are module-M (called M-codes above)
  // where M is the number of summands and N is the exponent power.

  list<int> code22;
  PowerOfSums::getCode(code22, 2, 2);
  stringstream ss22;
  for(list<int>::const_iterator it = code22.begin(); it != code22.end(); it++) ss22 << *it;
  CHECK_EQUAL("201102", ss22.str());

  list<int> code23;
  PowerOfSums::getCode(code23, 2, 3);
  stringstream ss23;
  for(list<int>::const_iterator it = code23.begin(); it != code23.end(); it++) ss23 << *it;
  CHECK_EQUAL("30211203", ss23.str());

  list<int> code31;
  PowerOfSums::getCode(code31, 3, 1);
  stringstream ss31;
  for(list<int>::const_iterator it = code31.begin(); it != code31.end(); it++) ss31 << *it;
  CHECK_EQUAL("100010001", ss31.str());

  list<int> code13;
  PowerOfSums::getCode(code13, 1, 3);
  stringstream ss13;
  for(list<int>::const_iterator it = code13.begin(); it != code13.end(); it++) ss13 << *it;
  CHECK_EQUAL("3", ss13.str());

  list<int> code34;
  PowerOfSums::getCode(code34, 3, 4);
  stringstream ss34;
  for(list<int>::const_iterator it = code34.begin(); it != code34.end(); it++) ss34 << *it;
  CHECK_EQUAL("400310301220211202130121112103040031022013004", ss34.str());
}

TEST(TestPowerOfSums2) {
  RicoPtr X = Rico::Normal();
  RicoPtr Y = Rico::Cauchy();

  vRicoPtr XY; XY.push_back(X); XY.push_back(Y);
  RicoPtr Q1 = PowerOfSums::Power(XY,2);
  CHECK_EQUAL(Rico::Add(Rico::Square(X), Rico::Multiply(Rico::Two, X, Y), Rico::Square(Y)), Q1); 

  vRicoPtr PX; PX.push_back(X);
  RicoPtr Q2 = PowerOfSums::Power(PX,1);
  CHECK_EQUAL(X, Q2);

  vRicoPtr PX2; PX2.push_back(X);
  RicoPtr Q3 = PowerOfSums::Power(PX2,2);
  CHECK_EQUAL(Rico::Square(X), Q3);

  vRicoPtr P3X; P3X.push_back(Rico::Three); P3X.push_back(X);
  RicoPtr Q4 = PowerOfSums::Power(P3X,2);
  CHECK_EQUAL(Rico::Add(Rico::Integer(9), Rico::Multiply(Rico::Integer(6),X), Rico::Square(X)), Q4);
}

//////////////////////////////////////////////////////////////////////////////////////////
