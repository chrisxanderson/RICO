#include <iostream>
#include <string>
#include <sstream>

#include "SumProductMatrix.h"

#include "RicoAdd.h"
#include "RicoMultiply.h"
#include "RicoPower.h"
#include "RicoNumber.h"
#include "RicoInteger.h"
#include "RicoOperators.h"

#include "UnitTest++/UnitTest++.h"

// SumProductMatrix /////////////////////////////////////////////////////////////////////

SumProductMatrix::SumProductMatrix() {}

SumProductMatrix::SumProductMatrix(RicoPtr x) {
  // Decode AB+BC+C
  if(x->isAdd()) {
    const RicoAdd& a = (const RicoAdd&) *x;
    for(vRicoPtr::const_iterator it = a.get_params().begin(); it != a.get_params().end(); it++) 
      updateRow(*it, newRow());
  } else 
    updateRow(x, newRow());
}


void SumProductMatrix::updateRow(RicoPtr x, iiMap& row) {
  // AB^2C
  if(x->isMultiply()) {
    const RicoMultiply& P = (const RicoMultiply&) *x;
    for(vRicoPtr::const_iterator it = P.get_params().begin(); it != P.get_params().end(); it++) {
      RicoPtr p = *it;
      updateElement(p, row);
    }
  } else 
    updateElement(x, row);
}

void SumProductMatrix::updateElement(RicoPtr x, iiMap& row) {
  if(x->isPower()) {
    const RicoPower& P = (const RicoPower&) *x;
    if(P.y()->isNumber()) {
      const RicoNumber& N = (const RicoNumber&) *(P.y());
      if(N.isInteger()) {
	const RicoInteger& I = (const RicoInteger&) N;
	int n = I.x();
	int i = insertHeader(P.x());
	incrementElement(row, i, n);
      }
    }
  } else {
    int i = insertHeader(x);
    incrementElement(row, i, 1);
  }
}

void SumProductMatrix::decrementColumn(int idx, int power) {
  // Assume every row has the given (idx) column.
  // Could create an empty row. We do not erase it.
  for(list<iiMap>::iterator it = m_matrix.begin(); it != m_matrix.end(); it++) {
    iiMap& row = *it;
    iiMap::iterator jt = row.find(idx);
    if(jt == row.end()) {
      cout << "ERROR: SumProductMatrix::decrementColumn invalid column: " << idx << endl;
      continue;
    }
    pair<const int, int>& kv = (*jt);
    kv.second -= power;
  }
}

void SumProductMatrix::incrementElement(iiMap& row, int idx, int power) {
  if(hasElement(row, idx)) row[idx] += power;
  else                     row[idx]  = power;
}

int SumProductMatrix::insertHeader(RicoPtr x) {
  int i = 0;
  for(vector<RicoPtr>::const_iterator it = getHeader().begin(); it != getHeader().end(); it++,i++) {
    RicoPtr s = *it; 
    if(s == x) return i;
  }
  m_header.push_back(x);
  return i;
}

RicoPtr SumProductMatrix::getParseTree() const {
  // AB^2C + CD + E  --- Product == row, sum == rows
  vRicoPtr S;
  for(list<iiMap>::const_iterator it = getMatrix().begin(); it != getMatrix().end(); it++) {
    vRicoPtr P;
    const iiMap& row = *it; 
    for(iiMap::const_iterator it = row.begin(); it != row.end(); it++) {
      RicoPtr expr  = getHeader()[it->first];
      RicoPtr power = Rico::Integer(it->second);
      P.push_back(Rico::Power(expr, power));
    }
    S.push_back(Rico::Multiply(P));
  }
  return Rico::Add(S);
}

string SumProductMatrix::getMatrixExpression() const {
  stringstream ss;
  ss << "{";
  for(vector<RicoPtr>::const_iterator it = getHeader().begin(); it != getHeader().end(); it++) {
    ss << *it;
    if(it + 1 != getHeader().end()) ss << " ";
  }
  ss << "}" << endl;
  for(list<iiMap>::const_iterator it = getMatrix().begin(); it != getMatrix().end(); it++) {
    ss << "|";
    for(int idx = 0; idx < getHeader().size(); idx++) {
      if(hasElement(*it, idx)) ss << it->find(idx)->second;
      else                     ss << "-";
      if(idx+1 != getHeader().size()) ss << " ";
    }
    ss << "|" << endl;
  }
  ss << "<";
  for(int i = 0; i < getHeader().size(); i++) {
    ss << getCount(i);
    if(i+1 < getHeader().size()) ss << " ";
  }
  ss << "> count" << endl;
  ss << "<";
  for(int i = 0; i < getHeader().size(); i++) {
    ss << getTally(i);
    if(i+1 < getHeader().size()) ss << " ";
  }
  ss << "> tally" << endl;
  return ss.str();
}

int SumProductMatrix::getMaxCount() const {
  int max_count = 0;
  for(int i = 0; i < getHeader().size(); i++) {
    int count = getCount(i);
    if(count > max_count) max_count = count;
  }
  return max_count;
}

int SumProductMatrix::getMajorColumn() const {
  // Find the maximum count
  int max_count = getMaxCount();

  if(max_count == 0) {
    cout << "ERROR: SumProductMatrix::getMajorColumn(), max_count == 0" << endl;
    return -1;
  }
  
  if(max_count == 1) {
    // Find first non-zero-count index
    for(int i = 0; i < getHeader().size(); i++) if(getCount(i) == max_count) return i;
  }

  // Gather indices that share the max count
  list<int> idx_list;
  for(int i = 0; i < getHeader().size(); i++) if(getCount(i) == max_count) idx_list.push_back(i);

  // Among the max count indices, find the max tally
  int max_tally = 0;
  for(list<int>::const_iterator it = idx_list.begin(); it != idx_list.end(); it++) {
    int tally = getTally(*it);
    if(tally > max_tally) max_tally = tally;
  }

  // Return the first max_count index that matches the max tally
  for(list<int>::const_iterator it = idx_list.begin(); it != idx_list.end(); it++)
    if(getTally(*it) == max_tally) return *it;
}

int SumProductMatrix::getMinPower(int idx) const {
  int power = 0;
  for(list<iiMap>::const_iterator it = getMatrix().begin(); it != getMatrix().end(); it++) {
    const iiMap& row = *it;
    iiMap::const_iterator jt = row.find(idx);
    if(jt == row.end()) continue;
    if(power == 0 || jt->second < power) power = jt->second;
  }
  return power;
}

int SumProductMatrix::getCount(int idx) const {
  int count = 0;
  for(list<iiMap>::const_iterator it = getMatrix().begin(); it != getMatrix().end(); it++) 
    if(hasPositiveElement(*it, idx)) count++;
  return count;
}

int SumProductMatrix::getTally(int idx) const {
  int tally = 0;
  for(list<iiMap>::const_iterator it = getMatrix().begin(); it != getMatrix().end(); it++) {
    const iiMap& row = *it;
    iiMap::const_iterator jt = row.find(idx);
    if(jt == row.end()) continue;
    tally += jt->second;
  }
  return tally;
}

void SumProductMatrix::splitRows(int idx, SumProductMatrix& A, SumProductMatrix& B) const {
  A.setHeader(getHeader());
  B.setHeader(getHeader());
  for(list<iiMap>::const_iterator it = getMatrix().begin(); it != getMatrix().end(); it++) {
    if(hasPositiveElement(*it, idx)) A.appendRow(*it);
    else                             B.appendRow(*it);
  }
}

RicoPtr SumProductMatrix::getFactored() const {
  if(getMaxCount() <= 1) {
    vRicoPtr S;
    for(list<iiMap>::const_iterator it = getMatrix().begin(); it != getMatrix().end(); it++)
      S.push_back(getRicoRow(*it));
    return Rico::Add(S);
  } else {
    // The whole factoring process boils down to this sequence
    int idx   = getMajorColumn();
    int power = getMinPower(idx);
    RicoPtr F = Rico::Power(getHeader()[idx], Rico::Integer(power));

    SumProductMatrix SPMA;
    SumProductMatrix SPMB;
    splitRows(idx, SPMA, SPMB);
    SPMA.decrementColumn(idx, power);
    if(SPMB.getRows() == 0) {
      RicoPtr A = SPMA.getFactored();
      return Rico::Multiply(F, A);
    } else {
      RicoPtr A = SPMA.getFactored();
      RicoPtr B = SPMB.getFactored();
      return Rico::Add(Rico::Multiply(F, A), B);
    }
  }
}

RicoPtr SumProductMatrix::getRicoRow(const iiMap& row) const {
  if(row.size() == 0) return Rico::One;
  vRicoPtr M;
  for(iiMap::const_iterator it = row.begin(); it != row.end(); it++) {
    const pair<const int, int>& ii = *it;
    RicoPtr E = getHeader()[ii.first];
    int power = ii.second;
    RicoPtr EP = Rico::Power(E, Rico::Integer(power));
    M.push_back(EP);
  }
  return Rico::Multiply(M);
}

// Testing //////////////////////////////////////////////////////////////////////////

TEST(TestSumProductMatrix) {
  // Given (3 * X^2 * Y) + (Y * Z) + Z
  RicoPtr X  = Rico::Normal();
  RicoPtr Y  = Rico::Triangular();
  RicoPtr Z  = Rico::Cauchy();
  RicoPtr E1 = Rico::Multiply(Rico::Three, Rico::Square(X), Y);
  RicoPtr E  = Rico::Add(E1, Rico::Multiply(Y,Z), Z);

  // 3X^2Y + YZ + Z
  RicoPtr X2Y = Rico::Multiply(Rico::Square(X),Y);
  RicoPtr YZ  = Rico::Multiply(Y,Z);
  CHECK_EQUAL(Rico::Add(Rico::Multiply(Rico::Three, X2Y), YZ, Z), E);

  SumProductMatrix SPM = SumProductMatrix(E);
  // 3X^2Y + YZ + Z == 3X^2Y + YZ + Z
  CHECK_EQUAL(E, SPM.getParseTree());

  stringstream ssM;
  ssM << "{3 X" << X->id() << " X" << Y->id() << " X" << Z->id() << "}" << endl;
  ssM << "|1 2 1 -|" << endl;
  ssM << "|- - 1 1|" << endl;
  ssM << "|- - - 1|" << endl;
  ssM << "<1 1 2 2> count" << endl;
  ssM << "<1 2 2 2> tally" << endl;
  CHECK_EQUAL(ssM.str(), SPM.getMatrixExpression());

  int idx       = SPM.getMajorColumn();
  int power     = SPM.getMinPower(idx);
  int count     = SPM.getCount(idx);
  int max_count = SPM.getMaxCount();
  CHECK_EQUAL(2, idx);
  CHECK_EQUAL(1, power);
  CHECK_EQUAL(2, count);
  CHECK_EQUAL(2, max_count);
  
  SumProductMatrix SPMA;
  SumProductMatrix SPMB;
  SPM.splitRows(idx, SPMA, SPMB);

  stringstream ssMA;
  ssMA << "{3 X" << X->id() << " X" << Y->id() << " X" << Z->id() << "}" << endl;
  ssMA << "|1 2 1 -|" << endl;
  ssMA << "|- - 1 1|" << endl;
  ssMA << "<1 1 2 1> count" << endl;
  ssMA << "<1 2 2 1> tally" << endl;
  CHECK_EQUAL(ssMA.str(), SPMA.getMatrixExpression());

  SPMA.decrementColumn(idx, power);
  stringstream ssMA2;
  ssMA2 << "{3 X" << X->id() << " X" << Y->id() << " X" << Z->id() << "}" << endl;
  ssMA2 << "|1 2 0 -|" << endl;
  ssMA2 << "|- - 0 1|" << endl;
  ssMA2 << "<1 1 0 1> count" << endl;
  ssMA2 << "<1 2 0 1> tally" << endl;
  CHECK_EQUAL(ssMA2.str(), SPMA.getMatrixExpression());
  int max_countA = SPMB.getMaxCount();
  CHECK_EQUAL(1, max_countA);

  stringstream ssMB;
  ssMB << "{3 X" << X->id() << " X" << Y->id() << " X" << Z->id() << "}" << endl;
  ssMB << "|- - - 1|" << endl;
  ssMB << "<0 0 0 1> count" << endl;
  ssMB << "<0 0 0 1> tally" << endl;
  CHECK_EQUAL(ssMB.str(), SPMB.getMatrixExpression());

  int idxB       = SPMB.getMajorColumn();
  int powerB     = SPMB.getMinPower(idxB);
  int countB     = SPMB.getCount(idxB);
  int max_countB = SPMB.getMaxCount();
  CHECK_EQUAL(3, idxB);
  CHECK_EQUAL(1, powerB);
  CHECK_EQUAL(1, countB);
  CHECK_EQUAL(1, max_countB);

  // Y*(3X^2 + Z) + Z
  RicoPtr T=Rico::Add(Rico::Multiply(Y,Rico::Add(Rico::Multiply(Rico::Three,Rico::Square(X)),Z)),Z);
  CHECK_EQUAL(T, SPM.getFactored());
}

TEST(TestSumProductMatrix2) {
  // X*5 + X*Y*5^2 + X*Y*Z*5^3  =>  X(5 + Y(25 + 125Z)
  RicoPtr X          = Rico::Normal();
  RicoPtr Y          = Rico::Triangular();
  RicoPtr Z          = Rico::Cauchy();
  RicoPtr One        = Rico::One;
  RicoPtr Two        = Rico::Two;
  RicoPtr Three      = Rico::Three;
  RicoPtr Five       = Rico::Integer(5);
  RicoPtr E1         = Rico::Multiply(X, Five);
  RicoPtr E2         = Rico::Multiply(X, Y, Rico::Power(Five, Two));
  vRicoPtr PE3;
  PE3.push_back(X); PE3.push_back(Y); PE3.push_back(Z) ;PE3.push_back(Rico::Power(Five, Three)); 
  RicoPtr E3    = Rico::Multiply(PE3);
  RicoPtr E     = Rico::Add(E1, E2, E3);

  SumProductMatrix SPM = SumProductMatrix(E);
  // 5X + 5^2XY + 5^3XYZ 
  RicoPtr T1 = Rico::Add(Rico::Multiply(Five, X), 
			 Rico::Multiply(Rico::Square(Five), X, Y), 
			 Rico::Multiply(Rico::Power(Five, Three), Rico::Multiply(X,Y,Z)));
  CHECK_EQUAL(T1, SPM.getParseTree());
  
  stringstream ssM;
  ssM << "{X" << X->id() << " 5 X" << Y->id() << " X" << Z->id() << "}" << endl;
  ssM << "|1 1 - -|" << endl;
  ssM << "|1 2 1 -|" << endl;
  ssM << "|1 3 1 1|" << endl;
  ssM << "<3 3 2 1> count" << endl;
  ssM << "<3 6 2 1> tally" << endl;
  CHECK_EQUAL(ssM.str(), SPM.getMatrixExpression());

  // 5X(1 + 5Y(1 + 5Z)
  RicoPtr T2 = Rico::Multiply(Five, X, Rico::Add(One, Rico::Multiply(Five, Y, Rico::Add(One, Rico::Multiply(Five, Z)))));
  CHECK_EQUAL(T2, SPM.getFactored());
}

/////////////////////////////////////////////////////////////////////////////////////////
