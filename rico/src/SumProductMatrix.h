#ifndef __SUM_PRODUCT_MATRIX_H__
#define __SUM_PRODUCT_MATRIX_H__

#include <vector>
#include <list>

#include "Rico.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////

class SumProductMatrix {
 private:
  vector<RicoPtr> m_header;  // parallel to header
  list<iiMap>     m_matrix;  // sparse rows of powers of corresponding header expression

 public:
  SumProductMatrix();
  SumProductMatrix(RicoPtr x);
  ~SumProductMatrix() {}

 public: // Access
  void setHeader(const vector<RicoPtr>& header) {m_header = header;}
  void appendRow(const iiMap& row)              {m_matrix.push_back(row);}
  void splitRows(int idx, SumProductMatrix& A, SumProductMatrix& B) const;

  const vector<RicoPtr>& getHeader() const {return m_header;}
  const list<iiMap>&     getMatrix() const {return m_matrix;}

  int     getRows        ()        const {return getMatrix().size();}
  int     getCount       (int idx) const;
  int     getTally       (int idx) const;

  int     getMaxCount    ()        const;
  int     getMajorColumn ()        const;
  int     getMinPower    (int idx) const;

 public: // Operations
  RicoPtr getParseTree () const;

  string  getExpression()       const {return getParseTree()->getExpression();}
  string  getMatrixExpression() const;

  void    decrementColumn(int idx, int power);

  RicoPtr getFactored()                 const; // main result
  RicoPtr getRicoRow (const iiMap& row) const;

 private: // utility
  iiMap& newRow() {iiMap m; m_matrix.push_back(m); return m_matrix.back();}

  int  insertHeader  (RicoPtr x);
  void updateRow     (RicoPtr x, iiMap& row);
  void updateElement (RicoPtr x, iiMap& row);

  void  incrementElement(      iiMap& row, int idx, int power);
  static bool hasElement(const iiMap& row, int idx) {return row.find(idx) != row.end();}
  static bool hasPositiveElement(const iiMap& row, int idx) {
    iiMap::const_iterator it = row.find(idx);
    return it != row.end() && it->second > 0;
  }

};

////////////////////////////////////////////////////////////////////////////////////
#endif // __SUM_PRODUCT_MATRIX_H__
