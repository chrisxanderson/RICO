#ifndef __RICO_CHI_SQUARED_H__
#define __RICO_CHI_SQUARED_H__

#include "RicoBasicRV.h"

/////////////////////////////////////////////////////////////////////////

class RicoChiSquared : public RicoBasicRV {
  // From: http://people.richland.edu/james/lecture/m170/tbl-chi.html
  // The x-values at the .01 probability level for the first N df starting at 1.
  // Beyond the table, we just fit a line.

  static const double RANGE_INFLATE = 0.2;
  static const double RANGE_DF[30];
  
 private:
  int m_df;

 public:
 RicoChiSquared(int df) : RicoBasicRV(), m_df(df) {}
  ~RicoChiSquared() {}

 public: // Interface
  ddPair getSuggestedPlotRange() const;
  double CDF(double x) const;

 public: // C++ Only
  virtual string getDefinition() const;

}; // RicoChiSquared

/////////////////////////////////////////////////////////////////////////
#endif //__RICO_CHI_SQUARED_H__
