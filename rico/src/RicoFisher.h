#ifndef __RICO_FISHER_H__
#define __RICO_FISHER_H__

#include "RicoBasicRV.h"

/////////////////////////////////////////////////////////////////////////

class RicoFisher : public RicoBasicRV {

 private:
  int m_df1;
  int m_df2;

 public:
 RicoFisher(int df1, int df2) : RicoBasicRV(), m_df1(df1), m_df2(df2) {}
  ~RicoFisher() {}

 public: // Interface
  ddPair getSuggestedPlotRange() const;
  double CDF(double x) const;

 public: // C++ Only
  virtual string getDefinition() const;

}; // RicoFisher

/////////////////////////////////////////////////////////////////////////
#endif //__RICO_FISHER_H__
