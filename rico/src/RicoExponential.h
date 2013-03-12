#ifndef __RICO_EXPONENTIAL_H__
#define __RICO_EXPONENTIAL_H__

#include "RicoBasicRV.h"

/////////////////////////////////////////////////////////////////////////

class RicoExponential : public RicoBasicRV {
  static const double RANGE_LIMIT = 0.001; // Probability not shown.

 private:
  double m_lambda;

 public:
 RicoExponential(double lambda) : RicoBasicRV(), m_lambda(lambda) {}
  ~RicoExponential() {}

 public: // Interface
  ddPair getSuggestedPlotRange() const;
  double CDF(double x) const;

 public: // C++ Only
  virtual string getDefinition() const;

}; // RicoExponential

/////////////////////////////////////////////////////////////////////////
#endif //__RICO_EXPONENTIAL_H__

