#ifndef __RICO_LOGISTIC_H__
#define __RICO_LOGISTIC_H__

#include "RicoBasicRV.h"

/////////////////////////////////////////////////////////////////////////

class RicoLogistic : public RicoBasicRV {
  static const double RANGE_LIMIT = 0.001; // Probability not shown.

 private:
  double m_mu;
  double m_s;

 public:
  RicoLogistic(double mu, double s) : RicoBasicRV(), m_mu(mu), m_s(s) {}
  ~RicoLogistic() {}

 public: // Interface
  ddPair getSuggestedPlotRange() const;
  double CDF(double x) const;

 public: // C++ Only
  virtual string getDefinition() const;

}; // RicoLogistic

/////////////////////////////////////////////////////////////////////////
#endif //__RICO_LOGISTIC_H__

