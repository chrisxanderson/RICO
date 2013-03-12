#ifndef __RICO_CAUCHY_H__
#define __RICO_CAUCHY_H__

#include "RicoBasicRV.h"

/////////////////////////////////////////////////////////////////////////

class RicoCauchy : public RicoBasicRV {
  static const double SIGMA_FACTOR = 4.0;
  static const double PI;

 private:
  double m_mu;
  double m_sigma;

 public:
 RicoCauchy() :                        RicoBasicRV(), m_mu(0), m_sigma(1)      {}
 RicoCauchy(double mu) :               RicoBasicRV(), m_mu(mu), m_sigma(1)     {}
 RicoCauchy(double mu, double sigma) : RicoBasicRV(), m_mu(mu), m_sigma(sigma) {}
  ~RicoCauchy() {}

 public: // Interface
  ddPair getSuggestedPlotRange() const;
  double CDF(double x) const;

 public: // C++ Only
  virtual string getDefinition() const;

}; // RicoCauchy

/////////////////////////////////////////////////////////////////////////
#endif //__RICO_CAUCHY_H__
