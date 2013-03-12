#ifndef __RICO_LOGNORMAL_H__
#define __RICO_LOGNORMAL_H__

#include "RicoBasicRV.h"

/////////////////////////////////////////////////////////////////////////

class RicoLogNormal : public RicoBasicRV {
  static const double SIGMA_FACTOR = 4.0;

 private:
  double m_mu;
  double m_sigma;

 public:
 RicoLogNormal() :                        RicoBasicRV(), m_mu(0), m_sigma(1)  {}
 RicoLogNormal(double mu) :               RicoBasicRV(), m_mu(mu), m_sigma(1) {}
 RicoLogNormal(double mu, double sigma) : RicoBasicRV(), m_mu(mu), m_sigma(sigma) {}
  ~RicoLogNormal() {}

 public: // Interface
  ddPair getSuggestedPlotRange() const;
  double CDF(double x) const;

 public: // C++ Only
  virtual string getDefinition() const;

}; // RicoLogNormal

/////////////////////////////////////////////////////////////////////////
#endif //__RICO_LOGNORMAL_H__
