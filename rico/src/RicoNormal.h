#ifndef __RICO_NORMAL_H__
#define __RICO_NORMAL_H__

#include "RicoBasicRV.h"

/////////////////////////////////////////////////////////////////////////

class RicoNormal : public RicoBasicRV {
  static const double SIGMA_FACTOR = 4.0;

 private:
  double m_mu;
  double m_sigma;

 public:
 RicoNormal() :                        RicoBasicRV(), m_mu(0), m_sigma(1)      {}
 RicoNormal(double mu) :               RicoBasicRV(), m_mu(mu), m_sigma(1)     {}
 RicoNormal(double mu, double sigma) : RicoBasicRV(), m_mu(mu), m_sigma(sigma) {}
  ~RicoNormal() {}

 public: // Access
  virtual ddPair getSuggestedPlotRange() const;
  virtual double CDF(double x) const;
  virtual string getDefinition() const;

}; // RicoNormal

/////////////////////////////////////////////////////////////////////////
#endif //__RICO_NORMAL_H__

