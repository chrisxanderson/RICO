#ifndef __RICO_BETA_H__
#define __RICO_BETA_H__

#include "RicoBasicRV.h"

/////////////////////////////////////////////////////////////////////////

class RicoBeta : public RicoBasicRV {

 private:
  double m_alpha;
  double m_beta;

 public:
 RicoBeta(double alpha, double beta) : RicoBasicRV(), m_alpha(alpha), m_beta(beta) {}
  ~RicoBeta() {}


 public: // Interface
  ddPair getSuggestedPlotRange() const {return make_pair(0,1);}
  double CDF(double x) const;

 public: // C++ Only
  virtual string getDefinition() const;

}; // RicoBeta

/////////////////////////////////////////////////////////////////////////
#endif //__RICO_BETA_H__

