#ifndef __RICO_WEIBULL_H__
#define __RICO_WEIBULL_H__

#include "RicoBasicRV.h"

/////////////////////////////////////////////////////////////////////////

class RicoWeibull : public RicoBasicRV {
  static const double RANGE_VALUE = 0.001;  // Amount of probability ti right of limit

 private:
  double m_shape;
  double m_scale;

 public:
 RicoWeibull(double shape) :               RicoBasicRV(), m_shape(shape), m_scale(1) {}
 RicoWeibull(double shape, double scale) : RicoBasicRV(), m_shape(shape), m_scale(scale) {}
  ~RicoWeibull() {}

 public: // Interface
  ddPair getSuggestedPlotRange() const;
  double CDF(double x) const;

 public: // C++ Only
  virtual string getDefinition() const;

}; // RicoWeibull

/////////////////////////////////////////////////////////////////////////
#endif //__RICO_WEIBULL_H__
