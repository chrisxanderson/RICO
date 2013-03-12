#ifndef __RICO_UNIFORM_H__
#define __RICO_UNIFORM_H__

#include "RicoBasicRV.h"

/////////////////////////////////////////////////////////////////////////

class RicoUniform : public RicoBasicRV {
  static const double RANGE_FACTOR = 0.05;  // Pct of range expansion for suggested range.

 private:
  double m_left;
  double m_right;

 public:
 RicoUniform() :                          RicoBasicRV(), m_left(0),    m_right(1)     {}
 RicoUniform(double left, double right) : RicoBasicRV(), m_left(left), m_right(right) {}
  ~RicoUniform() {}

 public: // Interface
  ddPair getSuggestedPlotRange() const;
  double CDF(double x) const;

 public: // C++ Only
  virtual string getDefinition() const;

}; // RicoUniform

/////////////////////////////////////////////////////////////////////////
#endif //__RICO_UNIFORM_H__
