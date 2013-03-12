#ifndef __RICO_TRIANGULAR_H__
#define __RICO_TRIANGULAR_H__

#include "RicoBasicRV.h"

/////////////////////////////////////////////////////////////////////////

class RicoTriangular : public RicoBasicRV {
 private:
  double m_left;
  double m_middle;
  double m_right;

 public:
 RicoTriangular() : 
  RicoBasicRV(), m_left(-1.0), m_middle(0), m_right(1.0) {}
 RicoTriangular(double middle) : 
  RicoBasicRV(), m_left(middle-1), m_middle(middle), m_right(middle+1) {}
 RicoTriangular(double left, double right) : 
  RicoBasicRV(), m_left(left), m_middle((left + right) / 2), m_right(right) {}
 RicoTriangular(double left, double middle, double right) : 
  RicoBasicRV(), m_left(left), m_middle(middle), m_right(right) {}

  ~RicoTriangular() {}

 public: // Interface
  ddPair getSuggestedPlotRange() const;
  double CDF(double x) const;

 public: // C++ Only
  virtual string getDefinition() const;

}; // RicoTriangular

/////////////////////////////////////////////////////////////////////////
#endif //__RICO_TRIANGULAR_H__

