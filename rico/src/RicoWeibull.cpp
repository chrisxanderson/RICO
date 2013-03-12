#include "RicoWeibull.h"

#include <iostream>
#include <sstream>
#include <string>

#include <cmath>

// RicoWeibull //////////////////////////////////////////////////////////////////

ddPair RicoWeibull::getSuggestedPlotRange() const {
  double R = m_scale*pow(-log(RANGE_VALUE), 1/m_shape);
  return make_pair(0, R);
}

double RicoWeibull::CDF(double x) const {
  // Pr(X <= x)
  if(x <= 0) return 0;
  return 1 - exp( - pow(x/m_scale, m_shape));
}

string RicoWeibull::getDefinition() const {
  stringstream sstr;
  sstr << "Weibull[" << m_shape << ", " << m_scale << "]";
  return sstr.str();
}

///////////////////////////////////////////////////////////////////////////////////
