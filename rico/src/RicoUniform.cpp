#include "RicoUniform.h"

#include <iostream>
#include <sstream>
#include <string>

// RicoUniform //////////////////////////////////////////////////////////////////

ddPair RicoUniform::getSuggestedPlotRange() const {
  double extra = (m_right - m_left)*RANGE_FACTOR;

  return make_pair(m_left-extra, m_right+extra);
}

double RicoUniform::CDF(double x) const {
  // Pr(X <= x)
  if(x <= m_left) return 0;
  if(x >= m_right) return 1;
  return (x-m_left) / (m_right - m_left);
}

string RicoUniform::getDefinition() const {
  stringstream sstr;
  sstr << "U[" << m_left << ", " << m_right << "]";
  return sstr.str();
}

///////////////////////////////////////////////////////////////////////////////////
