#include <iostream>
#include <sstream>
#include <string>

#include "RicoTriangular.h"

// RicoTriangular //////////////////////////////////////////////////////////////////

ddPair RicoTriangular::getSuggestedPlotRange() const {
  double span = m_right-m_left;
  return make_pair(m_left-0.05*span, m_right+0.05*span);
}

double RicoTriangular::CDF(double x) const {
  // See notes: 9/27/12 p2. Just integrate the triangle.
  // Pr(X <= x)
  double L = m_left; double M = m_middle; double R = m_right;
  double cdf;
  if(x <= L)      cdf = 0;
  else if(x >= R) cdf = 1;
  else if(x < M)  cdf = (x-L)*(x-L)/(M-L)/(R-L);
  else            cdf = (M-L)/(R-L) + (M-x)*(x+M-2*R)/(R-L)/(R-M);
  return cdf;
}

string RicoTriangular::getDefinition() const {
  stringstream sstr;
  sstr << "Tri[" << m_left << ", " << m_middle << ", " << m_right << "]";
  return sstr.str();
}

///////////////////////////////////////////////////////////////////////////////////////
