#include "RicoLogNormal.h"

#include <iostream>
#include <sstream>
#include <string>

#include <cmath>

// RicoLogNormal //////////////////////////////////////////////////////////////////

ddPair RicoLogNormal::getSuggestedPlotRange() const {
  return make_pair(0, m_mu+SIGMA_FACTOR*m_sigma/2); // NEEDS WORK
}

double RicoLogNormal::CDF(double x) const {
  // Pr(X <= x)
  if(x <= 0) return 0;
  return 0.5 + 0.5 * erf((log(x)-m_mu)/sqrt(2)/m_sigma);
}

string RicoLogNormal::getDefinition() const {
  stringstream sstr;
  sstr << "LogN[" << m_mu << ", " << m_sigma << "]";
  return sstr.str();
}

///////////////////////////////////////////////////////////////////////
