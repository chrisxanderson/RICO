#include "RicoNormal.h"

#include <iostream>
#include <sstream>
#include <string>

#include <cmath>

// RicoNormal //////////////////////////////////////////////////////////////////

ddPair RicoNormal::getSuggestedPlotRange() const {
  return make_pair(m_mu-SIGMA_FACTOR*m_sigma, m_mu+SIGMA_FACTOR*m_sigma);
}

double RicoNormal::CDF(double x) const {
  // Pr(X <= x)
  return 0.5 * erfc((m_mu - x)/sqrt(2)/m_sigma);
}

string RicoNormal::getDefinition() const {
  stringstream sstr;
  sstr << "Norm[u=" << m_mu << ", s=" << m_sigma << "]";
  return sstr.str();
}

//////////////////////////////////////////////////////////////////////////////////
