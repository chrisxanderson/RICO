#include "RicoLogistic.h"

#include <iostream>
#include <sstream>
#include <string>

#include <cmath>

// RicoLogistic //////////////////////////////////////////////////////////////////////

ddPair RicoLogistic::getSuggestedPlotRange() const {
  // solve CDF = RANGE_LIMIT and CDF = 1-RANGE_LIMIT
  double left  = m_mu - m_s * log(1/   RANGE_LIMIT  - 1);
  double right = m_mu - m_s * log(1/(1-RANGE_LIMIT) - 1);
  return make_pair(left, right);
}

double RicoLogistic::CDF(double x) const {
  // Pr(X <= x)
  return 1/(1+exp(-(x-m_mu)/m_s));
}

string RicoLogistic::getDefinition() const {
  stringstream sstr;
  sstr << "Logistic[" << m_mu << ", " << m_s << "]";
  return sstr.str();
}

/////////////////////////////////////////////////////////////////////////////////////
