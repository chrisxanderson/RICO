#include "RicoCauchy.h"

#include <iostream>
#include <sstream>
#include <string>

#include <cmath>

// RicoCauchy //////////////////////////////////////////////////////////////////

const double RicoCauchy::PI = 4 * atan(1);

// Interface ///////////////////////////////////////////////////////////////////////

ddPair RicoCauchy::getSuggestedPlotRange() const {
  return make_pair(m_mu-SIGMA_FACTOR*m_sigma, m_mu+SIGMA_FACTOR*m_sigma);
}

double RicoCauchy::CDF(double x) const {
  // Pr(X <= x)
  return 1/PI * atan((x-m_mu)/m_sigma) + 0.5;
}

string RicoCauchy::getDefinition() const {
  stringstream sstr;
  sstr << "Cauchy[" << m_mu << ", " << m_sigma << "]";
  return sstr.str();
}

/////////////////////////////////////////////////////////////////////////////////
