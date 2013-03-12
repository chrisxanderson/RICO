#include "RicoExponential.h"

#include <iostream>
#include <sstream>
#include <string>

#include <cmath>

// RicoExponential //////////////////////////////////////////////////////////////////

ddPair RicoExponential::getSuggestedPlotRange() const {
  return make_pair(0, -log(RANGE_LIMIT)/m_lambda);
}

double RicoExponential::CDF(double x) const {
  // Pr(X <= x)
  if(x <= 0) return 0;
  return 1 - exp(-m_lambda * x);
}

string RicoExponential::getDefinition() const {
  stringstream sstr;
  sstr << "Exp[" << m_lambda << "]";
  return sstr.str();
}

////////////////////////////////////////////////////////////////////////////////////
