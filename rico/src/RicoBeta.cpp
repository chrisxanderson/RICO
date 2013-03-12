#include "RicoBeta.h"

#include <iostream>
#include <sstream>
#include <string>

#include "gsl.h"

// RicoBeta ////////////////////////////////////////////////////////////////////////

double RicoBeta::CDF(double x) const {
  // Pr(X <= x)
  if(x <= 0) return 0;
  if(x >= 1) return 1;
  return gsl_sf_beta_inc(m_alpha, m_beta, x);
}

string RicoBeta::getDefinition() const {
  stringstream sstr;
  sstr << "Beta[" << m_alpha << ", " << m_beta << "]";
  return sstr.str();
}

/////////////////////////////////////////////////////////////////////////////////////
