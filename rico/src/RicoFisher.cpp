#include "RicoFisher.h"

#include <iostream>
#include <sstream>
#include <string>

#include "gsl.h"

// RicoFisher //////////////////////////////////////////////////////////////////

ddPair RicoFisher::getSuggestedPlotRange() const {
  return make_pair(0,3);
}

double RicoFisher::CDF(double x) const {
  // Pr(X <= x)
  if(x <= 0) return 0;
  double z   = m_df1 * x / (m_df1 * x + m_df2);
  return gsl_sf_beta_inc(m_df1/2.0, m_df2/2.0, z);
}

string RicoFisher::getDefinition() const {
  stringstream sstr;
  sstr << "F[" << m_df1 << ", " << m_df2 << "]";
  return sstr.str();
}

/////////////////////////////////////////////////////////////////////////////////
