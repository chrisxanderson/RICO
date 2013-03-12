#include "RicoStudents_t.h"

#include <iostream>
#include <sstream>
#include <string>

#include "gsl.h"

// RicoStudents_t //////////////////////////////////////////////////////////////////

ddPair RicoStudents_t::getSuggestedPlotRange() const {
  return make_pair(-SIGMA_FACTOR, SIGMA_FACTOR);
}

double RicoStudents_t::CDF(double x) const {
  // Pr(X <= x)
   if (x == 0) return 0.5;
   double x2 = x * x;
   double z  = m_df / (m_df + x2);
   double probability = gsl_sf_beta_inc(m_df/2.0, 0.5, z) / 2;
   return (x > 0 ? 1 - probability : probability);
}

string RicoStudents_t::getDefinition() const {
  stringstream sstr;
  sstr << "Student[" << m_df << "]";
  return sstr.str();
}

////////////////////////////////////////////////////////////////////////
