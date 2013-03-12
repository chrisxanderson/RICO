#include "RicoChiSquared.h"

#include <iostream>
#include <sstream>
#include <string>

#include "gsl.h"

// RicoChiSquared Static ////////////////////////////////////////////////////////////

const double RicoChiSquared::RANGE_DF[] = {6.635,9.21,11.345,13.277,15.086,16.812,18.475,20.09, \
					   21.666,23.209,24.725,26.217,27.688,29.141,30.578,32,\
					   33.409,34.805,36.191,37.566,38.932,40.289,41.638,42.98,\
					   44.314,45.642,46.963,48.278,49.588,50.892};

// RicoChiSquared //////////////////////////////////////////////////////////////////

ddPair RicoChiSquared::getSuggestedPlotRange() const {
    int N = sizeof(RANGE_DF)/sizeof(double);
    if(m_df <= N) return make_pair(0,RANGE_DF[m_df]);
    double y1 = RANGE_DF[N-2];
    double y2 = RANGE_DF[N-1];
    double R  = (y2-y1)*(m_df - N) + RANGE_DF[N-1];
    return make_pair(0, R*(1+RANGE_INFLATE));
  }

double RicoChiSquared::CDF(double x) const {
  // Pr(X <= x)
  if(x <= 0) return 0;
  return gsl_sf_gamma_inc_P (m_df/2.0, x/2.0);
}

string RicoChiSquared::getDefinition() const {
  stringstream sstr;
  sstr << "Chisq[" << m_df << "]";
  return sstr.str();
}

/////////////////////////////////////////////////////////////////////////////////////////
