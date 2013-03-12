#ifndef __RICO_STUDENTS_T_H__
#define __RICO_STUDENTS_T_H__

#include "RicoBasicRV.h"

/////////////////////////////////////////////////////////////////////////

class RicoStudents_t : public RicoBasicRV {
  static const double SIGMA_FACTOR = 4.0;

 private:
  int m_df;

 public:
 RicoStudents_t(int df) : RicoBasicRV(), m_df(df) {}
  ~RicoStudents_t() {}

 public: // Interface
  ddPair getSuggestedPlotRange() const;
  double CDF(double x) const;

 public: // C++ Only
  virtual string getDefinition() const;

}; // RicoStudents_t

/////////////////////////////////////////////////////////////////////////
#endif //__RICO_STUDENTS_T_H__

