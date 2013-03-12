#ifndef __RICO_FUNCTION_H__
#define __RICO_FUNCTION_H__

#include "Rico.h"

#include <vector>

/////////////////////////////////////////////////////////////////////////

class RicoFunction : public Rico {
 protected:
  enum FTYPE {NEG, ADD, SUB, MUL, DIV, POW, LOG, EXP, SQRT};

 private:
  FTYPE    m_ftype;
  vRicoPtr m_params;

 public:
 RicoFunction(FTYPE ftype, RicoPtr x) : Rico(Rico::FUNCTION), m_ftype(ftype) {
    m_params.push_back(x);
  }
 RicoFunction(FTYPE ftype, RicoPtr x, RicoPtr y) : Rico(Rico::FUNCTION), m_ftype(ftype) {
    m_params.push_back(x);
    m_params.push_back(y);
  }
 RicoFunction(FTYPE ftype, RicoPtr x, RicoPtr y, RicoPtr z) : Rico(Rico::FUNCTION), m_ftype(ftype) {
    m_params.push_back(x);
    m_params.push_back(y);
    m_params.push_back(z);
  }
  RicoFunction(FTYPE ftype, const vRicoPtr& v) : Rico(Rico::FUNCTION),m_params(v),m_ftype(ftype) {}
  ~RicoFunction() {}

  static RicoFunction* Function(const RicoFunction& f, const vRicoPtr& v);

 public: // Access
  FTYPE           ftype()      const {return m_ftype;}
  RicoPtr         x()          const {return m_params.front();}
  RicoPtr         y()          const {return *(++m_params.begin());}
  RicoPtr         z()          const {return m_params.back();}
  int             size()       const {return m_params.size();}
  const vRicoPtr& get_params() const {return m_params;}

  void params(vRicoPtr& P)        const;
  void params(vRicoPtr& P, int i) const; // all params except i'th

  bool isNegate()   const {return ftype() == NEG;}
  bool isAdd()      const {return ftype() == ADD;}
  bool isSubtract() const {return ftype() == SUB;}
  bool isMultiply() const {return ftype() == MUL;}
  bool isDivide()   const {return ftype() == DIV;}
  bool isPower()    const {return ftype() == POW;}
  bool isLog()      const {return ftype() == LOG;}
  bool isExp()      const {return ftype() == EXP;}
  bool isSqrt()     const {return ftype() == SQRT;}

  virtual ppPair getPair() const {return ppPair(x(), y());}

 public: // Operations
  virtual bool operator==(RicoPtr That) const;

 public: // Utility
  virtual void collectBasicRVs(set<RicoID>& rv_set) const;

  vRicoPtr differentiate_params(RicoID Wrt) const;

 public: // Testing
  bool operator==(const RicoFunction& that) const;
  static bool test_membership(RicoPtr ai, vRicoPtr& B);

}; // RicoFunction

/////////////////////////////////////////////////////////////////////////
#endif //__RICO_FUNCTION_H__
