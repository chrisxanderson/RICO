#ifndef __RICO_H__
#define __RICO_H__

#include <iostream>
#include <sstream>
#include <string>

#include <map>
#include <set>
#include <vector>
#include <list>
#include <utility>  

#include "SmartPtr.h"

using namespace std;

// Forwards //////////////////////////////////////////////////////////////////////////////

class Rico;
class RicoFunction;

typedef SmartPtr<const Rico>   RicoPtr;
typedef list<RicoPtr>         vRicoPtr; 
typedef int                    RicoID;
typedef pair<int, int>         iiPair;
typedef pair<double, double>   ddPair;
typedef pair<RicoPtr, RicoPtr> ppPair;
typedef map<int, int>          iiMap;
//////////////////////////////////////////////////////////////////////////////////////////

class Rico {
 protected:
  enum RICO_TYPE {BASIC_RV, NUMBER, FUNCTION};

 private:
  RICO_TYPE m_type;

 public:
  Rico(RICO_TYPE type) : m_type(type) {}
  virtual ~Rico() {}

 public: // Static Constructors
  static int registerRicoPtr(RicoPtr x);

  static RicoPtr NaN;
  static RicoPtr NegOne;
  static RicoPtr Zero;
  static RicoPtr Half;
  static RicoPtr One;
  static RicoPtr Two;
  static RicoPtr Three;

  static RicoPtr Integer(int x);
  static RicoPtr Double(); // NaN
  static RicoPtr Double(double x);
  static RicoPtr Fraction(int x, int y);

  static RicoPtr Negate          (RicoPtr x);
  static RicoPtr Sqrt            (RicoPtr x);
  static RicoPtr Log             (RicoPtr x);
  static RicoPtr Exp             (RicoPtr x);
  static RicoPtr Square          (RicoPtr x);
  static RicoPtr Inverse         (RicoPtr x);

  static RicoPtr Add(RicoPtr x, RicoPtr y);
  static RicoPtr Add(RicoPtr x, RicoPtr y, RicoPtr z);
  static RicoPtr Add(RicoPtr x, RicoPtr y, RicoPtr z, RicoPtr k);
  static RicoPtr Add(const vRicoPtr& v);
  static void    Add(vRicoPtr& V, RicoPtr A); // internal

  static RicoPtr Subtract(RicoPtr x, RicoPtr y);
  static RicoPtr Subtract(const vRicoPtr& v);

  static RicoPtr Multiply(RicoPtr x, RicoPtr y);
  static RicoPtr Multiply(RicoPtr x, RicoPtr y, RicoPtr z);
  static RicoPtr Multiply(RicoPtr x, RicoPtr y, RicoPtr z, RicoPtr k);
  static RicoPtr Multiply(const vRicoPtr& v);
  static void    Multiply(vRicoPtr& V, RicoPtr A); // internal
  
  static RicoPtr Function(const RicoFunction& f, const vRicoPtr& v);

  static RicoPtr Divide     (RicoPtr x, RicoPtr y);
  static RicoPtr Power      (RicoPtr x, RicoPtr y);

  static RicoPtr Triangular();
  static RicoPtr Triangular(double middle);
  static RicoPtr Triangular(double left, double right);
  static RicoPtr Triangular(double left, double middle, double right);

  static RicoPtr Normal();
  static RicoPtr Normal(double mu);
  static RicoPtr Normal(double mu, double sigma);

  static RicoPtr LogNormal();
  static RicoPtr LogNormal(double mu);
  static RicoPtr LogNormal(double mu, double sigma);

  static RicoPtr Cauchy();
  static RicoPtr Cauchy(double mu);
  static RicoPtr Cauchy(double mu, double sigma);

  static RicoPtr ChiSquared(int df);
  static RicoPtr Students_t(int df);
  static RicoPtr Fisher(int df1, int df2);
  static RicoPtr Beta(double alpha, double beta);
  static RicoPtr Exponential(double lambda);
  static RicoPtr Logistic(double mu, double s);
  static RicoPtr Weibull(double shape);
  static RicoPtr Weibull(double shape, double scale);
  static RicoPtr Uniform();
  static RicoPtr Uniform(double left, double right);

 public:
  virtual RicoID id()  const {return 0;}

  RICO_TYPE type()  const {return m_type;}
  bool isBasicRV()  const {return m_type == BASIC_RV;}
  bool isNumber()   const {return m_type == NUMBER;}
  bool isFunction() const {return m_type == FUNCTION;}

  bool isInteger () const;
  bool isFraction() const;
  bool isDouble  () const;
  bool isZero    () const;
  bool isPos     () const;
  bool isNeg     () const;
  bool isOne     () const;
  bool isNaN     () const;
  bool isPosInf  () const;
  bool isNegInf  () const;

  bool isNegate  () const;
  bool isAdd     () const;
  bool isSubtract() const;
  bool isMultiply() const;
  bool isDivide  () const;
  bool isPower   () const;
  bool isLog     () const;
  bool isExp     () const;
  bool isSqrt    () const;

 public: // Parameter access
  RicoPtr x() const;
  RicoPtr y() const;

  virtual ppPair getPair() const = 0;

 public: // R/C-interface
  string getDefinedExpression() const; // Assumed: top-level call.

  vector<double> CDF    (vector<double> X) const;
  vector<double> Pr     (vector<double> x) const;
  vector<double> Density(vector<double> x) const;  // the plot function

  const vRicoPtr differentiate(       RicoPtr This) const;  // i.e. with respect to each
         RicoPtr differentiate_single(RicoPtr This) const;  // i.e. with respect to (implied)
	 RicoPtr differentiate(RicoPtr This, RicoPtr Wrt) const;

 public: // Operations
  virtual string  getExpression()                         const = 0;
  virtual ddPair  getSuggestedPlotRange()                 const {return make_pair(0,0);} // STUB
  virtual double  CDF(double x)                           const {return 0;} // STUB
  virtual RicoPtr differentiate(RicoPtr This, RicoID Wrt) const = 0;

 public: // utility
  virtual void  collectBasicRVs(set<RicoID>& rv_set) const {}

  static vector<double> PointsToPartition(vector<double>);

};

////////////////////////////////////////////////////////////////////////////////////
#endif // __RICO_H__
