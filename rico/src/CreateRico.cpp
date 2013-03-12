#include <vector>
#include <utility>  
#include <iostream>
using std::cout; using std::endl;
using std::left; using std::fixed; using std::right; using std::scientific;
#include <iomanip>
using std::setw;
using std::setprecision;

#include "Rico.h"
#include "RicoRegistry.h"

#include "RicoInteger.h"
#include "RicoDouble.h"
#include "RicoFraction.h"
#include "RicoNegate.h"
#include "RicoSqrt.h"
#include "RicoLog.h"
#include "RicoExp.h"
#include "RicoAdd.h"
#include "RicoSubtract.h"
#include "RicoMultiply.h"
#include "RicoDivide.h"
#include "RicoPower.h"
#include "RicoTriangular.h"
#include "RicoNormal.h"
#include "RicoChiSquared.h"
#include "RicoStudents_t.h"
#include "RicoFisher.h"
#include "RicoCauchy.h"
#include "RicoLogNormal.h"
#include "RicoBeta.h"
#include "RicoExponential.h"
#include "RicoLogistic.h"
#include "RicoWeibull.h"
#include "RicoUniform.h"

#include "RicoSimplify.h"

#include <R.h>
#include <Rdefines.h>

#include "UnitTest++/UnitTest++.h"

//#include "Util.h"

// R interface ///////////////////////////////////////////////////////////////

extern "C" {

  int R_to_int(SEXP i_R) {
    int i;

    PROTECT(i_R = AS_INTEGER(i_R));
    i = INTEGER_POINTER(i_R)[0];
    UNPROTECT(1);

    return i;
  }

  SEXP int_to_R(int i) {
    SEXP i_R;

    PROTECT(i_R = NEW_INTEGER(1));
    INTEGER_POINTER(i_R)[0] = i;
    UNPROTECT(1);

    return i_R;
  }

  SEXP double_to_R(double d) {
    SEXP d_R;

    PROTECT(d_R = NEW_NUMERIC(1));
    NUMERIC_POINTER(d_R)[0] = d;
    UNPROTECT(1);

    return d_R;
  }

  SEXP bool_to_R(bool b) {
    SEXP b_R;

    PROTECT(b_R = NEW_LOGICAL(1));
    LOGICAL_POINTER(b_R)[0] = b;
    UNPROTECT(1);

    return b_R;    
  }

  double R_to_double(SEXP d_R) {
    double d;
    
    PROTECT(d_R = AS_NUMERIC(d_R)); 
    d = NUMERIC_POINTER(d_R)[0];
    UNPROTECT(1);

    return d;
  }

  vector<double> R_to_double_vector(SEXP v_R) {
    vector<double> v;
    R_len_t i,N = length(v_R);

    PROTECT(v_R = AS_NUMERIC(v_R)); 
    for(i = 0; i < N; i++) 
      v.push_back( NUMERIC_POINTER(v_R)[i] );
    UNPROTECT(1);

    return v;
  }
  
  SEXP double_vector_to_R(vector<double> x) {
    SEXP x_R;
    int N = x.size();
    PROTECT(x_R = allocVector(REALSXP,N));
    for(int i = 0; i < N; i++) REAL(x_R)[i] = x[i];
    UNPROTECT(1); 
    return x_R; 
  }

  SEXP Pair_to_R(ddPair pair) {
    SEXP x;
    PROTECT(x = allocVector(REALSXP,2));
    REAL(x)[0] = pair.first;
    REAL(x)[1] = pair.second;
    UNPROTECT(1); 
    return x; 
  }

  /// Creation and Util ///

  static RicoPtr RICO(SEXP id) {
    return RicoRegistry::getRicoPtr(R_to_int(id));
  }

  SEXP getIdPtr(RicoPtr rico) {
    return int_to_R(RicoRegistry::registerRicoPtr(rico)); // This is where we issue the id
  }

  // Numbers //////////////////////////////////////////////////////////

  SEXP create_Integer(SEXP i_R) {
    return getIdPtr(Rico::Integer(R_to_int(i_R)));
  }

  SEXP create_Double(SEXP x_R) {
    return getIdPtr(Rico::Double(R_to_double(x_R)));
  }

  // Distributions ///////////////////////////////////////////////////

  SEXP create_Triangular() {
    return getIdPtr(Rico::Triangular());
  }

  SEXP create_TriangularM(SEXP Middle) {
    double middle = R_to_double(Middle);
    return getIdPtr(Rico::Triangular(middle));
  }

  SEXP create_TriangularLR(SEXP Left, SEXP Right) {
    double left  = R_to_double(Left);
    double right = R_to_double(Right);
    return getIdPtr(Rico::Triangular(left, right));
  }

  SEXP create_TriangularLMR(SEXP Left, SEXP Middle, SEXP Right) {
    double left   = R_to_double(Left);
    double middle = R_to_double(Middle);
    double right  = R_to_double(Right);
    return getIdPtr(Rico::Triangular(left, middle, right));
  }

  SEXP create_Normal() {
    return getIdPtr(Rico::Normal());
  }

  SEXP create_Normal_Mean(SEXP Mean) {
    return getIdPtr(Rico::Normal(R_to_double(Mean)));
  }

  SEXP create_Normal_Mean_Stdev(SEXP Mean, SEXP StdDev) {
    return getIdPtr(Rico::Normal(R_to_double(Mean), R_to_double(StdDev)));
  }

  SEXP create_ChiSquared(SEXP DF) {
    return getIdPtr(Rico::ChiSquared(R_to_int(DF)));
  }

  SEXP create_Students_t(SEXP DF) {
    return getIdPtr(Rico::Students_t(R_to_int(DF)));
  }

  SEXP create_Fisher_F(SEXP DF1, SEXP DF2) {
    return getIdPtr(Rico::Fisher(R_to_int(DF1), R_to_int(DF2)));
  }

  SEXP create_Cauchy() {
    return getIdPtr(Rico::Cauchy());
  }

  SEXP create_Cauchy_Mean(SEXP Mean) {
    return getIdPtr(Rico::Cauchy(R_to_double(Mean)));
  }

  SEXP create_Cauchy_Mean_Stdev(SEXP Mean, SEXP StdDev) {
    return getIdPtr(Rico::Cauchy(R_to_double(Mean), R_to_double(StdDev)));
  }

  SEXP create_LogNormal() {
    return getIdPtr(Rico::LogNormal());
  }

  SEXP create_LogNormal_Mean(SEXP Mean) {
    return getIdPtr(Rico::LogNormal(R_to_double(Mean)));
  }

  SEXP create_LogNormal_Mean_Stdev(SEXP Mean, SEXP StdDev) {
    return getIdPtr(Rico::LogNormal(R_to_double(Mean), R_to_double(StdDev)));
  }

  SEXP create_Beta(SEXP Alpha, SEXP Beta) {
    return getIdPtr(Rico::Beta(R_to_double(Alpha), R_to_double(Beta)));
  }

  SEXP create_Exponential(SEXP Lambda) {
    return getIdPtr(Rico::Exponential(R_to_double(Lambda)));
  }

  SEXP create_Logistic(SEXP Mean, SEXP Scale) {
    return getIdPtr(Rico::Logistic(R_to_double(Mean), R_to_double(Scale)));
  }

  SEXP create_Weibull_Shape(SEXP Shape) {
    return getIdPtr(Rico::Weibull(R_to_double(Shape)));
  }

  SEXP create_Weibull_Shape_Scale(SEXP Shape, SEXP Scale) {
    return getIdPtr(Rico::Weibull(R_to_double(Shape), R_to_double(Scale)));
  }

  SEXP create_Uniform_Std() {
    return getIdPtr(Rico::Uniform());
  }

  SEXP create_Uniform_Left_Right(SEXP Left, SEXP Right) {
    return getIdPtr(Rico::Uniform(R_to_double(Left), R_to_double(Right)));
  }

  // R/C-Operations ///////////////////////////////////////////////////

  SEXP createRico_Number(SEXP X) {
    switch(TYPEOF(X)) {
    case REALSXP:
      {
	double x = R_to_double(X);
	int    i = (int) x;
	if( (double)i == x ) {
	  //Rprintf("(int) %d\n", i);
	  return getIdPtr(Rico::Integer(i));
	} else {
	  //Rprintf("(real) %f\n", x);
	  return getIdPtr(Rico::Double(x));
	}
      }
    case LGLSXP:
    case INTSXP:
      //Rprintf("(int) %d\n", INTEGER(X)[0]);
      return getIdPtr(Rico::Integer(R_to_int(X)));
    case CPLXSXP:
      {
	Rcomplex cpl = COMPLEX(X)[0];
	Rprintf("(complex) %f + %fi\n", cpl.r, cpl.i);
      }
      break;
    case STRSXP:
      Rprintf("(string) %s\n", CHAR(STRING_ELT(X, 0)));
      break;
    default:
      Rprintf("(unknown) R type\n");
    }
    return(R_NilValue);
  }

  SEXP createRico_Negate(SEXP Xid) {
    return getIdPtr(Rico::Negate(RICO(Xid)));
  }

  SEXP createRico_Sqrt(SEXP Xid) {
    return getIdPtr(Rico::Sqrt(RICO(Xid)));
  }

  SEXP createRico_Log(SEXP Xid) {
    return getIdPtr(Rico::Log(RICO(Xid)));
  }

  SEXP createRico_Exp(SEXP Xid) {
    return getIdPtr(Rico::Exp(RICO(Xid)));
  }

  SEXP createRico_Add(SEXP Xid, SEXP Yid) {
    return getIdPtr(Rico::Add(RICO(Xid), RICO(Yid)));
  }

  SEXP createRico_Subtract(SEXP Xid, SEXP Yid) {
    return getIdPtr(Rico::Subtract(RICO(Xid), RICO(Yid)));
  }

  SEXP createRico_Multiply(SEXP Xid, SEXP Yid) {
    return getIdPtr(Rico::Multiply(RICO(Xid), RICO(Yid)));
  }

  SEXP createRico_Divide(SEXP Xid, SEXP Yid) {
    return getIdPtr(Rico::Divide(RICO(Xid), RICO(Yid)));
  }

  SEXP createRico_Power(SEXP Xid, SEXP Yid) {
    return getIdPtr(Rico::Power(RICO(Xid), RICO(Yid)));
  }

  // R/C-Interface ///////////////////////////////////////////////////

  void RicoTestRun() {
    UnitTest::RunAllTests();
  }

  SEXP getSuggestedPlotRange(SEXP id_R) {
    return Pair_to_R( RICO(id_R)->getSuggestedPlotRange() );
  }

  SEXP getPlotPoints(SEXP id_R, SEXP X) {
    vector<double> x = R_to_double_vector(X);
    vector<double> h = RICO(id_R)->Density(x) ;
    return double_vector_to_R(h);
  }

  void show_Rico(SEXP id) {
    string str = RICO(id)->getDefinedExpression();
    Rprintf("%s", str.data());
  }

  SEXP differentiateRico_X(SEXP EXP_id) {
    RicoPtr This = RICO(EXP_id);
    return getIdPtr(This->differentiate_single(This));
  }

  SEXP differentiateRico_WRT(SEXP EXP_id, SEXP WRT_id) {
    RicoPtr This = RICO(EXP_id);
    RicoPtr Wrt  = RICO(WRT_id);
    return getIdPtr(This->differentiate(This, Wrt));
  }

  SEXP compareRico(SEXP Xid, SEXP Yid) {
    RicoPtr X = RICO(Xid);
    RicoPtr Y = RICO(Yid);
    return bool_to_R(X == Y);
  }

  SEXP simplifyRico(SEXP Xid) {
    RicoPtr X = RICO(Xid);
    return getIdPtr(RicoSimplify::Simplify(X));
  }

} // extern "C"

