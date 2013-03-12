#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>

#include "Rico.h"
#include "Util.h"

#include <R.h>
#include <Rdefines.h>

using namespace std;

/// R/C Interface /////////////////////////////////////////////////////////////////////////////////

extern "C" {

  // Testing Only /////////////////////////////////////////////////////

  SEXP showArgs(SEXP args)
  {
    args = CDR(args); /* skip ’name’ */
    for(int i = 0; args != R_NilValue; i++, args = CDR(args)) {
      const char *name = isNull(TAG(args)) ? "" : CHAR(PRINTNAME(TAG(args)));
      SEXP el = CAR(args);
      if (length(el) == 0) {
	Rprintf("[%d](empty) ’%s’ R type, length 0\n", i+1, name);
	continue;
      }
      switch(TYPEOF(el)) {
      case REALSXP:
	Rprintf("[%d](real) ’%s’ %f\n", i+1, name, REAL(el)[0]);
	break;
      case LGLSXP:
      case INTSXP:
	Rprintf("[%d](int) ’%s’ %d\n", i+1, name, INTEGER(el)[0]);
	break;
      case CPLXSXP:
	{
	  Rcomplex cpl = COMPLEX(el)[0];
	  Rprintf("[%d] ’%s’ %f + %fi\n", i+1, name, cpl.r, cpl.i);
	}
	break;
      case STRSXP:
	Rprintf("[%d](string) ’%s’ %s\n", i+1, name, CHAR(STRING_ELT(el, 0)));
	break;
      default:
	Rprintf("[%d](unknown) ’%s’ R type\n", i+1, name);
      }
    }
    return(R_NilValue);
  }

  SEXP showArgs1(SEXP largs)
  {
    int i, nargs = LENGTH(largs);
    Rcomplex cpl;
    SEXP el, names = getAttrib(largs, R_NamesSymbol);
    const char *name;
    
    for(i = 0; i < nargs; i++) {
      el = VECTOR_ELT(largs, i);
      name = isNull(names) ? "" : CHAR(STRING_ELT(names, i));
      switch(TYPEOF(el)) {
      case REALSXP:
	Rprintf("[%d](real) '%s' %f\n", i+1, name, REAL(el)[0]);
	break;
      case LGLSXP:
      case INTSXP:
	Rprintf("[%d](int) '%s' %d\n", i+1, name, INTEGER(el)[0]);
	break;
      case CPLXSXP:
	cpl = COMPLEX(el)[0];
	Rprintf("[%d](complex) '%s' %f + %fi\n", i+1, name, cpl.r, cpl.i);
	break;
      case STRSXP:
	Rprintf("[%d](string) '%s' %s\n", i+1, name, CHAR(STRING_ELT(el, 0)));
	break;
      default:
	Rprintf("[%d](unknown) '%s' R type\n", i+1, name);
      }
    }
    return(R_NilValue);
  }
  
  SEXP double_to_R_test(double d) {
    SEXP d_R;

    PROTECT(d_R = NEW_NUMERIC(1));
    NUMERIC_POINTER(d_R)[0] = d;
    UNPROTECT(1);

    return d_R;
  }

  SEXP int_to_R_test(int i) {
    SEXP i_R;

    PROTECT(i_R = NEW_INTEGER(1));
    INTEGER_POINTER(i_R)[0] = i;
    UNPROTECT(1);

    return i_R;
  }

  SEXP testGetList() {
    SEXP x, i;
    PROTECT(x = i = allocList(2));
    SETCAR(i, double_to_R_test(7)); i = CDR(i);
    SETCAR(i, int_to_R_test(5));    i = CDR(i);
    UNPROTECT(1); 
    return x; 
  }

  SEXP testGetArray() {
    static int N = 100;
    SEXP x;
    PROTECT(x = allocVector(REALSXP,N));
    for(int i = 0; i < N; i++) {
      REAL(x)[i] = i*i;
    }
    UNPROTECT(1); 
    return x; 
  }

  SEXP testGetX() {
    static int N = 100;
    SEXP x;
    PROTECT(x = allocVector(REALSXP,N));
    for(int i = 0; i < N; i++) {
      REAL(x)[i] = i;
    }
    UNPROTECT(1); 
    return x; 
  }

  SEXP testGetY() {
    static int N = 100;
    SEXP x;
    PROTECT(x = allocVector(REALSXP,N));
    for(int i = 0; i < N; i++) {
      REAL(x)[i] = i*i;
    }
    UNPROTECT(1); 
    return x; 
  }

  SEXP testGetArrayList() {
    SEXP x, i;
    PROTECT(x = i = allocList(2));
    SETCAR(i, testGetX());  i = CDR(i);
    SETCAR(i, testGetY());  i = CDR(i);
    UNPROTECT(1); 
    return x; 
  }

  // This is pulled from R/stats package file nls.c and from R-exts.pdf
  // It shows how to grab elements from a list by name
  SEXP getListElement(SEXP list, const char *str)
  {
    SEXP elmt = (SEXP) NULL, names = getAttrib(list, R_NamesSymbol);

    for (int i = 0; i < length(list); i++) {
      const char *tempChar = CHAR(STRING_ELT(names, i));     // ASCII only 
      if( strcmp(tempChar,str) == 0) {
	elmt = VECTOR_ELT(list, i);
	break;
      }
    }
    return elmt;
  }

} // extern "C"

//////////////////////////////////////////////////////////////////////////////////////////////////
