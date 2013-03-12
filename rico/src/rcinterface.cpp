#include <R.h>
#include <Rdefines.h>

/*  This collection of functions is just for example purposes.
 *  I found it in R-C-interface.ppt
 */

extern "C" {

  SEXP getInt(SEXP myint, SEXP myintVar) {         // SEXP; Simple EXPression. myint is int*
    int Imyint, n;
    int *Pmyint; 
    PROTECT(myint = AS_INTEGER(myint));	 // R objects create in C must be protected.
    Imyint = INTEGER_POINTER(myint)[0];         // grab first value
    Pmyint = INTEGER_POINTER(myint);
    n      = INTEGER_VALUE(myintVar);
    printf(" TRF: Printed from C: \n");
    printf(" TRF: Imyint: %d \n", Imyint);
    printf(" TRF: n: %d \n", n);
    printf(" TRF: Pmyint[0], Pmyint[1]: %d %d \n", Pmyint[0], Pmyint[1]);
    Pmyint[0] = 5;                              // No effect.
    UNPROTECT(1);                               // un-protect one object back. i.e. pop stack.
    return(R_NilValue);
  }

  SEXP getChar(SEXP mychar) {
    char *Pmychar[2];
    
    PROTECT(mychar = AS_CHARACTER(mychar));
    
    Pmychar[0] = R_alloc(strlen(CHAR(STRING_ELT(mychar, 0))), sizeof(char)); 
    Pmychar[1] = R_alloc(strlen(CHAR(STRING_ELT(mychar, 1))), sizeof(char)); 
    
    strcpy(Pmychar[0], CHAR(STRING_ELT(mychar, 0))); 
    strcpy(Pmychar[1], CHAR(STRING_ELT(mychar, 1))); 
    
    printf(" Printed from C:");
    printf(" %s %s \n",Pmychar[0],Pmychar[1]);
    
    UNPROTECT(1);
    
    return(R_NilValue); 
  }

  SEXP setInt() {
    SEXP myint;
    int *p_myint;
    int len = 5;
    
    PROTECT(myint = NEW_INTEGER(len));
    
    p_myint    = INTEGER_POINTER(myint);
    p_myint[0] = 7;
    
    UNPROTECT(1);
    
    return myint;
  }

  SEXP setChar() {
    SEXP mychar;
    PROTECT(mychar = allocVector(STRSXP, 5));
    SET_STRING_ELT(mychar, 0, mkChar("A")); 
    UNPROTECT(1);
    return mychar;
  }

  SEXP setList() {
    int    *p_myint, i; 
    double *p_double;
    char   *names[2] = {"integer", "numeric"};
    
    SEXP mydouble, myint, list, list_names;   
    
    PROTECT(myint = NEW_INTEGER(5)); 
    p_myint = INTEGER_POINTER(myint);
    
    PROTECT(mydouble = NEW_NUMERIC(5)); 
    p_double = NUMERIC_POINTER(mydouble);
    
    for(i = 0; i < 5; i++) {
      p_double[i] = 1/(double)(i + 1);
      p_myint[i] = i + 1;
    }
    
    PROTECT(list_names = allocVector(STRSXP,2));
    for(i = 0; i < 2; i++)   
      SET_STRING_ELT(list_names,i,mkChar(names[i])); 
    
    PROTECT(list = allocVector(VECSXP, 2)); 
    SET_VECTOR_ELT(list, 0, myint); 
    SET_VECTOR_ELT(list, 1, mydouble); 
    
    setAttrib(list, R_NamesSymbol, list_names); 
    
    UNPROTECT(4);
    return list;
  }

} // extern "C"
