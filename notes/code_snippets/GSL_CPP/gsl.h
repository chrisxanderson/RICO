/* gsl.h 
 * This file is intended to be the unified header file for the subset of GSL files I need for Rico.
 */

#include <math.h>
#include <limits.h>
#include <float.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
//#include <ieeefp.h>

// acconfig.h ////////////////////////////////////////////////////////////////////
//#include "acconfig.h"         // DONE was <config.h>
/* Defined if this is an official release */
#undef RELEASED

/* Define if you have inline */
#undef HAVE_INLINE

/* Define if you need to hide the static definitions of inline functions */
#undef HIDE_INLINE_STATIC

/* Defined if you have ansi EXIT_SUCCESS and EXIT_FAILURE in stdlib.h */
#undef HAVE_EXIT_SUCCESS_AND_FAILURE

/* Use 0 and 1 for EXIT_SUCCESS and EXIT_FAILURE if we don't have them */
#if !HAVE_EXIT_SUCCESS_AND_FAILURE
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1
#endif

/* Define this if printf can handle %Lf for long double */
#undef HAVE_PRINTF_LONGDOUBLE

/* Define one of these if you have a known IEEE arithmetic interface */
#undef HAVE_GNUSPARC_IEEE_INTERFACE
#undef HAVE_GNUM68K_IEEE_INTERFACE
#undef HAVE_GNUPPC_IEEE_INTERFACE
#undef HAVE_GNUX86_IEEE_INTERFACE
#undef HAVE_SUNOS4_IEEE_INTERFACE
#undef HAVE_SOLARIS_IEEE_INTERFACE
#undef HAVE_HPUX11_IEEE_INTERFACE
#undef HAVE_HPUX_IEEE_INTERFACE
#undef HAVE_TRU64_IEEE_INTERFACE
#undef HAVE_IRIX_IEEE_INTERFACE
#undef HAVE_AIX_IEEE_INTERFACE
#undef HAVE_FREEBSD_IEEE_INTERFACE
#undef HAVE_OS2EMX_IEEE_INTERFACE
#undef HAVE_NETBSD_IEEE_INTERFACE
#undef HAVE_OPENBSD_IEEE_INTERFACE
#undef HAVE_DARWIN_IEEE_INTERFACE
#undef HAVE_DARWIN86_IEEE_INTERFACE

/* Define this if IEEE comparisons work correctly (e.g. NaN != NaN) */
//#undef HAVE_IEEE_COMPARISONS
#define HAVE_IEEE_COMPARISONS 1 // TRF: I hope this is true

/* Define this if IEEE denormalized numbers are available */
#undef HAVE_IEEE_DENORMALS

/* Define a rounding function which moves extended precision values
   out of registers and rounds them to double-precision. This should
   be used *sparingly*, in places where it is necessary to keep
   double-precision rounding for critical expressions while running in
   extended precision. For example, the following code should ensure
   exact equality, even when extended precision registers are in use,

      double q = GSL_COERCE_DBL(3.0/7.0) ;
      if (q == GSL_COERCE_DBL(3.0/7.0)) { ... } ;

   It carries a penalty even when the program is running in double
   precision mode unless you compile a separate version of the
   library with HAVE_EXTENDED_PRECISION_REGISTERS turned off. */

#undef HAVE_EXTENDED_PRECISION_REGISTERS

#if HAVE_EXTENDED_PRECISION_REGISTERS
#define GSL_COERCE_DBL(x) (gsl_coerce_double(x))
#else
#define GSL_COERCE_DBL(x) (x)
#endif

/* Substitute gsl functions for missing system functions */

/* TRF: Have this 
#if !HAVE_DECL_HYPOT
#define hypot gsl_hypot
#endif
*/

/* TRF: Have this 
#if !HAVE_DECL_LOG1P
#define log1p gsl_log1p
#endif
*/

#if !HAVE_DECL_EXPM1
#define expm1 gsl_expm1
#endif

/* TRF: Have this
#if !HAVE_DECL_ACOSH
#define acosh gsl_acosh
#endif
*/

#if !HAVE_DECL_ASINH
#define asinh gsl_asinh
#endif

/* TRF: Have this
#if !HAVE_DECL_ATANH
#define atanh gsl_atanh
#endif
*/

/* TRF: Have this
#if !HAVE_DECL_LDEXP
#define ldexp gsl_ldexp
#endif
*/

#if !HAVE_DECL_FREXP
#define frexp gsl_frexp
#endif

/* TRF: We have this
#if !HAVE_DECL_ISINF
#define isinf gsl_isinf
#endif
*/

#if !HAVE_DECL_FINITE
#if HAVE_DECL_ISFINITE
#define finite isfinite
#else
#define finite gsl_finite
#endif
#endif

/* TRF: Have this
#if !HAVE_DECL_ISNAN
#define isnan gsl_isnan
#endif
*/

#ifdef __GNUC__
#define DISCARD_POINTER(p) do { ; } while(p ? 0 : 0);
#else
#define DISCARD_POINTER(p) /* ignoring discarded pointer */
#endif

#if defined(GSL_RANGE_CHECK_OFF) || !defined(GSL_RANGE_CHECK)
#define GSL_RANGE_CHECK 0  /* turn off range checking by default internally */
#endif

/* Disable deprecated functions and enums while building */
#define GSL_DISABLE_DEPRECATED 1
// </accconfig.h> 

///////////////////////////////////////////////////////////////////////////////////////
//#include "gsl_sys.h"          // DONE
extern "C" {
double gsl_log1p (const double x);
double gsl_expm1 (const double x);
double gsl_hypot (const double x, const double y);
double gsl_acosh (const double x);
double gsl_asinh (const double x);
double gsl_atanh (const double x);

int gsl_isnan (const double x);
int gsl_isinf (const double x);
int gsl_finite (const double x);

double gsl_nan (void);
double gsl_posinf (void);
double gsl_neginf (void);
double gsl_fdiv (const double x, const double y);

double gsl_coerce_double (const double x);
float gsl_coerce_float (const float x);
long double gsl_coerce_long_double (const long double x);

double gsl_ldexp(const double x, const int e);
double gsl_frexp(const double x, int * e);

int gsl_fcmp (const double x1, const double x2, const double epsilon);
} // gsl_sys.h
///////////////////////////////////////////////////////////////////////////////////////
//#include "gsl_machine.h"      // DONE
/* -*-MACHINE CONSTANTS-*-
 *
 * PLATFORM: Whiz-O-Matic 9000
 * FP_PLATFORM: IEEE-Virtual
 * HOSTNAME: nnn.lanl.gov
 * DATE: Fri Nov 20 17:53:26 MST 1998
 */
#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define GSL_SQRT_DBL_EPSILON   1.4901161193847656e-08
#define GSL_ROOT3_DBL_EPSILON  6.0554544523933429e-06
#define GSL_ROOT4_DBL_EPSILON  1.2207031250000000e-04
#define GSL_ROOT5_DBL_EPSILON  7.4009597974140505e-04
#define GSL_ROOT6_DBL_EPSILON  2.4607833005759251e-03
#define GSL_LOG_DBL_EPSILON   (-3.6043653389117154e+01)

#define GSL_DBL_MIN        2.2250738585072014e-308
#define GSL_SQRT_DBL_MIN   1.4916681462400413e-154
#define GSL_ROOT3_DBL_MIN  2.8126442852362996e-103
#define GSL_ROOT4_DBL_MIN  1.2213386697554620e-77
#define GSL_ROOT5_DBL_MIN  2.9476022969691763e-62
#define GSL_ROOT6_DBL_MIN  5.3034368905798218e-52
#define GSL_LOG_DBL_MIN   (-7.0839641853226408e+02)

#define GSL_DBL_MAX        1.7976931348623157e+308
#define GSL_SQRT_DBL_MAX   1.3407807929942596e+154
#define GSL_ROOT3_DBL_MAX  5.6438030941222897e+102
#define GSL_ROOT4_DBL_MAX  1.1579208923731620e+77
#define GSL_ROOT5_DBL_MAX  4.4765466227572707e+61
#define GSL_ROOT6_DBL_MAX  2.3756689782295612e+51
#define GSL_LOG_DBL_MAX    7.0978271289338397e+02

#define GSL_FLT_EPSILON        1.1920928955078125e-07
#define GSL_SQRT_FLT_EPSILON   3.4526698300124393e-04
#define GSL_ROOT3_FLT_EPSILON  4.9215666011518501e-03
#define GSL_ROOT4_FLT_EPSILON  1.8581361171917516e-02
#define GSL_ROOT5_FLT_EPSILON  4.1234622211652937e-02
#define GSL_ROOT6_FLT_EPSILON  7.0153878019335827e-02
#define GSL_LOG_FLT_EPSILON   (-1.5942385152878742e+01)

#define GSL_FLT_MIN        1.1754943508222875e-38
#define GSL_SQRT_FLT_MIN   1.0842021724855044e-19
#define GSL_ROOT3_FLT_MIN  2.2737367544323241e-13
#define GSL_ROOT4_FLT_MIN  3.2927225399135965e-10
#define GSL_ROOT5_FLT_MIN  2.5944428542140822e-08
#define GSL_ROOT6_FLT_MIN  4.7683715820312542e-07
#define GSL_LOG_FLT_MIN   (-8.7336544750553102e+01)

#define GSL_FLT_MAX        3.4028234663852886e+38
#define GSL_SQRT_FLT_MAX   1.8446743523953730e+19
#define GSL_ROOT3_FLT_MAX  6.9814635196223242e+12
#define GSL_ROOT4_FLT_MAX  4.2949672319999986e+09
#define GSL_ROOT5_FLT_MAX  5.0859007855960041e+07
#define GSL_ROOT6_FLT_MAX  2.6422459233807749e+06
#define GSL_LOG_FLT_MAX    8.8722839052068352e+01

#define GSL_SFLT_EPSILON        4.8828125000000000e-04
#define GSL_SQRT_SFLT_EPSILON   2.2097086912079612e-02
#define GSL_ROOT3_SFLT_EPSILON  7.8745065618429588e-02
#define GSL_ROOT4_SFLT_EPSILON  1.4865088937534013e-01
#define GSL_ROOT5_SFLT_EPSILON  2.1763764082403100e-01
#define GSL_ROOT6_SFLT_EPSILON  2.8061551207734325e-01
#define GSL_LOG_SFLT_EPSILON   (-7.6246189861593985e+00)

/* !MACHINE CONSTANTS! */


/* a little internal backwards compatibility */
#define GSL_MACH_EPS  GSL_DBL_EPSILON



/* Here are the constants related to or derived from
 * machine constants. These are not to be confused with
 * the constants that define various precision levels
 * for the precision/error system.
 *
 * This information is determined at configure time
 * and is platform dependent. Edit at your own risk.
 *
 * PLATFORM: WHIZ-O-MATIC
 * CONFIG-DATE: Thu Nov 19 19:27:18 MST 1998
 * CONFIG-HOST: nnn.lanl.gov
 */

/* machine precision constants */
/* #define GSL_MACH_EPS         1.0e-15 */
#define GSL_SQRT_MACH_EPS       3.2e-08
#define GSL_ROOT3_MACH_EPS      1.0e-05
#define GSL_ROOT4_MACH_EPS      0.000178
#define GSL_ROOT5_MACH_EPS      0.00100
#define GSL_ROOT6_MACH_EPS      0.00316
#define GSL_LOG_MACH_EPS       (-34.54)
// gsl_machine.h
///////////////////////////////////////////////////////////////////////////////////////
//#include "gsl_types.h"        // DONE
#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif
///////////////////////////////////////////////////////////////////////////////////////
//#include "gsl_precision.h"    // DONE
extern "C" {
/* A type for the precision indicator.
 * This is mainly for pedagogy.
 */
typedef  unsigned int  gsl_prec_t;

/* The number of precision types.
 * Remember that precision-mode
 * can index an array.
 */
#define _GSL_PREC_T_NUM 3

/* Arrays containing derived
 * precision constants for the
 * different precision levels.
 */
GSL_VAR const double gsl_prec_eps[];
GSL_VAR const double gsl_prec_sqrt_eps[];
GSL_VAR const double gsl_prec_root3_eps[];
GSL_VAR const double gsl_prec_root4_eps[];
GSL_VAR const double gsl_prec_root5_eps[];
GSL_VAR const double gsl_prec_root6_eps[];

}
///////////////////////////////////////////////////////////////////////////////////////
//#include "gsl_nan.h"          // DONE
#ifdef INFINITY
# define GSL_POSINF INFINITY
# define GSL_NEGINF (-INFINITY)
#elif defined(HUGE_VAL)
# define GSL_POSINF HUGE_VAL
# define GSL_NEGINF (-HUGE_VAL)
#else
# define GSL_POSINF (gsl_posinf())
# define GSL_NEGINF (gsl_neginf())
#endif

#ifdef NAN
# define GSL_NAN NAN
#elif defined(INFINITY)
# define GSL_NAN (INFINITY/INFINITY)
#else
# define GSL_NAN (gsl_nan())
#endif

#define GSL_POSZERO (+0)
#define GSL_NEGZERO (-0)

///////////////////////////////////////////////////////////////////////////////////////
//#include "gsl_pow_int.h"      // DONE
extern "C" {
#ifdef HAVE_INLINE
extern inline double gsl_pow_2(const double x);
extern inline double gsl_pow_3(const double x);
extern inline double gsl_pow_4(const double x);
extern inline double gsl_pow_5(const double x);
extern inline double gsl_pow_6(const double x);
extern inline double gsl_pow_7(const double x);
extern inline double gsl_pow_8(const double x);
extern inline double gsl_pow_9(const double x);

extern inline double gsl_pow_2(const double x) { return x*x;   }
extern inline double gsl_pow_3(const double x) { return x*x*x; }
extern inline double gsl_pow_4(const double x) { double x2 = x*x;   return x2*x2;    }
extern inline double gsl_pow_5(const double x) { double x2 = x*x;   return x2*x2*x;  }
extern inline double gsl_pow_6(const double x) { double x2 = x*x;   return x2*x2*x2; }
extern inline double gsl_pow_7(const double x) { double x3 = x*x*x; return x3*x3*x;  }
extern inline double gsl_pow_8(const double x) { double x2 = x*x;   double x4 = x2*x2; return x4*x4; }
extern inline double gsl_pow_9(const double x) { double x3 = x*x*x; return x3*x3*x3; }
#else
double gsl_pow_2(const double x);
double gsl_pow_3(const double x);
double gsl_pow_4(const double x);
double gsl_pow_5(const double x);
double gsl_pow_6(const double x);
double gsl_pow_7(const double x);
double gsl_pow_8(const double x);
double gsl_pow_9(const double x);
#endif

double gsl_pow_int(double x, int n);

} // gsl_pow_int.h
///////////////////////////////////////////////////////////////////////////////////////
//#include "gsl_math.h"         // DONE
#ifndef M_E
#define M_E        2.71828182845904523536028747135      /* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E    1.44269504088896340735992468100      /* log_2 (e) */
#endif

#ifndef M_LOG10E
#define M_LOG10E   0.43429448190325182765112891892      /* log_10 (e) */
#endif

#ifndef M_SQRT2
#define M_SQRT2    1.41421356237309504880168872421      /* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2  0.70710678118654752440084436210      /* sqrt(1/2) */
#endif


#ifndef M_SQRT3
#define M_SQRT3    1.73205080756887729352744634151      /* sqrt(3) */
#endif

#ifndef M_PI
#define M_PI       3.14159265358979323846264338328      /* pi */
#endif

#ifndef M_PI_2
#define M_PI_2     1.57079632679489661923132169164      /* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4     0.78539816339744830961566084582     /* pi/4 */
#endif

#ifndef M_SQRTPI
#define M_SQRTPI   1.77245385090551602729816748334      /* sqrt(pi) */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257389615890312      /* 2/sqrt(pi) */
#endif

#ifndef M_1_PI
#define M_1_PI     0.31830988618379067153776752675      /* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI     0.63661977236758134307553505349      /* 2/pi */
#endif

#ifndef M_LN10
#define M_LN10     2.30258509299404568401799145468      /* ln(10) */
#endif

#ifndef M_LN2
#define M_LN2      0.69314718055994530941723212146      /* ln(2) */
#endif

#ifndef M_LNPI
#define M_LNPI     1.14472988584940017414342735135      /* ln(pi) */
#endif

#ifndef M_EULER
#define M_EULER    0.57721566490153286060651209008      /* Euler constant */
#endif

extern "C" {
/* other needlessly compulsive abstractions */

#define GSL_IS_ODD(n)  ((n) & 1)
#define GSL_IS_EVEN(n) (!(GSL_IS_ODD(n)))
#define GSL_SIGN(x)    ((x) >= 0.0 ? 1 : -1)

/* Return nonzero if x is a real number, i.e. non NaN or infinite. */
#define GSL_IS_REAL(x) (gsl_finite(x))

/* Define MAX and MIN macros/functions if they don't exist. */

/* plain old macros for general use */
#define GSL_MAX(a,b) ((a) > (b) ? (a) : (b))
#define GSL_MIN(a,b) ((a) < (b) ? (a) : (b))

/* function versions of the above, in case they are needed */
double gsl_max (double a, double b);
double gsl_min (double a, double b);

/* inline-friendly strongly typed versions */
#ifdef HAVE_INLINE

extern inline int GSL_MAX_INT (int a, int b);
extern inline int GSL_MIN_INT (int a, int b);
extern inline double GSL_MAX_DBL (double a, double b);
extern inline double GSL_MIN_DBL (double a, double b);
extern inline long double GSL_MAX_LDBL (long double a, long double b);
extern inline long double GSL_MIN_LDBL (long double a, long double b);

extern inline int
GSL_MAX_INT (int a, int b)
{
  return GSL_MAX (a, b);
}

extern inline int
GSL_MIN_INT (int a, int b)
{
  return GSL_MIN (a, b);
}

extern inline double
GSL_MAX_DBL (double a, double b)
{
  return GSL_MAX (a, b);
}

extern inline double
GSL_MIN_DBL (double a, double b)
{
  return GSL_MIN (a, b);
}

extern inline long double
GSL_MAX_LDBL (long double a, long double b)
{
  return GSL_MAX (a, b);
}

extern inline long double
GSL_MIN_LDBL (long double a, long double b)
{
  return GSL_MIN (a, b);
}
#else
#define GSL_MAX_INT(a,b)   GSL_MAX(a,b)
#define GSL_MIN_INT(a,b)   GSL_MIN(a,b)
#define GSL_MAX_DBL(a,b)   GSL_MAX(a,b)
#define GSL_MIN_DBL(a,b)   GSL_MIN(a,b)
#define GSL_MAX_LDBL(a,b)  GSL_MAX(a,b)
#define GSL_MIN_LDBL(a,b)  GSL_MIN(a,b)
#endif /* HAVE_INLINE */

/* Definition of an arbitrary function with parameters */

struct gsl_function_struct 
{
  double (* function) (double x, void * params);
  void * params;
};

typedef struct gsl_function_struct gsl_function ;

#define GSL_FN_EVAL(F,x) (*((F)->function))(x,(F)->params)

/* Definition of an arbitrary function returning two values, r1, r2 */

struct gsl_function_fdf_struct 
{
  double (* f) (double x, void * params);
  double (* df) (double x, void * params);
  void (* fdf) (double x, void * params, double * f, double * df);
  void * params;
};

typedef struct gsl_function_fdf_struct gsl_function_fdf ;

#define GSL_FN_FDF_EVAL_F(FDF,x) (*((FDF)->f))(x,(FDF)->params)
#define GSL_FN_FDF_EVAL_DF(FDF,x) (*((FDF)->df))(x,(FDF)->params)
#define GSL_FN_FDF_EVAL_F_DF(FDF,x,y,dy) (*((FDF)->fdf))(x,(FDF)->params,(y),(dy))


/* Definition of an arbitrary vector-valued function with parameters */

struct gsl_function_vec_struct 
{
  int (* function) (double x, double y[], void * params);
  void * params;
};

typedef struct gsl_function_vec_struct gsl_function_vec ;

#define GSL_FN_VEC_EVAL(F,x,y) (*((F)->function))(x,y,(F)->params)

} // gsl_math.h
///////////////////////////////////////////////////////////////////////////////////////
//#include "gsl_complex.h"      // DONE
extern "C" {
/* two consecutive built-in types as a complex number */
typedef double *       gsl_complex_packed ;
typedef float *        gsl_complex_packed_float  ;
typedef long double *  gsl_complex_packed_long_double ;

typedef const double *       gsl_const_complex_packed ;
typedef const float *        gsl_const_complex_packed_float  ;
typedef const long double *  gsl_const_complex_packed_long_double ;


/* 2N consecutive built-in types as N complex numbers */
typedef double *       gsl_complex_packed_array ;
typedef float *        gsl_complex_packed_array_float  ;
typedef long double *  gsl_complex_packed_array_long_double ;

typedef const double *       gsl_const_complex_packed_array ;
typedef const float *        gsl_const_complex_packed_array_float  ;
typedef const long double *  gsl_const_complex_packed_array_long_double ;


/* Yes... this seems weird. Trust us. The point is just that
   sometimes you want to make it obvious that something is
   an output value. The fact that it lacks a 'const' may not
   be enough of a clue for people in some contexts.
 */
typedef double *       gsl_complex_packed_ptr ;
typedef float *        gsl_complex_packed_float_ptr  ;
typedef long double *  gsl_complex_packed_long_double_ptr ;

typedef const double *       gsl_const_complex_packed_ptr ;
typedef const float *        gsl_const_complex_packed_float_ptr  ;
typedef const long double *  gsl_const_complex_packed_long_double_ptr ;


typedef struct
  {
    long double dat[2];
  }
gsl_complex_long_double;

typedef struct
  {
    double dat[2];
  }
gsl_complex;

typedef struct
  {
    float dat[2];
  }
gsl_complex_float;

#define GSL_REAL(z)     ((z).dat[0])
#define GSL_IMAG(z)     ((z).dat[1])
#define GSL_COMPLEX_P(zp) ((zp)->dat)
#define GSL_COMPLEX_P_REAL(zp)  ((zp)->dat[0])
#define GSL_COMPLEX_P_IMAG(zp)  ((zp)->dat[1])
#define GSL_COMPLEX_EQ(z1,z2) (((z1).dat[0] == (z2).dat[0]) && ((z1).dat[1] == (z2).dat[1]))

#define GSL_SET_COMPLEX(zp,x,y) do {(zp)->dat[0]=(x); (zp)->dat[1]=(y);} while(0)
#define GSL_SET_REAL(zp,x) do {(zp)->dat[0]=(x);} while(0)
#define GSL_SET_IMAG(zp,y) do {(zp)->dat[1]=(y);} while(0)

#define GSL_SET_COMPLEX_PACKED(zp,n,x,y) do {*((zp)+2*(n))=(x); *((zp)+(2*(n)+1))=(y);} while(0)

} // gsl_complex.h
///////////////////////////////////////////////////////////////////////////////////////
//#include "gsl_complex_math.h" // DONE
extern "C" {
/* Complex numbers */

gsl_complex gsl_complex_rect (double x, double y);  /* r= real+i*imag */
gsl_complex gsl_complex_polar (double r, double theta); /* r= r e^(i theta) */

#ifdef HAVE_INLINE
extern inline gsl_complex
gsl_complex_rect (double x, double y)
{                               /* return z = x + i y */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, x, y);
  return z;
}
#endif

#define GSL_COMPLEX_ONE (gsl_complex_rect(1.0,0.0))
#define GSL_COMPLEX_ZERO (gsl_complex_rect(0.0,0.0))
#define GSL_COMPLEX_NEGONE (gsl_complex_rect(-1.0,0.0))

/* Properties of complex numbers */

double gsl_complex_arg (gsl_complex z); /* return arg(z), -pi< arg(z) <=+pi */
double gsl_complex_abs (gsl_complex z);   /* return |z|   */
double gsl_complex_abs2 (gsl_complex z);  /* return |z|^2 */
double gsl_complex_logabs (gsl_complex z); /* return log|z| */

/* Complex arithmetic operators */

gsl_complex gsl_complex_add (gsl_complex a, gsl_complex b);  /* r=a+b */
gsl_complex gsl_complex_sub (gsl_complex a, gsl_complex b);  /* r=a-b */
gsl_complex gsl_complex_mul (gsl_complex a, gsl_complex b);  /* r=a*b */
gsl_complex gsl_complex_div (gsl_complex a, gsl_complex b);  /* r=a/b */
                                                           
gsl_complex gsl_complex_add_real (gsl_complex a, double x);  /* r=a+x */
gsl_complex gsl_complex_sub_real (gsl_complex a, double x);  /* r=a-x */
gsl_complex gsl_complex_mul_real (gsl_complex a, double x);  /* r=a*x */
gsl_complex gsl_complex_div_real (gsl_complex a, double x);  /* r=a/x */

gsl_complex gsl_complex_add_imag (gsl_complex a, double y);  /* r=a+iy */
gsl_complex gsl_complex_sub_imag (gsl_complex a, double y);  /* r=a-iy */
gsl_complex gsl_complex_mul_imag (gsl_complex a, double y);  /* r=a*iy */
gsl_complex gsl_complex_div_imag (gsl_complex a, double y);  /* r=a/iy */

gsl_complex gsl_complex_conjugate (gsl_complex z);  /* r=conj(z) */
gsl_complex gsl_complex_inverse (gsl_complex a);    /* r=1/a */
gsl_complex gsl_complex_negative (gsl_complex a);    /* r=-a */

/* Elementary Complex Functions */

gsl_complex gsl_complex_sqrt (gsl_complex z);  /* r=sqrt(z) */
gsl_complex gsl_complex_sqrt_real (double x);  /* r=sqrt(x) (x<0 ok) */

gsl_complex gsl_complex_pow (gsl_complex a, gsl_complex b);  /* r=a^b */
gsl_complex gsl_complex_pow_real (gsl_complex a, double b);  /* r=a^b */

gsl_complex gsl_complex_exp (gsl_complex a);    /* r=exp(a) */
gsl_complex gsl_complex_log (gsl_complex a);    /* r=log(a) (base e) */
gsl_complex gsl_complex_log10 (gsl_complex a);  /* r=log10(a) (base 10) */
gsl_complex gsl_complex_log_b (gsl_complex a, gsl_complex b);   /* r=log_b(a) (base=b) */

/* Complex Trigonometric Functions */

gsl_complex gsl_complex_sin (gsl_complex a);  /* r=sin(a) */
gsl_complex gsl_complex_cos (gsl_complex a);  /* r=cos(a) */
gsl_complex gsl_complex_sec (gsl_complex a);  /* r=sec(a) */
gsl_complex gsl_complex_csc (gsl_complex a);  /* r=csc(a) */
gsl_complex gsl_complex_tan (gsl_complex a);  /* r=tan(a) */
gsl_complex gsl_complex_cot (gsl_complex a);  /* r=cot(a) */

/* Inverse Complex Trigonometric Functions */

gsl_complex gsl_complex_arcsin (gsl_complex a);  /* r=arcsin(a) */
gsl_complex gsl_complex_arcsin_real (double a);  /* r=arcsin(a) */
gsl_complex gsl_complex_arccos (gsl_complex a);  /* r=arccos(a) */
gsl_complex gsl_complex_arccos_real (double a);  /* r=arccos(a) */
gsl_complex gsl_complex_arcsec (gsl_complex a);  /* r=arcsec(a) */
gsl_complex gsl_complex_arcsec_real (double a);  /* r=arcsec(a) */
gsl_complex gsl_complex_arccsc (gsl_complex a);  /* r=arccsc(a) */
gsl_complex gsl_complex_arccsc_real (double a);  /* r=arccsc(a) */
gsl_complex gsl_complex_arctan (gsl_complex a);  /* r=arctan(a) */
gsl_complex gsl_complex_arccot (gsl_complex a);  /* r=arccot(a) */

/* Complex Hyperbolic Functions */

gsl_complex gsl_complex_sinh (gsl_complex a);  /* r=sinh(a) */
gsl_complex gsl_complex_cosh (gsl_complex a);  /* r=coshh(a) */
gsl_complex gsl_complex_sech (gsl_complex a);  /* r=sech(a) */
gsl_complex gsl_complex_csch (gsl_complex a);  /* r=csch(a) */
gsl_complex gsl_complex_tanh (gsl_complex a);  /* r=tanh(a) */
gsl_complex gsl_complex_coth (gsl_complex a);  /* r=coth(a) */

/* Inverse Complex Hyperbolic Functions */

gsl_complex gsl_complex_arcsinh (gsl_complex a);  /* r=arcsinh(a) */
gsl_complex gsl_complex_arccosh (gsl_complex a);  /* r=arccosh(a) */
gsl_complex gsl_complex_arccosh_real (double a);  /* r=arccosh(a) */
gsl_complex gsl_complex_arcsech (gsl_complex a);  /* r=arcsech(a) */
gsl_complex gsl_complex_arccsch (gsl_complex a);  /* r=arccsch(a) */
gsl_complex gsl_complex_arctanh (gsl_complex a);  /* r=arctanh(a) */
gsl_complex gsl_complex_arctanh_real (double a);  /* r=arctanh(a) */
gsl_complex gsl_complex_arccoth (gsl_complex a);  /* r=arccoth(a) */
} // gsl_complex_math.h
///////////////////////////////////////////////////////////////////////////////////////
//#include "gsl_errno.h"        // DONE
extern "C" {
enum { 
  GSL_SUCCESS  = 0, 
  GSL_FAILURE  = -1,
  GSL_CONTINUE = -2,  /* iteration has not converged */
  GSL_EDOM     = 1,   /* input domain error, e.g sqrt(-1) */
  GSL_ERANGE   = 2,   /* output range error, e.g. exp(1e100) */
  GSL_EFAULT   = 3,   /* invalid pointer */
  GSL_EINVAL   = 4,   /* invalid argument supplied by user */
  GSL_EFAILED  = 5,   /* generic failure */
  GSL_EFACTOR  = 6,   /* factorization failed */
  GSL_ESANITY  = 7,   /* sanity check failed - shouldn't happen */
  GSL_ENOMEM   = 8,   /* malloc failed */
  GSL_EBADFUNC = 9,   /* problem with user-supplied function */
  GSL_ERUNAWAY = 10,  /* iterative process is out of control */
  GSL_EMAXITER = 11,  /* exceeded max number of iterations */
  GSL_EZERODIV = 12,  /* tried to divide by zero */
  GSL_EBADTOL  = 13,  /* user specified an invalid tolerance */
  GSL_ETOL     = 14,  /* failed to reach the specified tolerance */
  GSL_EUNDRFLW = 15,  /* underflow */
  GSL_EOVRFLW  = 16,  /* overflow  */
  GSL_ELOSS    = 17,  /* loss of accuracy */
  GSL_EROUND   = 18,  /* failed because of roundoff error */
  GSL_EBADLEN  = 19,  /* matrix, vector lengths are not conformant */
  GSL_ENOTSQR  = 20,  /* matrix not square */
  GSL_ESING    = 21,  /* apparent singularity detected */
  GSL_EDIVERGE = 22,  /* integral or series is divergent */
  GSL_EUNSUP   = 23,  /* requested feature is not supported by the hardware */
  GSL_EUNIMPL  = 24,  /* requested feature not (yet) implemented */
  GSL_ECACHE   = 25,  /* cache limit exceeded */
  GSL_ETABLE   = 26,  /* table limit exceeded */
  GSL_ENOPROG  = 27,  /* iteration is not making progress towards solution */
  GSL_ENOPROGJ = 28,  /* jacobian evaluations are not improving the solution */
  GSL_ETOLF    = 29,  /* cannot reach the specified tolerance in F */
  GSL_ETOLX    = 30,  /* cannot reach the specified tolerance in X */
  GSL_ETOLG    = 31,  /* cannot reach the specified tolerance in gradient */
  GSL_EOF      = 32   /* end of file */
} ;

void gsl_error (const char * reason, const char * file, int line,
                int gsl_errno);

void gsl_stream_printf (const char *label, const char *file,
                        int line, const char *reason);

const char * gsl_strerror (const int gsl_errno);

typedef void gsl_error_handler_t (const char * reason, const char * file,
                                  int line, int gsl_errno);

typedef void gsl_stream_handler_t (const char * label, const char * file,
                                   int line, const char * reason);

gsl_error_handler_t * 
gsl_set_error_handler (gsl_error_handler_t * new_handler);

gsl_error_handler_t *
gsl_set_error_handler_off (void);

gsl_stream_handler_t * 
gsl_set_stream_handler (gsl_stream_handler_t * new_handler);

FILE * gsl_set_stream (FILE * new_stream);

/* GSL_ERROR: call the error handler, and return the error code */

#define GSL_ERROR(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return gsl_errno ; \
       } while (0)

/* GSL_ERROR_VAL: call the error handler, and return the given value */

#define GSL_ERROR_VAL(reason, gsl_errno, value) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return value ; \
       } while (0)

/* GSL_ERROR_VOID: call the error handler, and then return
   (for void functions which still need to generate an error) */

#define GSL_ERROR_VOID(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return ; \
       } while (0)

/* GSL_ERROR_NULL suitable for out-of-memory conditions */

#define GSL_ERROR_NULL(reason, gsl_errno) GSL_ERROR_VAL(reason, gsl_errno, 0)

/* Sometimes you have several status results returned from
 * function calls and you want to combine them in some sensible
 * way. You cannot produce a "total" status condition, but you can
 * pick one from a set of conditions based on an implied hierarchy.
 *
 * In other words:
 *    you have: status_a, status_b, ...
 *    you want: status = (status_a if it is bad, or status_b if it is bad,...)
 *
 * In this example you consider status_a to be more important and
 * it is checked first, followed by the others in the order specified.
 *
 * Here are some dumb macros to do this.
 */
#define GSL_ERROR_SELECT_2(a,b)       ((a) != GSL_SUCCESS ? (a) : ((b) != GSL_SUCCESS ? (b) : GSL_SUCCESS))
#define GSL_ERROR_SELECT_3(a,b,c)     ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_2(b,c))
#define GSL_ERROR_SELECT_4(a,b,c,d)   ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_3(b,c,d))
#define GSL_ERROR_SELECT_5(a,b,c,d,e) ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_4(b,c,d,e))

#define GSL_STATUS_UPDATE(sp, s) do { if ((s) != GSL_SUCCESS) *(sp) = (s);} while(0)
} // gsl_errno.h
///////////////////////////////////////////////////////////////////////////////////////
//#include "gsl_message.h"      // DONE
extern "C" {
/* Provide a general messaging service for client use.  Messages can
 * be selectively turned off at compile time by defining an
 * appropriate message mask. Client code which uses the GSL_MESSAGE()
 * macro must provide a mask which is or'ed with the GSL_MESSAGE_MASK.
 *
 * The messaging service can be completely turned off
 * by defining GSL_MESSAGING_OFF.  */

void gsl_message(const char * message, const char * file, int line,
                 unsigned int mask);

#ifndef GSL_MESSAGE_MASK
#define GSL_MESSAGE_MASK 0xffffffffu /* default all messages allowed */
#endif

GSL_VAR unsigned int gsl_message_mask ;

/* Provide some symolic masks for client ease of use. */

enum {
  GSL_MESSAGE_MASK_A = 1,
  GSL_MESSAGE_MASK_B = 2,
  GSL_MESSAGE_MASK_C = 4,
  GSL_MESSAGE_MASK_D = 8,
  GSL_MESSAGE_MASK_E = 16,
  GSL_MESSAGE_MASK_F = 32,
  GSL_MESSAGE_MASK_G = 64,
  GSL_MESSAGE_MASK_H = 128
} ;

#ifdef GSL_MESSAGING_OFF        /* throw away messages */ 
#define GSL_MESSAGE(message, mask) do { } while(0)
#else                           /* output all messages */
#define GSL_MESSAGE(message, mask) \
       do { \
       if (mask & GSL_MESSAGE_MASK) \
         gsl_message (message, __FILE__, __LINE__, mask) ; \
       } while (0)
#endif
} // gsl_message.h
///////////////////////////////////////////////////////////////////////////////////////
//#include "gsl_sf_result.h"    // DONE
extern "C" {
struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);
} // gsl_sf_result.h
///////////////////////////////////////////////////////////////////////////////////////
//#include "gsl_sf_erf.h"       // DONE
extern "C" {
/* Complementary Error Function
 * erfc(x) := 2/Sqrt[Pi] Integrate[Exp[-t^2], {t,x,Infinity}]
 *
 * exceptions: none
 */
int gsl_sf_erfc_e(double x, gsl_sf_result * result);
double gsl_sf_erfc(double x);


/* Log Complementary Error Function
 *
 * exceptions: none
 */
int gsl_sf_log_erfc_e(double x, gsl_sf_result * result);
double gsl_sf_log_erfc(double x);


/* Error Function
 * erf(x) := 2/Sqrt[Pi] Integrate[Exp[-t^2], {t,0,x}]
 *
 * exceptions: none
 */
int gsl_sf_erf_e(double x, gsl_sf_result * result);
double gsl_sf_erf(double x);


/* Probability functions:
 * Z(x) :  Abramowitz+Stegun 26.2.1
 * Q(x) :  Abramowitz+Stegun 26.2.3
 *
 * exceptions: none
 */
int gsl_sf_erf_Z_e(double x, gsl_sf_result * result);
int gsl_sf_erf_Q_e(double x, gsl_sf_result * result);
double gsl_sf_erf_Z(double x);
double gsl_sf_erf_Q(double x);


/* Hazard function, also known as the inverse Mill's ratio.
 *
 *   H(x) := Z(x)/Q(x)
 *         = Sqrt[2/Pi] Exp[-x^2 / 2] / Erfc[x/Sqrt[2]]
 *
 * exceptions: GSL_EUNDRFLW
 */
int gsl_sf_hazard_e(double x, gsl_sf_result * result);
double gsl_sf_hazard(double x);
} // gsl_sf_erf.h
///////////////////////////////////////////////////////////////////////////////////////
//#include "gsl_sf_exp.h"       // DONE
extern "C" {
/* Provide an exp() function with GSL semantics,
 * i.e. with proper error checking, etc.
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_e(const double x, gsl_sf_result * result);
double gsl_sf_exp(const double x);


/* Exp(x)
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_e10_e(const double x, gsl_sf_result_e10 * result);


/* Exponentiate and multiply by a given factor:  y * Exp(x)
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_mult_e(const double x, const double y, gsl_sf_result * result);
double gsl_sf_exp_mult(const double x, const double y);


/* Exponentiate and multiply by a given factor:  y * Exp(x)
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_mult_e10_e(const double x, const double y, gsl_sf_result_e10 * result);


/* exp(x)-1
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_expm1_e(const double x, gsl_sf_result * result);
double gsl_sf_expm1(const double x);


/* (exp(x)-1)/x = 1 + x/2 + x^2/(2*3) + x^3/(2*3*4) + ...
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_exprel_e(const double x, gsl_sf_result * result);
double gsl_sf_exprel(const double x);


/* 2(exp(x)-1-x)/x^2 = 1 + x/3 + x^2/(3*4) + x^3/(3*4*5) + ...
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_exprel_2_e(double x, gsl_sf_result * result);
double gsl_sf_exprel_2(const double x);


/* Similarly for the N-th generalization of
 * the above. The so-called N-relative exponential
 *
 * exprel_N(x) = N!/x^N (exp(x) - Sum[x^k/k!, {k,0,N-1}])
 *             = 1 + x/(N+1) + x^2/((N+1)(N+2)) + ...
 *             = 1F1(1,1+N,x)
 */
int gsl_sf_exprel_n_e(const int n, const double x, gsl_sf_result * result);
double gsl_sf_exprel_n(const int n, const double x);


/* Exponentiate a quantity with an associated error.
 */
int gsl_sf_exp_err_e(const double x, const double dx, gsl_sf_result * result);

/* Exponentiate a quantity with an associated error.
 */
int gsl_sf_exp_err_e10_e(const double x, const double dx, gsl_sf_result_e10 * result);


/* Exponentiate and multiply by a given factor:  y * Exp(x),
 * for quantities with associated errors.
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_mult_err_e(const double x, const double dx, const double y, const double dy, gsl_sf_result * result);


/* Exponentiate and multiply by a given factor:  y * Exp(x),
 * for quantities with associated errors.
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_mult_err_e10_e(const double x, const double dx, const double y, const double dy, gsl_sf_result_e10 * result);
} // gsl_sf_exp.h
///////////////////////////////////////////////////////////////////////////////////////
//#include "gsl_sf_log.h"       // DONE
extern "C" {
/* Provide a logarithm function with GSL semantics.
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_log_e(const double x, gsl_sf_result * result);
double gsl_sf_log(const double x);


/* Log(|x|)
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_log_abs_e(const double x, gsl_sf_result * result);
double gsl_sf_log_abs(const double x);


/* Complex Logarithm
 *   exp(lnr + I theta) = zr + I zi
 * Returns argument in [-pi,pi].
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_complex_log_e(const double zr, const double zi, gsl_sf_result * lnr, gsl_sf_result * theta);


/* Log(1 + x)
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_log_1plusx_e(const double x, gsl_sf_result * result);
double gsl_sf_log_1plusx(const double x);


/* Log(1 + x) - x
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_log_1plusx_mx_e(const double x, gsl_sf_result * result);
double gsl_sf_log_1plusx_mx(const double x);
} // gsl_sf_log.h
///////////////////////////////////////////////////////////////////////////////////////
//#include "gsl_sf_psi.h"       // DONE
extern "C" {
/* Poly-Gamma Functions
 *
 * psi(m,x) := (d/dx)^m psi(0,x) = (d/dx)^{m+1} log(gamma(x))
 */


/* Di-Gamma Function  psi(n) = psi(0,n)
 *
 * n > 0
 * exceptions: GSL_EDOM
 */
int     gsl_sf_psi_int_e(const int n, gsl_sf_result * result);
double  gsl_sf_psi_int(const int n);


/* Di-Gamma Function psi(x) = psi(0, x)
 *
 * x != 0.0, -1.0, -2.0, ...
 * exceptions: GSL_EDOM, GSL_ELOSS
 */
int     gsl_sf_psi_e(const double x, gsl_sf_result * result);
double  gsl_sf_psi(const double x);


/* Di-Gamma Function Re[psi(1 + I y)]
 *
 * exceptions: none
 */
int     gsl_sf_psi_1piy_e(const double y, gsl_sf_result * result);
double  gsl_sf_psi_1piy(const double y);


/* Di-Gamma Function psi(z) for general complex argument z = x + iy
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_complex_psi_e(
  const double x,
  const double y,
  gsl_sf_result * result_re,
  gsl_sf_result * result_im
  );


/* Tri-Gamma Function psi^(1)(n)
 *
 * n > 0
 * exceptions: GSL_EDOM
 */
int     gsl_sf_psi_1_int_e(const int n, gsl_sf_result * result);
double  gsl_sf_psi_1_int(const int n);


/* Tri-Gamma Function psi^(1)(x)
 *
 * x != 0.0, -1.0, -2.0, ...
 * exceptions: GSL_EDOM, GSL_ELOSS
 */
int     gsl_sf_psi_1_e(const double x, gsl_sf_result * result);
double  gsl_sf_psi_1(const double x);


/* Poly-Gamma Function psi^(n)(x)
 *
 * n >= 0, x > 0.0
 * exceptions: GSL_EDOM
 */
int     gsl_sf_psi_n_e(const int n, const double x, gsl_sf_result * result);
double  gsl_sf_psi_n(const int n, const double x);
} // gsl_sf_psi.h
///////////////////////////////////////////////////////////////////////////////////////
//#include "gsl_sf_trig.h"      // DONE
extern "C" {
/* Sin(x) with GSL semantics. This is actually important
 * because we want to control the error estimate, and trying
 * to guess the error for the standard library implementation
 * every time it is used would be a little goofy.
 */
int gsl_sf_sin_e(double x, gsl_sf_result * result);
double gsl_sf_sin(const double x);


/* Cos(x) with GSL semantics.
 */
int gsl_sf_cos_e(double x, gsl_sf_result * result);
double gsl_sf_cos(const double x);


/* Hypot(x,y) with GSL semantics.
 */
int gsl_sf_hypot_e(const double x, const double y, gsl_sf_result * result);
double gsl_sf_hypot(const double x, const double y);


/* Sin(z) for complex z
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_complex_sin_e(const double zr, const double zi, gsl_sf_result * szr, gsl_sf_result * szi);


/* Cos(z) for complex z
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_complex_cos_e(const double zr, const double zi, gsl_sf_result * czr, gsl_sf_result * czi);


/* Log(Sin(z)) for complex z
 *
 * exceptions: GSL_EDOM, GSL_ELOSS
 */
int gsl_sf_complex_logsin_e(const double zr, const double zi, gsl_sf_result * lszr, gsl_sf_result * lszi);


/* Sinc(x) = sin(pi x) / (pi x)
 *
 * exceptions: none
 */
int gsl_sf_sinc_e(double x, gsl_sf_result * result);
double gsl_sf_sinc(const double x);


/* Log(Sinh(x)), x > 0
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_lnsinh_e(const double x, gsl_sf_result * result);
double gsl_sf_lnsinh(const double x);


/* Log(Cosh(x))
 *
 * exceptions: none
 */
int gsl_sf_lncosh_e(const double x, gsl_sf_result * result);
double gsl_sf_lncosh(const double x);


/* Convert polar to rectlinear coordinates.
 *
 * exceptions: GSL_ELOSS
 */
int gsl_sf_polar_to_rect(const double r, const double theta, gsl_sf_result * x, gsl_sf_result * y);

/* Convert rectilinear to polar coordinates.
 * return argument in range [-pi, pi]
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_rect_to_polar(const double x, const double y, gsl_sf_result * r, gsl_sf_result * theta);

/* Sin(x) for quantity with an associated error.
 */
int gsl_sf_sin_err_e(const double x, const double dx, gsl_sf_result * result);


/* Cos(x) for quantity with an associated error.
 */
int gsl_sf_cos_err_e(const double x, const double dx, gsl_sf_result * result);


/* Force an angle to lie in the range (-pi,pi].
 *
 * exceptions: GSL_ELOSS
 */
int gsl_sf_angle_restrict_symm_e(double * theta);
double gsl_sf_angle_restrict_symm(const double theta);


/* Force an angle to lie in the range [0, 2pi)
 *
 * exceptions: GSL_ELOSS
 */
int gsl_sf_angle_restrict_pos_e(double * theta);
double gsl_sf_angle_restrict_pos(const double theta);


int gsl_sf_angle_restrict_symm_err_e(const double theta, gsl_sf_result * result);

int gsl_sf_angle_restrict_pos_err_e(const double theta, gsl_sf_result * result);
} // gsl_sf_trig.h
///////////////////////////////////////////////////////////////////////////////////////
//#include "gsl_sf_gamma.h"     // DONE
extern "C" {
/* Log[Gamma(x)], x not a negative integer
 * Uses real Lanczos method.
 * Returns the real part of Log[Gamma[x]] when x < 0,
 * i.e. Log[|Gamma[x]|].
 *
 * exceptions: GSL_EDOM, GSL_EROUND
 */
int gsl_sf_lngamma_e(double x, gsl_sf_result * result);
double gsl_sf_lngamma(const double x);


/* Log[Gamma(x)], x not a negative integer
 * Uses real Lanczos method. Determines
 * the sign of Gamma[x] as well as Log[|Gamma[x]|] for x < 0.
 * So Gamma[x] = sgn * Exp[result_lg].
 *
 * exceptions: GSL_EDOM, GSL_EROUND
 */
int gsl_sf_lngamma_sgn_e(double x, gsl_sf_result * result_lg, double *sgn);


/* Gamma(x), x not a negative integer
 * Uses real Lanczos method.
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EROUND
 */
int gsl_sf_gamma_e(const double x, gsl_sf_result * result);
double gsl_sf_gamma(const double x);


/* Regulated Gamma Function, x > 0
 * Gamma^*(x) = Gamma(x)/(Sqrt[2Pi] x^(x-1/2) exp(-x))
 *            = (1 + 1/(12x) + ...),  x->Inf
 * A useful suggestion of Temme.
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gammastar_e(const double x, gsl_sf_result * result);
double gsl_sf_gammastar(const double x);


/* 1/Gamma(x)
 * Uses real Lanczos method.
 *
 * exceptions: GSL_EUNDRFLW, GSL_EROUND
 */
int gsl_sf_gammainv_e(const double x, gsl_sf_result * result);
double gsl_sf_gammainv(const double x);


/* Log[Gamma(z)] for z complex, z not a negative integer
 * Uses complex Lanczos method. Note that the phase part (arg)
 * is not well-determined when |z| is very large, due
 * to inevitable roundoff in restricting to (-Pi,Pi].
 * This will raise the GSL_ELOSS exception when it occurs.
 * The absolute value part (lnr), however, never suffers.
 *
 * Calculates:
 *   lnr = log|Gamma(z)|
 *   arg = arg(Gamma(z))  in (-Pi, Pi]
 *
 * exceptions: GSL_EDOM, GSL_ELOSS
 */
int gsl_sf_lngamma_complex_e(double zr, double zi, gsl_sf_result * lnr, gsl_sf_result * arg);


/* x^n / n!
 *
 * x >= 0.0, n >= 0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_taylorcoeff_e(const int n, const double x, gsl_sf_result * result);
double gsl_sf_taylorcoeff(const int n, const double x);


/* n!
 *
 * exceptions: GSL_EDOM, GSL_OVRFLW
 */
int gsl_sf_fact_e(const unsigned int n, gsl_sf_result * result);
double gsl_sf_fact(const unsigned int n);


/* n!! = n(n-2)(n-4) ... 
 *
 * exceptions: GSL_EDOM, GSL_OVRFLW
 */
int gsl_sf_doublefact_e(const unsigned int n, gsl_sf_result * result);
double gsl_sf_doublefact(const unsigned int n);


/* log(n!) 
 * Faster than ln(Gamma(n+1)) for n < 170; defers for larger n.
 *
 * exceptions: none
 */
int gsl_sf_lnfact_e(const unsigned int n, gsl_sf_result * result);
double gsl_sf_lnfact(const unsigned int n);


/* log(n!!) 
 *
 * exceptions: none
 */
int gsl_sf_lndoublefact_e(const unsigned int n, gsl_sf_result * result);
double gsl_sf_lndoublefact(const unsigned int n);


/* log(n choose m)
 *
 * exceptions: GSL_EDOM 
 */
int gsl_sf_lnchoose_e(unsigned int n, unsigned int m, gsl_sf_result * result);
double gsl_sf_lnchoose(unsigned int n, unsigned int m);


/* n choose m
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_choose_e(unsigned int n, unsigned int m, gsl_sf_result * result);
double gsl_sf_choose(unsigned int n, unsigned int m);


/* Logarithm of Pochhammer (Apell) symbol
 *   log( (a)_x )
 *   where (a)_x := Gamma[a + x]/Gamma[a]
 *
 * a > 0, a+x > 0
 *
 * exceptions:  GSL_EDOM
 */
int gsl_sf_lnpoch_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_lnpoch(const double a, const double x);


/* Logarithm of Pochhammer (Apell) symbol, with sign information.
 *   result = log( |(a)_x| )
 *   sgn    = sgn( (a)_x )
 *   where (a)_x := Gamma[a + x]/Gamma[a]
 *
 * a != neg integer, a+x != neg integer
 *
 * exceptions:  GSL_EDOM
 */
int gsl_sf_lnpoch_sgn_e(const double a, const double x, gsl_sf_result * result, double * sgn);


/* Pochhammer (Apell) symbol
 *   (a)_x := Gamma[a + x]/Gamma[x]
 *
 * a != neg integer, a+x != neg integer
 *
 * exceptions:  GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_poch_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_poch(const double a, const double x);


/* Relative Pochhammer (Apell) symbol
 *   ((a,x) - 1)/x
 *   where (a,x) = (a)_x := Gamma[a + x]/Gamma[a]
 *
 * exceptions:  GSL_EDOM
 */
int gsl_sf_pochrel_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_pochrel(const double a, const double x);


/* Normalized Incomplete Gamma Function
 *
 * Q(a,x) = 1/Gamma(a) Integral[ t^(a-1) e^(-t), {t,x,Infinity} ]
 *
 * a >= 0, x >= 0
 *   Q(a,0) := 1
 *   Q(0,x) := 0, x != 0
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gamma_inc_Q_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_gamma_inc_Q(const double a, const double x);


/* Complementary Normalized Incomplete Gamma Function
 *
 * P(a,x) = 1/Gamma(a) Integral[ t^(a-1) e^(-t), {t,0,x} ]
 *
 * a > 0, x >= 0
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gamma_inc_P_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_gamma_inc_P(const double a, const double x);


/* Non-normalized Incomplete Gamma Function
 *
 * Gamma(a,x) := Integral[ t^(a-1) e^(-t), {t,x,Infinity} ]
 *
 * x >= 0.0
 *   Gamma(a, 0) := Gamma(a)
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gamma_inc_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_gamma_inc(const double a, const double x);


/* Logarithm of Beta Function
 * Log[B(a,b)]
 *
 * a > 0, b > 0
 * exceptions: GSL_EDOM
 */
int gsl_sf_lnbeta_e(const double a, const double b, gsl_sf_result * result);
double gsl_sf_lnbeta(const double a, const double b);

int gsl_sf_lnbeta_sgn_e(const double x, const double y, gsl_sf_result * result, double * sgn);


/* Beta Function
 * B(a,b)
 *
 * a > 0, b > 0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_beta_e(const double a, const double b, gsl_sf_result * result);
double gsl_sf_beta(const double a, const double b);


/* Normalized Incomplete Beta Function
 * B_x(a,b)/B(a,b)
 *
 * a > 0, b > 0, 0 <= x <= 1
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int gsl_sf_beta_inc_e(const double a, const double b, const double x, gsl_sf_result * result);
double gsl_sf_beta_inc(const double a, const double b, const double x);


/* The maximum x such that gamma(x) is not
 * considered an overflow.
 */
#define GSL_SF_GAMMA_XMAX  171.0

/* The maximum n such that gsl_sf_fact(n) does not give an overflow. */
#define GSL_SF_FACT_NMAX 170

/* The maximum n such that gsl_sf_doublefact(n) does not give an overflow. */
#define GSL_SF_DOUBLEFACT_NMAX 297
} // gsl_sf_gamma.h
///////////////////////////////////////////////////////////////////////////////////////
//#include "gsl_sf_expint.h"    // DONE
extern "C" {
/* E_1(x) := Re[ Integrate[ Exp[-xt]/t, {t,1,Infinity}] ]
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_expint_E1_e(const double x, gsl_sf_result * result);
double  gsl_sf_expint_E1(const double x);


/* E_2(x) := Re[ Integrate[ Exp[-xt]/t^2, {t,1,Infinity}] ]
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_expint_E2_e(const double x, gsl_sf_result * result);
double  gsl_sf_expint_E2(const double x);


/* E_1_scaled(x) := exp(x) E_1(x)
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_expint_E1_scaled_e(const double x, gsl_sf_result * result);
double  gsl_sf_expint_E1_scaled(const double x);


/* E_2_scaled(x) := exp(x) E_2(x)
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_expint_E2_scaled_e(const double x, gsl_sf_result * result);
double  gsl_sf_expint_E2_scaled(const double x);


/* Ei(x) := - PV Integrate[ Exp[-t]/t, {t,-x,Infinity}]
 *       :=   PV Integrate[ Exp[t]/t, {t,-Infinity,x}]
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_expint_Ei_e(const double x, gsl_sf_result * result);
double  gsl_sf_expint_Ei(const double x);


/* Ei_scaled(x) := exp(-x) Ei(x)
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_expint_Ei_scaled_e(const double x, gsl_sf_result * result);
double  gsl_sf_expint_Ei_scaled(const double x);


/* Shi(x) := Integrate[ Sinh[t]/t, {t,0,x}]
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_Shi_e(const double x, gsl_sf_result * result);
double  gsl_sf_Shi(const double x);


/* Chi(x) := Re[ M_EULER + log(x) + Integrate[(Cosh[t]-1)/t, {t,0,x}] ]
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_Chi_e(const double x, gsl_sf_result * result);
double  gsl_sf_Chi(const double x);


/* Ei_3(x) := Integral[ Exp[-t^3], {t,0,x}]
 *
 * x >= 0.0
 * exceptions: GSL_EDOM
 */
int     gsl_sf_expint_3_e(const double x, gsl_sf_result * result);
double  gsl_sf_expint_3(double x);


/* Si(x) := Integrate[ Sin[t]/t, {t,0,x}]
 *
 * exceptions: none
 */
int     gsl_sf_Si_e(const double x, gsl_sf_result * result);
double  gsl_sf_Si(const double x);


/* Ci(x) := -Integrate[ Cos[t]/t, {t,x,Infinity}]
 *
 * x > 0.0
 * exceptions: GSL_EDOM 
 */
int     gsl_sf_Ci_e(const double x, gsl_sf_result * result);
double  gsl_sf_Ci(const double x);


/* AtanInt(x) := Integral[ Arctan[t]/t, {t,0,x}]
 *
 *
 * exceptions:
 */
int     gsl_sf_atanint_e(const double x, gsl_sf_result * result);
double  gsl_sf_atanint(const double x);
} // gsl_sf_expint.h
///////////////////////////////////////////////////////////////////////////////////////
//#include "gsl_sf_pow_int.h"   // DONE
extern "C" {
/* Calculate x^n.
 * Does not check for overflow/underflow.
 */
int     gsl_sf_pow_int_e(double x, int n, gsl_sf_result * result);
double  gsl_sf_pow_int(const double x, const int n);
} // gsl_sf_pow_int.h
///////////////////////////////////////////////////////////////////////////////////////
//#include "gsl_sf_zeta.h"      // DONE
extern "C" {
/* Riemann Zeta Function
 * zeta(n) = Sum[ k^(-n), {k,1,Infinity} ]
 *
 * n=integer, n != 1
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_zeta_int_e(const int n, gsl_sf_result * result);
double gsl_sf_zeta_int(const int n);


/* Riemann Zeta Function
 * zeta(x) = Sum[ k^(-s), {k,1,Infinity} ], s != 1.0
 *
 * s != 1.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_zeta_e(const double s, gsl_sf_result * result);
double gsl_sf_zeta(const double s);


/* Riemann Zeta Function minus 1
 *   useful for evaluating the fractional part
 *   of Riemann zeta for large argument
 *
 * s != 1.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_zetam1_e(const double s, gsl_sf_result * result);
double gsl_sf_zetam1(const double s);


/* Riemann Zeta Function minus 1 for integer arg
 *   useful for evaluating the fractional part
 *   of Riemann zeta for large argument
 *
 * s != 1.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_zetam1_int_e(const int s, gsl_sf_result * result);
double gsl_sf_zetam1_int(const int s);


/* Hurwitz Zeta Function
 * zeta(s,q) = Sum[ (k+q)^(-s), {k,0,Infinity} ]
 *
 * s > 1.0, q > 0.0
 * exceptions: GSL_EDOM, GSL_EUNDRFLW, GSL_EOVRFLW
 */
int gsl_sf_hzeta_e(const double s, const double q, gsl_sf_result * result);
double gsl_sf_hzeta(const double s, const double q);


/* Eta Function
 * eta(n) = (1-2^(1-n)) zeta(n)
 *
 * exceptions: GSL_EUNDRFLW, GSL_EOVRFLW
 */
int gsl_sf_eta_int_e(int n, gsl_sf_result * result);
double gsl_sf_eta_int(const int n);


/* Eta Function
 * eta(s) = (1-2^(1-s)) zeta(s)
 *
 * exceptions: GSL_EUNDRFLW, GSL_EOVRFLW
 */
int gsl_sf_eta_e(const double s, gsl_sf_result * result);
double gsl_sf_eta(const double s);
} // gsl_sf_zeta.h
///////////////////////////////////////////////////////////////////////////////////////

//#include "error.h"            // DONE
#define OVERFLOW_ERROR(result) do { (result)->val = GSL_POSINF; (result)->err = GSL_POSINF; GSL_ERROR ("overflow", GSL_EOVRFLW); } while(0)

#define UNDERFLOW_ERROR(result) do { (result)->val = 0.0; (result)->err = GSL_DBL_MIN; GSL_ERROR ("underflow", GSL_EUNDRFLW); } while(0)

#define INTERNAL_OVERFLOW_ERROR(result) do { (result)->val = GSL_POSINF; (result)->err = GSL_POSINF; return GSL_EOVRFLW; } while(0)

#define INTERNAL_UNDERFLOW_ERROR(result) do { (result)->val = 0.0; (result)->err = GSL_DBL_MIN; return GSL_EUNDRFLW; } while(0)

#define DOMAIN_ERROR(result) do { (result)->val = GSL_NAN; (result)->err = GSL_NAN; GSL_ERROR ("domain error", GSL_EDOM); } while(0)

#define DOMAIN_ERROR_MSG(msg, result) do { (result)->val = GSL_NAN; (result)->err = GSL_NAN; GSL_ERROR ((msg), GSL_EDOM); } while(0)

#define DOMAIN_ERROR_E10(result) do { (result)->val = GSL_NAN; (result)->err = GSL_NAN; (result)->e10 = 0 ; GSL_ERROR ("domain error", GSL_EDOM); } while(0)

#define OVERFLOW_ERROR_E10(result) do { (result)->val = GSL_POSINF; (result)->err = GSL_POSINF; (result)->e10 = 0; GSL_ERROR ("overflow", GSL_EOVRFLW); } while(0)
#define UNDERFLOW_ERROR_E10(result) do { (result)->val = 0.0; (result)->err = GSL_DBL_MIN; (result)->e10 = 0; GSL_ERROR ("underflow", GSL_EUNDRFLW); } while(0)


#define OVERFLOW_ERROR_2(r1,r2) do { (r1)->val = GSL_POSINF; (r1)->err = GSL_POSINF; (r2)->val = GSL_POSINF ; (r2)->err=GSL_POSINF; GSL_ERROR ("overflow", GSL_EOVRFLW); } while(0)

#define UNDERFLOW_ERROR_2(r1,r2) do { (r1)->val = 0.0; (r1)->err = GSL_DBL_MIN; (r2)->val = 0.0 ; (r2)->err = GSL_DBL_MIN; GSL_ERROR ("underflow", GSL_EUNDRFLW); } while(0)

#define DOMAIN_ERROR_2(r1,r2) do { (r1)->val = GSL_NAN; (r1)->err = GSL_NAN;  (r2)->val = GSL_NAN; (r2)->err = GSL_NAN;  GSL_ERROR ("domain error", GSL_EDOM); } while(0)
// \error.h
///////////////////////////////////////////////////////////////////////////////////////
//#include "eval.h"             // DONE
#define EVAL_RESULT(fn) \
   gsl_sf_result result; \
   int status = fn; \
   if (status != GSL_SUCCESS) { \
     GSL_ERROR_VAL(#fn, status, result.val); \
   } ; \
   return result.val;

#define EVAL_DOUBLE(fn) \
   int status = fn; \
   if (status != GSL_SUCCESS) { \
     GSL_ERROR_VAL(#fn, status, result); \
   } ; \
   return result;
///////////////////////////////////////////////////////////////////////////////////////
//#include "check.h"            // DONE
#define CHECK_UNDERFLOW(r) if (fabs((r)->val) < GSL_DBL_MIN) GSL_ERROR("underflow", GSL_EUNDRFLW);
///////////////////////////////////////////////////////////////////////////////////////
//#include "gsl_sf_elementary.h" // DONE
extern "C" {
/* Multiplication.
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_multiply_e(const double x, const double y, gsl_sf_result * result);
double gsl_sf_multiply(const double x, const double y);


/* Multiplication of quantities with associated errors.
 */
int gsl_sf_multiply_err_e(const double x, const double dx, const double y, const double dy, gsl_sf_result * result);
} // gsl_sf_elementary.h
///////////////////////////////////////////////////////////////////////////////////////
//#include "chebyshev.h"        // DONE
/* data for a Chebyshev series over a given interval */

struct cheb_series_struct {
  double * c;   /* coefficients                */
  int order;    /* order of expansion          */
  double a;     /* lower interval point        */
  double b;     /* upper interval point        */
  int order_sp; /* effective single precision order */
};
typedef struct cheb_series_struct cheb_series;
///////////////////////////////////////////////////////////////////////////////////////
//#include "cheb_eval.c.in"     // DONE
extern "C" {
static inline int cheb_eval_e(const cheb_series * cs, const double x, gsl_sf_result * result) {
  int j;
  double d  = 0.0;
  double dd = 0.0;

  double y  = (2.0*x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2.0 * y;

  double e = 0.0;

  for(j = cs->order; j>=1; j--) {
    double temp = d;
    d = y2*d - dd + cs->c[j];
    e += fabs(y2*temp) + fabs(dd) + fabs(cs->c[j]);
    dd = temp;
  }

  { 
    double temp = d;
    d = y*d - dd + 0.5 * cs->c[0];
    e += fabs(y*temp) + fabs(dd) + 0.5 * fabs(cs->c[0]);
  }

  result->val = d;
  result->err = GSL_DBL_EPSILON * e + fabs(cs->c[cs->order]);

  return GSL_SUCCESS;
}
} // cheb_eval.c.in
///////////////////////////////////////////////////////////////////////////////////////
