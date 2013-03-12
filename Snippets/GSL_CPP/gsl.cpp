/* gsl.cpp
 * This file is intended to be the unified C++ file for the subset of GSL files I need for Rico.
 */

/* Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * Author:  G. Jungman 
 */
/* Basic complex arithmetic functions

 * Original version by Jorma Olavi Tähtinen <jotahtin@cc.hut.fi>
 *
 * Modified for GSL by Brian Gough, 3/2000
 */

/* The following references describe the methods used in these
 * functions,
 *
 *   T. E. Hull and Thomas F. Fairgrieve and Ping Tak Peter Tang,
 *   "Implementing Complex Elementary Functions Using Exception
 *   Handling", ACM Transactions on Mathematical Software, Volume 20
 *   (1994), pp 215-244, Corrigenda, p553
 *
 *   Hull et al, "Implementing the complex arcsin and arccosine
 *   functions using exception handling", ACM Transactions on
 *   Mathematical Software, Volume 23 (1997) pp 299-335
 *
 *   Abramowitz and Stegun, Handbook of Mathematical Functions, "Inverse
 *   Circular Functions in Terms of Real and Imaginary Parts", Formulas
 *   4.4.37, 4.4.38, 4.4.39
 */

#include "gsl.h"

// beta.c ///////////////////////////////////////////////////////////////////////
extern "C" { 

  static double isnegint (const double x) {
    return (x < 0) && (x == floor(x));
  }

  int gsl_sf_lnbeta_e(const double x, const double y, gsl_sf_result * result) {
    double sgn;
    int status = gsl_sf_lnbeta_sgn_e(x,y,result,&sgn);
    if (sgn == -1) {
      DOMAIN_ERROR(result);
    }
    return status;
  }
  
  int gsl_sf_lnbeta_sgn_e(const double x, const double y, gsl_sf_result * result, double * sgn) {
    /* CHECK_POINTER(result) */
    
    if(x == 0.0 || y == 0.0) {
      *sgn = 0.0;
      DOMAIN_ERROR(result);
    } else if (isnegint(x) || isnegint(y)) {
      *sgn = 0.0;
      DOMAIN_ERROR(result); /* not defined for negative integers */
    }
    
    /* See if we can handle the postive case with min/max < 0.2 */
    
    if (x > 0 && y > 0) {
      const double max = GSL_MAX(x,y);
      const double min = GSL_MIN(x,y);
      const double rat = min/max;
      
      if(rat < 0.2) {
	/* min << max, so be careful
	 * with the subtraction
	 */
	double lnpre_val;
	double lnpre_err;
	double lnpow_val;
	double lnpow_err;
	double t1, t2, t3;
	gsl_sf_result lnopr;
	gsl_sf_result gsx, gsy, gsxy;
	gsl_sf_gammastar_e(x, &gsx);
	gsl_sf_gammastar_e(y, &gsy);
	gsl_sf_gammastar_e(x+y, &gsxy);
	gsl_sf_log_1plusx_e(rat, &lnopr);
	lnpre_val = log(gsx.val*gsy.val/gsxy.val * M_SQRT2*M_SQRTPI);
	lnpre_err = gsx.err/gsx.val + gsy.err/gsy.val + gsxy.err/gsxy.val;
	t1 = min*log(rat);
	t2 = 0.5*log(min);
	t3 = (x+y-0.5)*lnopr.val;
	lnpow_val  = t1 - t2 - t3;
	lnpow_err  = GSL_DBL_EPSILON * (fabs(t1) + fabs(t2) + fabs(t3));
	lnpow_err += fabs(x+y-0.5) * lnopr.err;
	result->val  = lnpre_val + lnpow_val;
	result->err  = lnpre_err + lnpow_err;
	result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
	*sgn = 1.0;
	return GSL_SUCCESS;
      }
    }
    
    /* General case - Fallback */
    {
      gsl_sf_result lgx, lgy, lgxy;
      double sgx, sgy, sgxy, xy = x+y;
      int stat_gx  = gsl_sf_lngamma_sgn_e(x, &lgx, &sgx);
      int stat_gy  = gsl_sf_lngamma_sgn_e(y, &lgy, &sgy);
      int stat_gxy = gsl_sf_lngamma_sgn_e(xy, &lgxy, &sgxy);
      *sgn = sgx * sgy * sgxy;
      result->val  = lgx.val + lgy.val - lgxy.val;
      result->err  = lgx.err + lgy.err + lgxy.err;
      result->err += GSL_DBL_EPSILON * (fabs(lgx.val) + fabs(lgy.val) + fabs(lgxy.val));
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_ERROR_SELECT_3(stat_gx, stat_gy, stat_gxy);
    }
  }
  
  
  int gsl_sf_beta_e(const double x, const double y, gsl_sf_result * result) {
    if((x > 0 && y > 0) && x < 50.0 && y < 50.0) {
      /* Handle the easy case */
      gsl_sf_result gx, gy, gxy;
      gsl_sf_gamma_e(x, &gx);
      gsl_sf_gamma_e(y, &gy);
      gsl_sf_gamma_e(x+y, &gxy);
      result->val  = (gx.val*gy.val)/gxy.val;
      result->err  = gx.err * fabs(gy.val/gxy.val);
      result->err += gy.err * fabs(gx.val/gxy.val);
      result->err += fabs((gx.val*gy.val)/(gxy.val*gxy.val)) * gxy.err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if (isnegint(x) || isnegint(y)) {
      DOMAIN_ERROR(result);
    } else if (isnegint(x+y)) {  /* infinity in the denominator */
      result->val = 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    } else {
      gsl_sf_result lb;
      double sgn;
      int stat_lb = gsl_sf_lnbeta_sgn_e(x, y, &lb, &sgn);
      if(stat_lb == GSL_SUCCESS) {
	int status = gsl_sf_exp_err_e(lb.val, lb.err, result);
	result->val *= sgn;
	return status;
      }
      else {
	result->val = 0.0;
	result->err = 0.0;
	return stat_lb;
      }
    }
  }
  
  
  /*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
  
  double gsl_sf_lnbeta(const double x, const double y) {
    EVAL_RESULT(gsl_sf_lnbeta_e(x, y, &result));
  }
  
  double gsl_sf_beta(const double x, const double y) {
    EVAL_RESULT(gsl_sf_beta_e(x, y, &result));
  }

} // beta.c

// beta_inc.c ///////////////////////////////////////////////////////////////////////
extern "C" { 
  static int beta_cont_frac(const double a, const double b, const double x, gsl_sf_result * result)
  {
    const unsigned int max_iter = 512;        /* control iterations      */
    const double cutoff = 2.0 * GSL_DBL_MIN;  /* control the zero cutoff */
    unsigned int iter_count = 0;
    double cf;
    
    /* standard initialization for continued fraction */
    double num_term = 1.0;
    double den_term = 1.0 - (a+b)*x/(a+1.0);
    if (fabs(den_term) < cutoff) den_term = cutoff;
    den_term = 1.0/den_term;
    cf = den_term;
    
    while(iter_count < max_iter) {
      const int k  = iter_count + 1;
      double coeff = k*(b-k)*x/(((a-1.0)+2*k)*(a+2*k));
      double delta_frac;
      
      /* first step */
      den_term = 1.0 + coeff*den_term;
      num_term = 1.0 + coeff/num_term;
      if(fabs(den_term) < cutoff) den_term = cutoff;
      if(fabs(num_term) < cutoff) num_term = cutoff;
      den_term  = 1.0/den_term;
      
      delta_frac = den_term * num_term;
      cf *= delta_frac;
      
      coeff = -(a+k)*(a+b+k)*x/((a+2*k)*(a+2*k+1.0));
      
      /* second step */
      den_term = 1.0 + coeff*den_term;
      num_term = 1.0 + coeff/num_term;
      if(fabs(den_term) < cutoff) den_term = cutoff;
      if(fabs(num_term) < cutoff) num_term = cutoff;
      den_term = 1.0/den_term;
      
      delta_frac = den_term*num_term;
      cf *= delta_frac;
      
      if(fabs(delta_frac-1.0) < 2.0*GSL_DBL_EPSILON) break;
      
      ++iter_count;
    }
    
    result->val = cf;
    result->err = iter_count * 4.0 * GSL_DBL_EPSILON * fabs(cf);
    
    if(iter_count >= max_iter)
      GSL_ERROR ("error", GSL_EMAXITER);
    else
      return GSL_SUCCESS;
  }
  
  /*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

  int gsl_sf_beta_inc_e(const double a, const double b, const double x, gsl_sf_result * result)
  {
    if(a <= 0.0 || b <= 0.0 || x < 0.0 || x > 1.0) {
      DOMAIN_ERROR(result);
    }
    else if(x == 0.0) {
      result->val = 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
    else if(x == 1.0) {
      result->val = 1.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
    else {
      gsl_sf_result ln_beta;
      gsl_sf_result ln_x;
      gsl_sf_result ln_1mx;
      gsl_sf_result prefactor;
      const int stat_ln_beta = gsl_sf_lnbeta_e(a, b, &ln_beta);
      const int stat_ln_1mx = gsl_sf_log_1plusx_e(-x, &ln_1mx);
      const int stat_ln_x = gsl_sf_log_e(x, &ln_x);
      const int stat_ln = GSL_ERROR_SELECT_3(stat_ln_beta, stat_ln_1mx, stat_ln_x);
      
      const double ln_pre_val = -ln_beta.val + a * ln_x.val + b * ln_1mx.val;
      const double ln_pre_err =  ln_beta.err + fabs(a*ln_x.err) + fabs(b*ln_1mx.err);
      const int stat_exp = gsl_sf_exp_err_e(ln_pre_val, ln_pre_err, &prefactor);
      
      if(stat_ln != GSL_SUCCESS) {
	result->val = 0.0;
	result->err = 0.0;
	GSL_ERROR ("error", GSL_ESANITY);
      }
      
      if(x < (a + 1.0)/(a+b+2.0)) {
	/* Apply continued fraction directly. */
	gsl_sf_result cf;
	const int stat_cf = beta_cont_frac(a, b, x, &cf);
	int stat;
	result->val = prefactor.val * cf.val / a;
	result->err = (fabs(prefactor.err * cf.val) + fabs(prefactor.val * cf.err))/a;
	
	stat = GSL_ERROR_SELECT_2(stat_exp, stat_cf);
	if(stat == GSL_SUCCESS) {
	  CHECK_UNDERFLOW(result);
	}
	return stat;
      }
      else {
	/* Apply continued fraction after hypergeometric transformation. */
	gsl_sf_result cf;
	const int stat_cf = beta_cont_frac(b, a, 1.0-x, &cf);
	int stat;
	const double term = prefactor.val * cf.val / b;
	result->val  = 1.0 - term;
	result->err  = fabs(prefactor.err * cf.val)/b;
	result->err += fabs(prefactor.val * cf.err)/b;
	result->err += 2.0 * GSL_DBL_EPSILON * (1.0 + fabs(term));
	stat = GSL_ERROR_SELECT_2(stat_exp, stat_cf);
	if(stat == GSL_SUCCESS) {
	  CHECK_UNDERFLOW(result);
	}
	return stat;
      }
    }
  }

  /*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

  double gsl_sf_beta_inc(const double a, const double b, const double x) {
    EVAL_RESULT(gsl_sf_beta_inc_e(a, b, x, &result));
  }

} // beta_inc.c

// elementary.c ///////////////////////////////////////////////////////////////////
extern "C" { 
  int gsl_sf_multiply_e(const double x, const double y, gsl_sf_result * result) {
    const double ax = fabs(x);
    const double ay = fabs(y);
    
    if(x == 0.0 || y == 0.0) {
      /* It is necessary to eliminate this immediately.
       */
      result->val = 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
    else if((ax <= 1.0 && ay >= 1.0) || (ay <= 1.0 && ax >= 1.0)) {
      /* Straddling 1.0 is always safe.
       */
      result->val = x*y;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      const double f = 1.0 - 2.0 * GSL_DBL_EPSILON;
      const double min = GSL_MIN_DBL(fabs(x), fabs(y));
      const double max = GSL_MAX_DBL(fabs(x), fabs(y));
      if(max < 0.9 * GSL_SQRT_DBL_MAX || min < (f * DBL_MAX)/max) {
	result->val = GSL_COERCE_DBL(x*y);
	result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
	CHECK_UNDERFLOW(result);
	return GSL_SUCCESS;
      }
      else {
	OVERFLOW_ERROR(result);
      }
    }
  }

  int gsl_sf_multiply_err_e(const double x, const double dx,
			    const double y, const double dy,
			    gsl_sf_result * result)
  {
    int status = gsl_sf_multiply_e(x, y, result);
    result->err += fabs(dx*y) + fabs(dy*x);
    return status;
  }

  /*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

  double gsl_sf_multiply(const double x, const double y) {
    EVAL_RESULT(gsl_sf_multiply_e(x, y, &result));
  }
} // elementary.c

// erfc.c ///////////////////////////////////////////////////////////////////////////
extern "C" {
#define LogRootPi_  0.57236494292470008706
  
  static double erfc8_sum(double x)
{
  /* estimates erfc(x) valid for 8 < x < 100 */
  /* This is based on index 5725 in Hart et al */
  
  static double P[] = {
    2.97886562639399288862,
    7.409740605964741794425,
    6.1602098531096305440906,
    5.019049726784267463450058,
    1.275366644729965952479585264,
    0.5641895835477550741253201704
  };
  static double Q[] = {
    3.3690752069827527677,
    9.608965327192787870698,
    17.08144074746600431571095,
    12.0489519278551290360340491,
    9.396034016235054150430579648,
    2.260528520767326969591866945,
    1.0
  };
  double num=0.0, den=0.0;
  int i;
  
  num = P[5];
  for (i=4; i>=0; --i) {
    num = x*num + P[i];
  }
  den = Q[6];
  for (i=5; i>=0; --i) {
    den = x*den + Q[i];
  }
  
  return num/den;
}
  
  inline static double erfc8(double x) {
    double e;
    e = erfc8_sum(x);
    e *= exp(-x*x);
    return e;
  }

  inline static double log_erfc8(double x)
  {
    double e;
    e = erfc8_sum(x);
    e = log(e) - x*x;
    return e;
  }

  /* Abramowitz+Stegun, 7.1.5 */
  static int erfseries(double x, gsl_sf_result * result)
  {
    double coef = x;
    double e    = coef;
    double del;
    int k;
    for (k=1; k<30; ++k) {
      coef *= -x*x/k;
      del   = coef/(2.0*k+1.0);
      e += del;
    }
    result->val = 2.0 / M_SQRTPI * e;
    result->err = 2.0 / M_SQRTPI * (fabs(del) + GSL_DBL_EPSILON);
    return GSL_SUCCESS;
  }

  /* Chebyshev fit for erfc((t+1)/2), -1 < t < 1  */
  static double erfc_xlt1_data[20] = {
    1.06073416421769980345174155056,
    -0.42582445804381043569204735291,
    0.04955262679620434040357683080,
    0.00449293488768382749558001242,
    -0.00129194104658496953494224761,
    -0.00001836389292149396270416979,
    0.00002211114704099526291538556,
    -5.23337485234257134673693179020e-7,
    -2.78184788833537885382530989578e-7,
    1.41158092748813114560316684249e-8,
    2.72571296330561699984539141865e-9,
    -2.06343904872070629406401492476e-10,
    -2.14273991996785367924201401812e-11,
    2.22990255539358204580285098119e-12,
    1.36250074650698280575807934155e-13,
    -1.95144010922293091898995913038e-14,
    -6.85627169231704599442806370690e-16,
    1.44506492869699938239521607493e-16,
    2.45935306460536488037576200030e-18,
    -9.29599561220523396007359328540e-19
  };

  static cheb_series erfc_xlt1_cs = {
    erfc_xlt1_data,
    19,
    -1, 1,
    12
  };
  
  /* Chebyshev fit for erfc(x) exp(x^2), 1 < x < 5, x = 2t + 3, -1 < t < 1
   */
  static double erfc_x15_data[25] = {
    0.44045832024338111077637466616,
    -0.143958836762168335790826895326,
    0.044786499817939267247056666937,
    -0.013343124200271211203618353102,
    0.003824682739750469767692372556,
    -0.001058699227195126547306482530,
    0.000283859419210073742736310108,
    -0.000073906170662206760483959432,
    0.000018725312521489179015872934,
    -4.62530981164919445131297264430e-6,
    1.11558657244432857487884006422e-6,
    -2.63098662650834130067808832725e-7,
    6.07462122724551777372119408710e-8,
    -1.37460865539865444777251011793e-8,
    3.05157051905475145520096717210e-9,
    -6.65174789720310713757307724790e-10,
    1.42483346273207784489792999706e-10,
    -3.00141127395323902092018744545e-11,
    6.22171792645348091472914001250e-12,
    -1.26994639225668496876152836555e-12,
    2.55385883033257575402681845385e-13,
    -5.06258237507038698392265499770e-14,
    9.89705409478327321641264227110e-15,
    -1.90685978789192181051961024995e-15,
    3.50826648032737849245113757340e-16
  };

  static cheb_series erfc_x15_cs = {
    erfc_x15_data,
    24,
    -1, 1,
    16
  };
  
  /* Chebyshev fit for erfc(x) x exp(x^2), 5 < x < 10, x = (5t + 15)/2, -1 < t < 1  */
  static double erfc_x510_data[20] = {
    1.11684990123545698684297865808,
    0.003736240359381998520654927536,
    -0.000916623948045470238763619870,
    0.000199094325044940833965078819,
    -0.000040276384918650072591781859,
    7.76515264697061049477127605790e-6,
    -1.44464794206689070402099225301e-6,
    2.61311930343463958393485241947e-7,
    -4.61833026634844152345304095560e-8,
    8.00253111512943601598732144340e-9,
    -1.36291114862793031395712122089e-9,
    2.28570483090160869607683087722e-10,
    -3.78022521563251805044056974560e-11,
    6.17253683874528285729910462130e-12,
    -9.96019290955316888445830597430e-13,
    1.58953143706980770269506726000e-13,
    -2.51045971047162509999527428316e-14,
    3.92607828989125810013581287560e-15,
    -6.07970619384160374392535453420e-16,
    9.12600607264794717315507477670e-17
  };

  static cheb_series erfc_x510_cs = {
    erfc_x510_data,
    19,
    -1, 1,
    12
  };

  /*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

  int gsl_sf_erfc_e(double x, gsl_sf_result * result) {
    const double ax = fabs(x);
    double e_val, e_err;
    
    /* CHECK_POINTER(result) */
    
    if(ax <= 1.0) {
      double t = 2.0*ax - 1.0;
      gsl_sf_result c;
      cheb_eval_e(&erfc_xlt1_cs, t, &c);
      e_val = c.val;
      e_err = c.err;
    }
    else if(ax <= 5.0) {
      double ex2 = exp(-x*x);
      double t = 0.5*(ax-3.0);
      gsl_sf_result c;
      cheb_eval_e(&erfc_x15_cs, t, &c);
      e_val = ex2 * c.val;
      e_err = ex2 * (c.err + 2.0*fabs(x)*GSL_DBL_EPSILON);
    }
    else if(ax < 10.0) {
      double exterm = exp(-x*x) / ax;
      double t = (2.0*ax - 15.0)/5.0;
      gsl_sf_result c;
      cheb_eval_e(&erfc_x510_cs, t, &c);
      e_val = exterm * c.val;
      e_err = exterm * (c.err + 2.0*fabs(x)*GSL_DBL_EPSILON + GSL_DBL_EPSILON);
    }
    else {
      e_val = erfc8(ax);
      e_err = (x*x + 1.0) * GSL_DBL_EPSILON * fabs(e_val);
    }
    
    if(x < 0.0) {
      result->val  = 2.0 - e_val;
      result->err  = e_err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    }
    else {
      result->val  = e_val;
      result->err  = e_err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    }
    
    return GSL_SUCCESS;
  }

  int gsl_sf_log_erfc_e(double x, gsl_sf_result * result) {
    /* CHECK_POINTER(result) */
    
    if(x*x < 10.0*GSL_ROOT6_DBL_EPSILON) {
      const double y = x / M_SQRTPI;
      /* series for -1/2 Log[Erfc[Sqrt[Pi] y]] */
      const double c3 = (4.0 - M_PI)/3.0;
      const double c4 = 2.0*(1.0 - M_PI/3.0);
      const double c5 = -0.001829764677455021;  /* (96.0 - 40.0*M_PI + 3.0*M_PI*M_PI)/30.0  */
      const double c6 =  0.02629651521057465;   /* 2.0*(120.0 - 60.0*M_PI + 7.0*M_PI*M_PI)/45.0 */
      const double c7 = -0.01621575378835404;
      const double c8 =  0.00125993961762116;
      const double c9 =  0.00556964649138;
      const double c10 = -0.0045563339802;
      const double c11 =  0.0009461589032;
      const double c12 =  0.0013200243174;
      const double c13 = -0.00142906;
      const double c14 =  0.00048204;
      double series = c8 + y*(c9 + y*(c10 + y*(c11 + y*(c12 + y*(c13 + c14*y)))));
      series = y*(1.0 + y*(1.0 + y*(c3 + y*(c4 + y*(c5 + y*(c6 + y*(c7 + y*series)))))));
      result->val = -2.0 * series;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if(x > 8.0) {
      result->val = log_erfc8(x);
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      gsl_sf_result result_erfc;
      gsl_sf_erfc_e(x, &result_erfc);
      result->val  = log(result_erfc.val);
      result->err  = fabs(result_erfc.err / result_erfc.val);
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
  }

  int gsl_sf_erf_e(double x, gsl_sf_result * result) {
    /* CHECK_POINTER(result) */
    
    if(fabs(x) < 1.0) {
      return erfseries(x, result);
    }
    else {
      gsl_sf_result result_erfc;
      gsl_sf_erfc_e(x, &result_erfc);
      result->val  = 1.0 - result_erfc.val;
      result->err  = result_erfc.err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
  }

  int gsl_sf_erf_Z_e(double x, gsl_sf_result * result) {
    /* CHECK_POINTER(result) */
    {
      const double ex2 = exp(-x*x/2.0);
      result->val  = ex2 / (M_SQRT2 * M_SQRTPI);
      result->err  = fabs(x * result->val) * GSL_DBL_EPSILON;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      CHECK_UNDERFLOW(result);
      return GSL_SUCCESS;
    }
  }
  
  
  int gsl_sf_erf_Q_e(double x, gsl_sf_result * result) {
    /* CHECK_POINTER(result) */
    {
      gsl_sf_result result_erfc;
      int stat = gsl_sf_erfc_e(x/M_SQRT2, &result_erfc);
      result->val  = 0.5 * result_erfc.val;
      result->err  = 0.5 * result_erfc.err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return stat;
    }
  }

  int gsl_sf_hazard_e(double x, gsl_sf_result * result) {
    if(x < 25.0)
      {
	gsl_sf_result result_ln_erfc;
	const int stat_l = gsl_sf_log_erfc_e(x/M_SQRT2, &result_ln_erfc);
	const double lnc = -0.22579135264472743236; /* ln(sqrt(2/pi)) */
	const double arg = lnc - 0.5*x*x - result_ln_erfc.val;
	const int stat_e = gsl_sf_exp_e(arg, result);
	result->err += 3.0 * (1.0 + fabs(x)) * GSL_DBL_EPSILON * fabs(result->val);
	result->err += fabs(result_ln_erfc.err * result->val);
	return GSL_ERROR_SELECT_2(stat_l, stat_e);
      }
    else
      {
	const double ix2 = 1.0/(x*x);
	const double corrB = 1.0 - 9.0*ix2 * (1.0 - 11.0*ix2);
	const double corrM = 1.0 - 5.0*ix2 * (1.0 - 7.0*ix2 * corrB);
	const double corrT = 1.0 - ix2 * (1.0 - 3.0*ix2*corrM);
	result->val = x / corrT;
	result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
	return GSL_SUCCESS;
      }
  }

  /*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/


  double gsl_sf_erfc(double x) {
    EVAL_RESULT(gsl_sf_erfc_e(x, &result));
  }

  double gsl_sf_log_erfc(double x) {
    EVAL_RESULT(gsl_sf_log_erfc_e(x, &result));
  }

  double gsl_sf_erf(double x) {
    EVAL_RESULT(gsl_sf_erf_e(x, &result));
  }

  double gsl_sf_erf_Z(double x) {
    EVAL_RESULT(gsl_sf_erf_Z_e(x, &result));
  }

  double gsl_sf_erf_Q(double x) {
    EVAL_RESULT(gsl_sf_erf_Q_e(x, &result));
  }

  double gsl_sf_hazard(double x) {
    EVAL_RESULT(gsl_sf_hazard_e(x, &result));
  }

} // erfc.c

// error.c ////////////////////////////////////////////////////////////////////////
extern "C" {
  gsl_error_handler_t * gsl_error_handler = NULL;

  static void no_error_handler (const char *reason, const char *file, int line, int gsl_errno);

  void gsl_error (const char * reason, const char * file, int line, int gsl_errno) {
    if (gsl_error_handler) 
      {
	(*gsl_error_handler) (reason, file, line, gsl_errno);
	return ;
      }
    
    gsl_stream_printf ("ERROR", file, line, reason);
    
    fflush (stdout);
    fprintf (stderr, "Default GSL error handler invoked.\n");
    fflush (stderr);
    
    abort ();
  }

  gsl_error_handler_t * gsl_set_error_handler (gsl_error_handler_t * new_handler) {
    gsl_error_handler_t * previous_handler = gsl_error_handler;
    gsl_error_handler = new_handler;
    return previous_handler;
  }

  gsl_error_handler_t * gsl_set_error_handler_off (void) {
    gsl_error_handler_t * previous_handler = gsl_error_handler;
    gsl_error_handler = no_error_handler;
    return previous_handler;
  }

  static void no_error_handler (const char *reason, const char *file, int line, int gsl_errno) {
    /* do nothing */
    reason = 0;
    file = 0;
    line = 0;
    gsl_errno = 0;
    return;
  }

} // error.c

// exp.c ///////////////////////////////////////////////////////////////////////
extern "C" { 
  /* Evaluate the continued fraction for exprel.
   * [Abramowitz+Stegun, 4.2.41]
   */
  static int exprel_n_CF(const int N, const double x, gsl_sf_result * result) {
    const double RECUR_BIG = GSL_SQRT_DBL_MAX;
    const int maxiter = 5000;
    int n = 1;
    double Anm2 = 1.0;
    double Bnm2 = 0.0;
    double Anm1 = 0.0;
    double Bnm1 = 1.0;
    double a1 = 1.0;
    double b1 = 1.0;
    double a2 = -x;
    double b2 = N+1;
    double an, bn;
    
    double fn;
    
    double An = b1*Anm1 + a1*Anm2;   /* A1 */
    double Bn = b1*Bnm1 + a1*Bnm2;   /* B1 */
  
    /* One explicit step, before we get to the main pattern. */
    n++;
    Anm2 = Anm1;
    Bnm2 = Bnm1;
    Anm1 = An;
    Bnm1 = Bn;
    An = b2*Anm1 + a2*Anm2;   /* A2 */
    Bn = b2*Bnm1 + a2*Bnm2;   /* B2 */
    
    fn = An/Bn;
    
    while(n < maxiter) {
      double old_fn;
      double del;
      n++;
      Anm2 = Anm1;
      Bnm2 = Bnm1;
      Anm1 = An;
      Bnm1 = Bn;
      an = ( GSL_IS_ODD(n) ? ((n-1)/2)*x : -(N+(n/2)-1)*x );
      bn = N + n - 1;
      An = bn*Anm1 + an*Anm2;
      Bn = bn*Bnm1 + an*Bnm2;
      
      if(fabs(An) > RECUR_BIG || fabs(Bn) > RECUR_BIG) {
	An /= RECUR_BIG;
	Bn /= RECUR_BIG;
	Anm1 /= RECUR_BIG;
	Bnm1 /= RECUR_BIG;
	Anm2 /= RECUR_BIG;
	Bnm2 /= RECUR_BIG;
      }
      
      old_fn = fn;
      fn = An/Bn;
      del = old_fn/fn;
      
      if(fabs(del - 1.0) < 2.0*GSL_DBL_EPSILON) break;
    }
    
    result->val = fn;
    result->err = 2.0*(n+1.0)*GSL_DBL_EPSILON*fabs(fn);
    
    if(n == maxiter)
      GSL_ERROR ("error", GSL_EMAXITER);
    else
      return GSL_SUCCESS;
  }
  
  
  /*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

  int gsl_sf_exp_e(const double x, gsl_sf_result * result) {
    if(x > GSL_LOG_DBL_MAX) {
      OVERFLOW_ERROR(result);
    }
    else if(x < GSL_LOG_DBL_MIN) {
      UNDERFLOW_ERROR(result);
    }
    else {
      result->val = exp(x);
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
  }
  
  int gsl_sf_exp_e10_e(const double x, gsl_sf_result_e10 * result)
  {
    if(x > INT_MAX-1) {
      OVERFLOW_ERROR_E10(result);
    }
    else if(x < INT_MIN+1) {
      UNDERFLOW_ERROR_E10(result);
    }
    else {
      const int N = (int) floor(x/M_LN10);
      result->val = exp(x-N*M_LN10);
      result->err = 2.0 * (fabs(x)+1.0) * GSL_DBL_EPSILON * fabs(result->val);
      result->e10 = N;
      return GSL_SUCCESS;
    }
  }
  
  int gsl_sf_exp_mult_e(const double x, const double y, gsl_sf_result * result) {
    const double ay  = fabs(y);
    
    if(y == 0.0) {
      result->val = 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
    else if(   ( x < 0.5*GSL_LOG_DBL_MAX   &&   x > 0.5*GSL_LOG_DBL_MIN)
	       && (ay < 0.8*GSL_SQRT_DBL_MAX  &&  ay > 1.2*GSL_SQRT_DBL_MIN)
	       ) {
      const double ex = exp(x);
      result->val = y * ex;
      result->err = (2.0 + fabs(x)) * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      const double ly  = log(ay);
      const double lnr = x + ly;
      
      if(lnr > GSL_LOG_DBL_MAX - 0.01) {
	OVERFLOW_ERROR(result);
      }
      else if(lnr < GSL_LOG_DBL_MIN + 0.01) {
	UNDERFLOW_ERROR(result);
      }
      else {
	const double sy   = GSL_SIGN(y);
	const double M    = floor(x);
	const double N    = floor(ly);
	const double a    = x  - M;
	const double b    = ly - N;
	const double berr = 2.0 * GSL_DBL_EPSILON * (fabs(ly) + fabs(N));
	result->val  = sy * exp(M+N) * exp(a+b);
	result->err  = berr * fabs(result->val);
	result->err += 2.0 * GSL_DBL_EPSILON * (M + N + 1.0) * fabs(result->val);
	return GSL_SUCCESS;
      }
    }
  }

  int gsl_sf_exp_mult_e10_e(const double x, const double y, gsl_sf_result_e10 * result) {
    const double ay  = fabs(y);
    
    if(y == 0.0) {
      result->val = 0.0;
      result->err = 0.0;
      result->e10 = 0;
      return GSL_SUCCESS;
    }
    else if(   ( x < 0.5*GSL_LOG_DBL_MAX   &&   x > 0.5*GSL_LOG_DBL_MIN)
	       && (ay < 0.8*GSL_SQRT_DBL_MAX  &&  ay > 1.2*GSL_SQRT_DBL_MIN)
	       ) {
      const double ex = exp(x);
      result->val = y * ex;
      result->err = (2.0 + fabs(x)) * GSL_DBL_EPSILON * fabs(result->val);
      result->e10 = 0;
      return GSL_SUCCESS;
    }
    else {
      const double ly  = log(ay);
      const double l10_val = (x + ly)/M_LN10;
      
      if(l10_val > INT_MAX-1) {
	OVERFLOW_ERROR_E10(result);
      }
      else if(l10_val < INT_MIN+1) {
	UNDERFLOW_ERROR_E10(result);
      }
      else {
	const double sy  = GSL_SIGN(y);
	const int    N   = (int) floor(l10_val);
	const double arg_val = (l10_val - N) * M_LN10;
	const double arg_err = 2.0 * GSL_DBL_EPSILON * fabs(ly);
	
	result->val  = sy * exp(arg_val);
	result->err  = arg_err * fabs(result->val);
	result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
	result->e10 = N;
	
	return GSL_SUCCESS;
      }
    }
  }

  int gsl_sf_exp_mult_err_e(const double x, const double dx,
			    const double y, const double dy,
			    gsl_sf_result * result)
  {
    const double ay  = fabs(y);
    
    if(y == 0.0) {
      result->val = 0.0;
      result->err = fabs(dy * exp(x));
      return GSL_SUCCESS;
    }
    else if(   ( x < 0.5*GSL_LOG_DBL_MAX   &&   x > 0.5*GSL_LOG_DBL_MIN)
	       && (ay < 0.8*GSL_SQRT_DBL_MAX  &&  ay > 1.2*GSL_SQRT_DBL_MIN)
	       ) {
      double ex = exp(x);
      result->val  = y * ex;
      result->err  = ex * (fabs(dy) + fabs(y*dx));
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      const double ly  = log(ay);
      const double lnr = x + ly;
      
      if(lnr > GSL_LOG_DBL_MAX - 0.01) {
	OVERFLOW_ERROR(result);
      }
      else if(lnr < GSL_LOG_DBL_MIN + 0.01) {
	UNDERFLOW_ERROR(result);
      }
      else {
	const double sy  = GSL_SIGN(y);
	const double M   = floor(x);
	const double N   = floor(ly);
	const double a   = x  - M;
	const double b   = ly - N;
	const double eMN = exp(M+N);
	const double eab = exp(a+b);
	result->val  = sy * eMN * eab;
	result->err  = eMN * eab * 2.0*GSL_DBL_EPSILON;
	result->err += eMN * eab * fabs(dy/y);
	result->err += eMN * eab * fabs(dx);
	return GSL_SUCCESS;
      }
    }
  }

  int gsl_sf_exp_mult_err_e10_e(const double x, const double dx,
				const double y, const double dy,
				gsl_sf_result_e10 * result)
  {
    const double ay  = fabs(y);
    
    if(y == 0.0) {
      result->val = 0.0;
      result->err = fabs(dy * exp(x));
      result->e10 = 0;
      return GSL_SUCCESS;
    }
    else if(   ( x < 0.5*GSL_LOG_DBL_MAX   &&   x > 0.5*GSL_LOG_DBL_MIN)
	       && (ay < 0.8*GSL_SQRT_DBL_MAX  &&  ay > 1.2*GSL_SQRT_DBL_MIN)
	       ) {
      const double ex = exp(x);
      result->val  = y * ex;
      result->err  = ex * (fabs(dy) + fabs(y*dx));
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      result->e10 = 0;
      return GSL_SUCCESS;
    }
    else {
      const double ly  = log(ay);
      const double l10_val = (x + ly)/M_LN10;
      
      if(l10_val > INT_MAX-1) {
	OVERFLOW_ERROR_E10(result);
      }
      else if(l10_val < INT_MIN+1) {
	UNDERFLOW_ERROR_E10(result);
      }
      else {
	const double sy  = GSL_SIGN(y);
	const int    N   = (int) floor(l10_val);
	const double arg_val = (l10_val - N) * M_LN10;
	const double arg_err = dy/fabs(y) + dx + 2.0*GSL_DBL_EPSILON*fabs(arg_val);
	
	result->val  = sy * exp(arg_val);
	result->err  = arg_err * fabs(result->val);
	result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
	result->e10 = N;
	
	return GSL_SUCCESS;
      }
    }
  }

  int gsl_sf_expm1_e(const double x, gsl_sf_result * result) {
    const double cut = 0.002;
    
    if(x < GSL_LOG_DBL_MIN) {
      result->val = -1.0;
      result->err = GSL_DBL_EPSILON;
      return GSL_SUCCESS;
    }
    else if(x < -cut) {
      result->val = exp(x) - 1.0;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if(x < cut) {
      result->val = x * (1.0 + 0.5*x*(1.0 + x/3.0*(1.0 + 0.25*x*(1.0 + 0.2*x))));
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    } 
    else if(x < GSL_LOG_DBL_MAX) {
      result->val = exp(x) - 1.0;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      OVERFLOW_ERROR(result);
    }
  }

  int gsl_sf_exprel_e(const double x, gsl_sf_result * result) {
    const double cut = 0.002;
    
    if(x < GSL_LOG_DBL_MIN) {
      result->val = -1.0/x;
      result->err = GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if(x < -cut) {
      result->val = (exp(x) - 1.0)/x;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if(x < cut) {
      result->val = (1.0 + 0.5*x*(1.0 + x/3.0*(1.0 + 0.25*x*(1.0 + 0.2*x))));
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    } 
    else if(x < GSL_LOG_DBL_MAX) {
      result->val = (exp(x) - 1.0)/x;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      OVERFLOW_ERROR(result);
    }
  }
  
  int gsl_sf_exprel_2_e(double x, gsl_sf_result * result) {
    const double cut = 0.002;
    
    if(x < GSL_LOG_DBL_MIN) {
      result->val = -2.0/x*(1.0 + 1.0/x);
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if(x < -cut) {
      result->val = 2.0*(exp(x) - 1.0 - x)/(x*x);
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if(x < cut) {
      result->val = (1.0 + 1.0/3.0*x*(1.0 + 0.25*x*(1.0 + 0.2*x*(1.0 + 1.0/6.0*x))));
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    } 
    else if(x < GSL_LOG_DBL_MAX) {
      result->val = 2.0*(exp(x) - 1.0 - x)/(x*x);
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      OVERFLOW_ERROR(result);
    }
  }

  int gsl_sf_exprel_n_e(const int N, const double x, gsl_sf_result * result) {
    if(N < 0) {
      DOMAIN_ERROR(result);
    }
    else if(x == 0.0) {
      result->val = 1.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
    else if(fabs(x) < GSL_ROOT3_DBL_EPSILON * N) {
      result->val = 1.0 + x/(N+1) * (1.0 + x/(N+2));
      result->err = 2.0 * GSL_DBL_EPSILON;
      return GSL_SUCCESS;
    }
    else if(N == 0) {
      return gsl_sf_exp_e(x, result);
    }
    else if(N == 1) {
      return gsl_sf_exprel_e(x, result);
    }
    else if(N == 2) {
      return gsl_sf_exprel_2_e(x, result);
    }
    else {
      if(x > N && (-x + N*(1.0 + log(x/N)) < GSL_LOG_DBL_EPSILON)) {
	/* x is much larger than n.
	 * Ignore polynomial part, so
	 * exprel_N(x) ~= e^x N!/x^N
	 */
	gsl_sf_result lnf_N;
	double lnr_val;
	double lnr_err;
	double lnterm;
	gsl_sf_lnfact_e(N, &lnf_N);
	lnterm = N*log(x);
	lnr_val  = x + lnf_N.val - lnterm;
	lnr_err  = GSL_DBL_EPSILON * (fabs(x) + fabs(lnf_N.val) + fabs(lnterm));
	lnr_err += lnf_N.err;
	return gsl_sf_exp_err_e(lnr_val, lnr_err, result);
      }
      else if(x > N) {
	/* Write the identity
	 *   exprel_n(x) = e^x n! / x^n (1 - Gamma[n,x]/Gamma[n])
	 * then use the asymptotic expansion
	 * Gamma[n,x] ~ x^(n-1) e^(-x) (1 + (n-1)/x + (n-1)(n-2)/x^2 + ...)
	 */
	double ln_x = log(x);
	gsl_sf_result lnf_N;
	double lg_N;
	double lnpre_val;
	double lnpre_err;
	gsl_sf_lnfact_e(N, &lnf_N);    /* log(N!)       */
	lg_N  = lnf_N.val - log(N);       /* log(Gamma(N)) */
	lnpre_val  = x + lnf_N.val - N*ln_x;
	lnpre_err  = GSL_DBL_EPSILON * (fabs(x) + fabs(lnf_N.val) + fabs(N*ln_x));
	lnpre_err += lnf_N.err;
	if(lnpre_val < GSL_LOG_DBL_MAX - 5.0) {
	  int stat_eG;
	  gsl_sf_result bigG_ratio;
	  gsl_sf_result pre;
	  int stat_ex = gsl_sf_exp_err_e(lnpre_val, lnpre_err, &pre);
	  double ln_bigG_ratio_pre = -x + (N-1)*ln_x - lg_N;
	  double bigGsum = 1.0;
	  double term = 1.0;
	  int k;
	  for(k=1; k<N; k++) {
	    term *= (N-k)/x;
	    bigGsum += term;
	  }
	  stat_eG = gsl_sf_exp_mult_e(ln_bigG_ratio_pre, bigGsum, &bigG_ratio);
	  if(stat_eG == GSL_SUCCESS) {
	    result->val  = pre.val * (1.0 - bigG_ratio.val);
	    result->err  = pre.val * (2.0*GSL_DBL_EPSILON + bigG_ratio.err);
	    result->err += pre.err * fabs(1.0 - bigG_ratio.val);
	    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
	    return stat_ex;
	  }
	  else {
	    result->val = 0.0;
	    result->err = 0.0;
	    return stat_eG;
	  }
	}
	else {
	  OVERFLOW_ERROR(result);
	}
      }
      else if(x > -10.0*N) {
	return exprel_n_CF(N, x, result);
      }
      else {
	/* x -> -Inf asymptotic:
	 * exprel_n(x) ~ e^x n!/x^n - n/x (1 + (n-1)/x + (n-1)(n-2)/x + ...)
	 *             ~ - n/x (1 + (n-1)/x + (n-1)(n-2)/x + ...)
	 */
	double sum  = 1.0;
	double term = 1.0;
	int k;
	for(k=1; k<N; k++) {
	  term *= (N-k)/x;
	  sum  += term;
	}
	result->val = -N/x * sum;
	result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
	return GSL_SUCCESS;
      }
    }
  }

  int gsl_sf_exp_err_e(const double x, const double dx, gsl_sf_result * result)
  {
    const double adx = fabs(dx);
    
    /* CHECK_POINTER(result) */
    
    if(x + adx > GSL_LOG_DBL_MAX) {
      OVERFLOW_ERROR(result);
    }
    else if(x - adx < GSL_LOG_DBL_MIN) {
      UNDERFLOW_ERROR(result);
    }
    else {
      const double ex  = exp(x);
      const double edx = exp(adx);
      result->val  = ex;
      result->err  = ex * GSL_MAX_DBL(GSL_DBL_EPSILON, edx - 1.0/edx);
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
  }

  int gsl_sf_exp_err_e10_e(const double x, const double dx, gsl_sf_result_e10 * result)
  {
    const double adx = fabs(dx);
    
    /* CHECK_POINTER(result) */
    
    if(x + adx > INT_MAX - 1) {
      OVERFLOW_ERROR_E10(result);
    }
    else if(x - adx < INT_MIN + 1) {
      UNDERFLOW_ERROR_E10(result);
    }
    else {
      const int    N  = (int)floor(x/M_LN10);
      const double ex = exp(x-N*M_LN10);
      result->val = ex;
      result->err = ex * (2.0 * GSL_DBL_EPSILON * (fabs(x) + 1.0) + adx);
      result->e10 = N;
      return GSL_SUCCESS;
    }
  }
  
  /*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
  
  double gsl_sf_exp(const double x)
  {
    EVAL_RESULT(gsl_sf_exp_e(x, &result));
  }
  
  double gsl_sf_exp_mult(const double x, const double y)
  {
    EVAL_RESULT(gsl_sf_exp_mult_e(x, y, &result));
  }
  
  double gsl_sf_expm1(const double x)
  {
    EVAL_RESULT(gsl_sf_expm1_e(x, &result));
  }
  
  double gsl_sf_exprel(const double x)
  {
    EVAL_RESULT(gsl_sf_exprel_e(x, &result));
  }
  
  double gsl_sf_exprel_2(const double x)
  {
    EVAL_RESULT(gsl_sf_exprel_2_e(x, &result));
  }
  
  double gsl_sf_exprel_n(const int n, const double x)
  {
    EVAL_RESULT(gsl_sf_exprel_n_e(n, x, &result));
  }

} // exp.c

// expint.c ///////////////////////////////////////////////////////////////////////

extern "C" { 
  
  /*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

  /*
    Chebyshev expansions: based on SLATEC e1.f, W. Fullerton
    
    Series for AE11       on the interval -1.00000D-01 to  0.
    with weighted error   1.76E-17
    log weighted error  16.75
    significant figures required  15.70
    decimal places required  17.55
    

    Series for AE12       on the interval -2.50000D-01 to -1.00000D-01
    with weighted error   5.83E-17
    log weighted error  16.23
    significant figures required  15.76
    decimal places required  16.93
    
    
    Series for E11        on the interval -4.00000D+00 to -1.00000D+00
    with weighted error   1.08E-18
    log weighted error  17.97
    significant figures required  19.02
    decimal places required  18.61
    
    
    Series for E12        on the interval -1.00000D+00 to  1.00000D+00
    with weighted error   3.15E-18
    log weighted error  17.50
    approx significant figures required  15.8
    decimal places required  18.10
    
    
    Series for AE13       on the interval  2.50000D-01 to  1.00000D+00
    with weighted error   2.34E-17
    log weighted error  16.63
    significant figures required  16.14
    decimal places required  17.33
    
    
    Series for AE14       on the interval  0.          to  2.50000D-01
    with weighted error   5.41E-17
    log weighted error  16.27
    significant figures required  15.38
    decimal places required  16.97
  */

  static double AE11_data[39] = {
    0.121503239716065790,
    -0.065088778513550150,
    0.004897651357459670,
    -0.000649237843027216,
    0.000093840434587471,
    0.000000420236380882,
    -0.000008113374735904,
    0.000002804247688663,
    0.000000056487164441,
    -0.000000344809174450,
    0.000000058209273578,
    0.000000038711426349,
    -0.000000012453235014,
    -0.000000005118504888,
    0.000000002148771527,
    0.000000000868459898,
    -0.000000000343650105,
    -0.000000000179796603,
    0.000000000047442060,
    0.000000000040423282,
    -0.000000000003543928,
    -0.000000000008853444,
    -0.000000000000960151,
    0.000000000001692921,
    0.000000000000607990,
    -0.000000000000224338,
    -0.000000000000200327,
    -0.000000000000006246,
    0.000000000000045571,
    0.000000000000016383,
    -0.000000000000005561,
    -0.000000000000006074,
    -0.000000000000000862,
    0.000000000000001223,
    0.000000000000000716,
    -0.000000000000000024,
    -0.000000000000000201,
    -0.000000000000000082,
    0.000000000000000017
  };
  static cheb_series AE11_cs = {
    AE11_data,
    38,
    -1, 1,
    20
  };
  
  static double AE12_data[25] = {
    0.582417495134726740,
    -0.158348850905782750,
    -0.006764275590323141,
    0.005125843950185725,
    0.000435232492169391,
    -0.000143613366305483,
    -0.000041801320556301,
    -0.000002713395758640,
    0.000001151381913647,
    0.000000420650022012,
    0.000000066581901391,
    0.000000000662143777,
    -0.000000002844104870,
    -0.000000000940724197,
    -0.000000000177476602,
    -0.000000000015830222,
    0.000000000002905732,
    0.000000000001769356,
    0.000000000000492735,
    0.000000000000093709,
    0.000000000000010707,
    -0.000000000000000537,
    -0.000000000000000716,
    -0.000000000000000244,
    -0.000000000000000058
  };
  static cheb_series AE12_cs = {
    AE12_data,
    24,
    -1, 1,
    15
  };
  
  static double E11_data[19] = {
    -16.11346165557149402600,
    7.79407277874268027690,
    -1.95540581886314195070,
    0.37337293866277945612,
    -0.05692503191092901938,
    0.00721107776966009185,
    -0.00078104901449841593,
    0.00007388093356262168,
    -0.00000620286187580820,
    0.00000046816002303176,
    -0.00000003209288853329,
    0.00000000201519974874,
    -0.00000000011673686816,
    0.00000000000627627066,
    -0.00000000000031481541,
    0.00000000000001479904,
    -0.00000000000000065457,
    0.00000000000000002733,
    -0.00000000000000000108
  };
  static cheb_series E11_cs = {
    E11_data,
    18,
    -1, 1,
    13
  };
  
  static double E12_data[16] = {
    -0.03739021479220279500,
    0.04272398606220957700,
    -0.13031820798497005440,
    0.01441912402469889073,
    -0.00134617078051068022,
    0.00010731029253063780,
    -0.00000742999951611943,
    0.00000045377325690753,
    -0.00000002476417211390,
    0.00000000122076581374,
    -0.00000000005485141480,
    0.00000000000226362142,
    -0.00000000000008635897,
    0.00000000000000306291,
    -0.00000000000000010148,
    0.00000000000000000315
  };
  static cheb_series E12_cs = {
    E12_data,
    15,
    -1, 1,
    10
  };
  
  static double AE13_data[25] = {
    -0.605773246640603460,
    -0.112535243483660900,
    0.013432266247902779,
    -0.001926845187381145,
    0.000309118337720603,
    -0.000053564132129618,
    0.000009827812880247,
    -0.000001885368984916,
    0.000000374943193568,
    -0.000000076823455870,
    0.000000016143270567,
    -0.000000003466802211,
    0.000000000758754209,
    -0.000000000168864333,
    0.000000000038145706,
    -0.000000000008733026,
    0.000000000002023672,
    -0.000000000000474132,
    0.000000000000112211,
    -0.000000000000026804,
    0.000000000000006457,
    -0.000000000000001568,
    0.000000000000000383,
    -0.000000000000000094,
    0.000000000000000023
  };
  static cheb_series AE13_cs = {
    AE13_data,
    24,
    -1, 1,
    15
  };
  
  static double AE14_data[26] = {
    -0.18929180007530170,
    -0.08648117855259871,
    0.00722410154374659,
    -0.00080975594575573,
    0.00010999134432661,
    -0.00001717332998937,
    0.00000298562751447,
    -0.00000056596491457,
    0.00000011526808397,
    -0.00000002495030440,
    0.00000000569232420,
    -0.00000000135995766,
    0.00000000033846628,
    -0.00000000008737853,
    0.00000000002331588,
    -0.00000000000641148,
    0.00000000000181224,
    -0.00000000000052538,
    0.00000000000015592,
    -0.00000000000004729,
    0.00000000000001463,
    -0.00000000000000461,
    0.00000000000000148,
    -0.00000000000000048,
    0.00000000000000016,
    -0.00000000000000005
  };
  static cheb_series AE14_cs = {
    AE14_data,
    25,
    -1, 1,
    13
  };
  
 
  /* implementation for E1, allowing for scaling by exp(x) */
  static int expint_E1_impl(const double x, gsl_sf_result * result, const int scale) {
    const double xmaxt = -GSL_LOG_DBL_MIN;      /* XMAXT = -LOG (R1MACH(1)) */
    const double xmax  = xmaxt - log(xmaxt);    /* XMAX = XMAXT - LOG(XMAXT) */
    
    /* CHECK_POINTER(result) */
    
    if(x < -xmax && !scale) {
      OVERFLOW_ERROR(result);
    }
    else if(x <= -10.0) {
      const double s = 1.0/x * ( scale ? 1.0 : exp(-x) );
      gsl_sf_result result_c;
      cheb_eval_e(&AE11_cs, 20.0/x+1.0, &result_c);
      result->val  = s * (1.0 + result_c.val);
      result->err  = s * result_c.err;
      result->err += 2.0 * GSL_DBL_EPSILON * (fabs(x) + 1.0) * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if(x <= -4.0) {
      const double s = 1.0/x * ( scale ? 1.0 : exp(-x) );
      gsl_sf_result result_c;
      cheb_eval_e(&AE12_cs, (40.0/x+7.0)/3.0, &result_c);
      result->val  = s * (1.0 + result_c.val);
      result->err  = s * result_c.err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if(x <= -1.0) {
      const double ln_term = -log(fabs(x));
      const double scale_factor = ( scale ? exp(x) : 1.0 );
      gsl_sf_result result_c;
      cheb_eval_e(&E11_cs, (2.0*x+5.0)/3.0, &result_c);
      result->val  = scale_factor * (ln_term + result_c.val);
      result->err  = scale_factor * (result_c.err + GSL_DBL_EPSILON * fabs(ln_term));
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if(x == 0.0) {
      DOMAIN_ERROR(result);
    }
    else if(x <= 1.0) {
      const double ln_term = -log(fabs(x));
      const double scale_factor = ( scale ? exp(x) : 1.0 );
      gsl_sf_result result_c;
      cheb_eval_e(&E12_cs, x, &result_c);
      result->val  = scale_factor * (ln_term - 0.6875 + x + result_c.val);
      result->err  = scale_factor * (result_c.err + GSL_DBL_EPSILON * fabs(ln_term));
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if(x <= 4.0) {
      const double s = 1.0/x * ( scale ? 1.0 : exp(-x) );
      gsl_sf_result result_c;
      cheb_eval_e(&AE13_cs, (8.0/x-5.0)/3.0, &result_c);
      result->val  = s * (1.0 + result_c.val);
      result->err  = s * result_c.err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if(x <= xmax || scale) {
      const double s = 1.0/x * ( scale ? 1.0 : exp(-x) );
      gsl_sf_result result_c;
      cheb_eval_e(&AE14_cs, 8.0/x-1.0, &result_c);
      result->val  = s * (1.0 +  result_c.val);
      result->err  = s * (GSL_DBL_EPSILON + result_c.err);
      result->err += 2.0 * (x + 1.0) * GSL_DBL_EPSILON * fabs(result->val);
      if(result->val == 0.0)
	UNDERFLOW_ERROR(result);
      else
	return GSL_SUCCESS;
    }
    else {
      UNDERFLOW_ERROR(result);
    }
  }

  static int expint_E2_impl(const double x, gsl_sf_result * result, const int scale) {
    const double xmaxt = -GSL_LOG_DBL_MIN;
    const double xmax  = xmaxt - log(xmaxt);
    
    /* CHECK_POINTER(result) */
    
    if(x < -xmax && !scale) {
      OVERFLOW_ERROR(result);
    }
    else if (x == 0.0) {
      result->val = (scale ? 1.0 : 1.0);
      result->err = 0.0;
      return GSL_SUCCESS;
    } else if(x < 100.0) {
      const double ex = ( scale ? 1.0 : exp(-x) );
      gsl_sf_result result_E1;
      int stat_E1 = expint_E1_impl(x, &result_E1, scale);
      result->val  = ex - x*result_E1.val;
      result->err  = GSL_DBL_EPSILON*ex + fabs(x) * result_E1.err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return stat_E1;
    }
    else if(x < xmax || scale) {
      const double s = ( scale ? 1.0 : exp(-x) );
      const double c1  = -2.0;
      const double c2  =  6.0;
      const double c3  = -24.0;
      const double c4  =  120.0;
      const double c5  = -720.0;
      const double c6  =  5040.0;
      const double c7  = -40320.0;
      const double c8  =  362880.0;
      const double c9  = -3628800.0;
      const double c10 =  39916800.0;
      const double c11 = -479001600.0;
      const double c12 =  6227020800.0;
      const double c13 = -87178291200.0;
      const double y = 1.0/x;
      const double sum6 = c6+y*(c7+y*(c8+y*(c9+y*(c10+y*(c11+y*(c12+y*c13))))));
      const double sum  = y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*sum6)))));
      result->val = s * (1.0 + sum)/x;
      result->err = 2.0 * (x + 1.0) * GSL_DBL_EPSILON * result->val;
      if(result->val == 0.0)
	UNDERFLOW_ERROR(result);
      else
	return GSL_SUCCESS;
    }
    else {
      UNDERFLOW_ERROR(result);
    }
  }

  /*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

  int gsl_sf_expint_E1_e(const double x, gsl_sf_result * result)
  {
    return expint_E1_impl(x, result, 0);
  }
  
  int gsl_sf_expint_E1_scaled_e(const double x, gsl_sf_result * result)
  {
    return expint_E1_impl(x, result, 1);
  }
  
  int gsl_sf_expint_E2_e(const double x, gsl_sf_result * result)
  {
    return expint_E2_impl(x, result, 0);
  }

  int gsl_sf_expint_E2_scaled_e(const double x, gsl_sf_result * result)
  {
    return expint_E2_impl(x, result, 1);
  }

  int gsl_sf_expint_Ei_e(const double x, gsl_sf_result * result)
  {
    /* CHECK_POINTER(result) */
    
    {
      int status = gsl_sf_expint_E1_e(-x, result);
      result->val = -result->val;
      return status;
    }
  }
  
  int gsl_sf_expint_Ei_scaled_e(const double x, gsl_sf_result * result)
  {
    /* CHECK_POINTER(result) */
    
    {
      int status = gsl_sf_expint_E1_scaled_e(-x, result);
      result->val = -result->val;
      return status;
    }
  }
  
#if 0
  static double recurse_En(int n, double x, double E1)
  {
    int i;
    double En = E1;
    double ex = exp(-x);
    for(i=2; i<=n; i++) {
      En = 1./(i-1) * (ex - x * En);
    }
    return En;
  }
#endif
  
  /*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
  
  double gsl_sf_expint_E1(const double x)
  {
    EVAL_RESULT(gsl_sf_expint_E1_e(x, &result));
  }
  
  double gsl_sf_expint_E1_scaled(const double x)
  {
    EVAL_RESULT(gsl_sf_expint_E1_scaled_e(x, &result));
  }
  
  double gsl_sf_expint_E2(const double x)
  {
    EVAL_RESULT(gsl_sf_expint_E2_e(x, &result));
  }
  
  double gsl_sf_expint_E2_scaled(const double x)
  {
    EVAL_RESULT(gsl_sf_expint_E2_scaled_e(x, &result));
  }
  
  double gsl_sf_expint_Ei(const double x)
  {
    EVAL_RESULT(gsl_sf_expint_Ei_e(x, &result));
  }
  
  double gsl_sf_expint_Ei_scaled(const double x)
  {
    EVAL_RESULT(gsl_sf_expint_Ei_scaled_e(x, &result));
  }

} // expint.c

// fdiv.c /////////////////////////////////////////////////////
extern "C" {
  double gsl_fdiv (const double x, const double y);

  double gsl_fdiv (const double x, const double y) {
    return x / y;
  }
  
} // fdiv.c

// gamma.c ///////////////////////////////////////////////////////////////////////
extern "C" { 
#define LogRootTwoPi_  0.9189385332046727418

  /*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

  static struct {int n; double f; long i; } fact_table[GSL_SF_FACT_NMAX + 1] = {
    { 0,  1.0,     1L     },
    { 1,  1.0,     1L     },
    { 2,  2.0,     2L     },
    { 3,  6.0,     6L     },
    { 4,  24.0,    24L    },
    { 5,  120.0,   120L   },
    { 6,  720.0,   720L   },
    { 7,  5040.0,  5040L  },
    { 8,  40320.0, 40320L },

    { 9,  362880.0,     362880L    },
    { 10, 3628800.0,    3628800L   },
    { 11, 39916800.0,   39916800L  },
    { 12, 479001600.0,  479001600L },

    { 13, 6227020800.0,                               0 },
    { 14, 87178291200.0,                              0 },
    { 15, 1307674368000.0,                            0 },
    { 16, 20922789888000.0,                           0 },
    { 17, 355687428096000.0,                          0 },
    { 18, 6402373705728000.0,                         0 },
    { 19, 121645100408832000.0,                       0 },
    { 20, 2432902008176640000.0,                      0 },
    { 21, 51090942171709440000.0,                     0 },
    { 22, 1124000727777607680000.0,                   0 },
    { 23, 25852016738884976640000.0,                  0 },
    { 24, 620448401733239439360000.0,                 0 },
    { 25, 15511210043330985984000000.0,               0 },
    { 26, 403291461126605635584000000.0,              0 },
    { 27, 10888869450418352160768000000.0,            0 },
    { 28, 304888344611713860501504000000.0,           0 },
    { 29, 8841761993739701954543616000000.0,          0 },
    { 30, 265252859812191058636308480000000.0,        0 },
    { 31, 8222838654177922817725562880000000.0,       0 },
    { 32, 263130836933693530167218012160000000.0,     0 },
    { 33, 8683317618811886495518194401280000000.0,    0 },
    { 34, 2.95232799039604140847618609644e38,  0 },
    { 35, 1.03331479663861449296666513375e40,  0 },
    { 36, 3.71993326789901217467999448151e41,  0 },
    { 37, 1.37637530912263450463159795816e43,  0 },
    { 38, 5.23022617466601111760007224100e44,  0 },
    { 39, 2.03978820811974433586402817399e46,  0 },
    { 40, 8.15915283247897734345611269600e47,  0 },
    { 41, 3.34525266131638071081700620534e49,  0 },
    { 42, 1.40500611775287989854314260624e51,  0 },
    { 43, 6.04152630633738356373551320685e52,  0 },
    { 44, 2.65827157478844876804362581101e54,  0 },
    { 45, 1.19622220865480194561963161496e56,  0 },
    { 46, 5.50262215981208894985030542880e57,  0 },
    { 47, 2.58623241511168180642964355154e59,  0 },
    { 48, 1.24139155925360726708622890474e61,  0 },
    { 49, 6.08281864034267560872252163321e62,  0 },
    { 50, 3.04140932017133780436126081661e64,  0 },
    { 51, 1.55111875328738228022424301647e66,  0 },
    { 52, 8.06581751709438785716606368564e67,  0 },
    { 53, 4.27488328406002556429801375339e69,  0 },
    { 54, 2.30843697339241380472092742683e71,  0 },
    { 55, 1.26964033536582759259651008476e73,  0 },
    { 56, 7.10998587804863451854045647464e74,  0 },
    { 57, 4.05269195048772167556806019054e76,  0 },
    { 58, 2.35056133128287857182947491052e78,  0 },
    { 59, 1.38683118545689835737939019720e80,  0 },
    { 60, 8.32098711274139014427634118320e81,  0 },
    { 61, 5.07580213877224798800856812177e83,  0 },
    { 62, 3.14699732603879375256531223550e85,  0 },
    { 63, 1.982608315404440064116146708360e87,  0 },
    { 64, 1.268869321858841641034333893350e89,  0 },
    { 65, 8.247650592082470666723170306800e90,  0 },
    { 66, 5.443449390774430640037292402480e92,  0 },
    { 67, 3.647111091818868528824985909660e94,  0 },
    { 68, 2.480035542436830599600990418570e96,  0 },
    { 69, 1.711224524281413113724683388810e98,  0 },
    { 70, 1.197857166996989179607278372170e100,  0 },
    { 71, 8.504785885678623175211676442400e101,  0 },
    { 72, 6.123445837688608686152407038530e103,  0 },
    { 73, 4.470115461512684340891257138130e105,  0 },
    { 74, 3.307885441519386412259530282210e107,  0 },
    { 75, 2.480914081139539809194647711660e109,  0 },
    { 76, 1.885494701666050254987932260860e111,  0 },
    { 77, 1.451830920282858696340707840860e113,  0 },
    { 78, 1.132428117820629783145752115870e115,  0 },
    { 79, 8.946182130782975286851441715400e116,  0 },
    { 80, 7.156945704626380229481153372320e118,  0 },
    { 81, 5.797126020747367985879734231580e120,  0 },
    { 82, 4.753643337012841748421382069890e122,  0 },
    { 83, 3.945523969720658651189747118010e124,  0 },
    { 84, 3.314240134565353266999387579130e126,  0 },
    { 85, 2.817104114380550276949479442260e128,  0 },
    { 86, 2.422709538367273238176552320340e130,  0 },
    { 87, 2.107757298379527717213600518700e132,  0 },
    { 88, 1.854826422573984391147968456460e134,  0 },
    { 89, 1.650795516090846108121691926250e136,  0 },
    { 90, 1.485715964481761497309522733620e138,  0 },
    { 91, 1.352001527678402962551665687590e140,  0 },
    { 92, 1.243841405464130725547532432590e142,  0 },
    { 93, 1.156772507081641574759205162310e144,  0 },
    { 94, 1.087366156656743080273652852570e146,  0 },
    { 95, 1.032997848823905926259970209940e148,  0 },
    { 96, 9.916779348709496892095714015400e149,  0 },
    { 97, 9.619275968248211985332842594960e151,  0 },
    { 98, 9.426890448883247745626185743100e153,  0 },
    { 99, 9.332621544394415268169923885600e155,  0 },
    { 100, 9.33262154439441526816992388563e157,  0 },
    { 101, 9.42594775983835942085162312450e159,  0 },
    { 102, 9.61446671503512660926865558700e161,  0 },
    { 103, 9.90290071648618040754671525458e163,  0 },
    { 104, 1.02990167451456276238485838648e166,  0 },
    { 105, 1.08139675824029090050410130580e168,  0 },
    { 106, 1.146280563734708354534347384148e170,  0 },
    { 107, 1.226520203196137939351751701040e172,  0 },
    { 108, 1.324641819451828974499891837120e174,  0 },
    { 109, 1.443859583202493582204882102460e176,  0 },
    { 110, 1.588245541522742940425370312710e178,  0 },
    { 111, 1.762952551090244663872161047110e180,  0 },
    { 112, 1.974506857221074023536820372760e182,  0 },
    { 113, 2.231192748659813646596607021220e184,  0 },
    { 114, 2.543559733472187557120132004190e186,  0 },
    { 115, 2.925093693493015690688151804820e188,  0 },
    { 116, 3.393108684451898201198256093590e190,  0 },
    { 117, 3.96993716080872089540195962950e192,  0 },
    { 118, 4.68452584975429065657431236281e194,  0 },
    { 119, 5.57458576120760588132343171174e196,  0 },
    { 120, 6.68950291344912705758811805409e198,  0 },
    { 121, 8.09429852527344373968162284545e200,  0 },
    { 122, 9.87504420083360136241157987140e202,  0 },
    { 123, 1.21463043670253296757662432419e205,  0 },
    { 124, 1.50614174151114087979501416199e207,  0 },
    { 125, 1.88267717688892609974376770249e209,  0 },
    { 126, 2.37217324288004688567714730514e211,  0 },
    { 127, 3.01266001845765954480997707753e213,  0 },
    { 128, 3.85620482362580421735677065923e215,  0 },
    { 129, 4.97450422247728744039023415041e217,  0 },
    { 130, 6.46685548922047367250730439554e219,  0 },
    { 131, 8.47158069087882051098456875820e221,  0 },
    { 132, 1.11824865119600430744996307608e224,  0 },
    { 133, 1.48727070609068572890845089118e226,  0 },
    { 134, 1.99294274616151887673732419418e228,  0 },
    { 135, 2.69047270731805048359538766215e230,  0 },
    { 136, 3.65904288195254865768972722052e232,  0 },
    { 137, 5.01288874827499166103492629211e234,  0 },
    { 138, 6.91778647261948849222819828311e236,  0 },
    { 139, 9.61572319694108900419719561353e238,  0 },
    { 140, 1.34620124757175246058760738589e241,  0 },
    { 141, 1.89814375907617096942852641411e243,  0 },
    { 142, 2.69536413788816277658850750804e245,  0 },
    { 143, 3.85437071718007277052156573649e247,  0 },
    { 144, 5.55029383273930478955105466055e249,  0 },
    { 145, 8.04792605747199194484902925780e251,  0 },
    { 146, 1.17499720439091082394795827164e254,  0 },
    { 147, 1.72724589045463891120349865931e256,  0 },
    { 148, 2.55632391787286558858117801578e258,  0 },
    { 149, 3.80892263763056972698595524351e260,  0 },
    { 150, 5.71338395644585459047893286526e262,  0 },
    { 151, 8.62720977423324043162318862650e264,  0 },
    { 152, 1.31133588568345254560672467123e267,  0 },
    { 153, 2.00634390509568239477828874699e269,  0 },
    { 154, 3.08976961384735088795856467036e271,  0 },
    { 155, 4.78914290146339387633577523906e273,  0 },
    { 156, 7.47106292628289444708380937294e275,  0 },
    { 157, 1.17295687942641442819215807155e278,  0 },
    { 158, 1.85327186949373479654360975305e280,  0 },
    { 159, 2.94670227249503832650433950735e282,  0 },
    { 160, 4.71472363599206132240694321176e284,  0 },
    { 161, 7.59070505394721872907517857094e286,  0 },
    { 162, 1.22969421873944943411017892849e289,  0 },
    { 163, 2.00440157654530257759959165344e291,  0 },
    { 164, 3.28721858553429622726333031164e293,  0 },
    { 165, 5.42391066613158877498449501421e295,  0 },
    { 166, 9.00369170577843736647426172359e297,  0 },
    { 167, 1.50361651486499904020120170784e300,  0 },
    { 168, 2.52607574497319838753801886917e302,  0 },
    { 169, 4.26906800900470527493925188890e304,  0 },
    { 170, 7.25741561530799896739672821113e306,  0 },

    /*
    { 171, 1.24101807021766782342484052410e309,  0 },
    { 172, 2.13455108077438865629072570146e311,  0 },
    { 173, 3.69277336973969237538295546352e313,  0 },
    { 174, 6.42542566334706473316634250653e315,  0 },
    { 175, 1.12444949108573632830410993864e318,  0 },
    { 176, 1.97903110431089593781523349201e320,  0 },
    { 177, 3.50288505463028580993296328086e322,  0 },
    { 178, 6.23513539724190874168067463993e324,  0 },
    { 179, 1.11608923610630166476084076055e327,  0 },
    { 180, 2.00896062499134299656951336898e329,  0 },
    { 181, 3.63621873123433082379081919786e331,  0 },
    { 182, 6.61791809084648209929929094011e333,  0 },
    { 183, 1.21107901062490622417177024204e336,  0 },
    { 184, 2.22838537954982745247605724535e338,  0 },
    { 185, 4.12251295216718078708070590390e340,  0 },
    { 186, 7.66787409103095626397011298130e342,  0 },
    { 187, 1.43389245502278882136241112750e345,  0 },
    { 188, 2.69571781544284298416133291969e347,  0 },
    { 189, 5.09490667118697324006491921822e349,  0 },
    { 190, 9.68032267525524915612334651460e351,  0 },
    { 191, 1.84894163097375258881955918429e354,  0 },
    { 192, 3.54996793146960497053355363384e356,  0 },
    { 193, 6.85143810773633759312975851330e358,  0 },
    { 194, 1.32917899290084949306717315158e361,  0 },
    { 195, 2.59189903615665651148098764559e363,  0 },
    { 196, 5.08012211086704676250273578535e365,  0 },
    { 197, 1.00078405584080821221303894971e368,  0 },
    { 198, 1.98155243056480026018181712043e370,  0 },
    { 199, 3.94328933682395251776181606966e372,  0 },
    { 200, 7.88657867364790503552363213932e374,  0 }
    */
  };

  static struct {int n; double f; long i; } doub_fact_table[GSL_SF_DOUBLEFACT_NMAX + 1] = {
    { 0,  1.000000000000000000000000000,    1L    },
    { 1,  1.000000000000000000000000000,    1L    },
    { 2,  2.000000000000000000000000000,    2L    },
    { 3,  3.000000000000000000000000000,    3L    },
    { 4,  8.000000000000000000000000000,    8L    },
    { 5,  15.00000000000000000000000000,    15L   },
    { 6,  48.00000000000000000000000000,    48L   },
    { 7,  105.0000000000000000000000000,    105L  },
    { 8,  384.0000000000000000000000000,    384L  },
    { 9,  945.0000000000000000000000000,    945L  },
    { 10, 3840.000000000000000000000000,    3840L   },
    { 11, 10395.00000000000000000000000,    10395L  },
    { 12, 46080.00000000000000000000000,       46080L       },
    { 13, 135135.0000000000000000000000,       135135L      },
    { 14, 645120.00000000000000000000000,      645120L      },
    { 15, 2.02702500000000000000000000000e6,   2027025L     },
    { 16, 1.03219200000000000000000000000e7,   10321920L    },
    { 17, 3.4459425000000000000000000000e7,    34459425L    },
    { 18, 1.85794560000000000000000000000e8,   185794560L   },
    { 19, 6.5472907500000000000000000000e8,            0 },
    { 20, 3.7158912000000000000000000000e9,            0 },
    { 21, 1.37493105750000000000000000000e10,          0 },
    { 22, 8.1749606400000000000000000000e10,           0 },
    { 23, 3.1623414322500000000000000000e11,           0 },
    { 24, 1.96199055360000000000000000000e12,          0 },
    { 25, 7.9058535806250000000000000000e12,           0 },
    { 26, 5.1011754393600000000000000000e13,           0 },
    { 27, 2.13458046676875000000000000000e14,          0 },
    { 28, 1.42832912302080000000000000000e15,          0 },
    { 29, 6.1902833536293750000000000000e15,           0 },
    { 30, 4.2849873690624000000000000000e16,           0 },
    { 31, 1.91898783962510625000000000000e17,          0 },
    { 32, 1.37119595809996800000000000000e18,          0 },
    { 33, 6.3326598707628506250000000000e18,           0 },
    { 34, 4.6620662575398912000000000000e19,           0 },
    { 35, 2.21643095476699771875000000000e20,          0 },
    { 36, 1.67834385271436083200000000000e21,          0 },
    { 37, 8.2007945326378915593750000000e21,           0 },
    { 38, 6.3777066403145711616000000000e22,           0 },
    { 39, 3.1983098677287777081562500000e23,           0 },
    { 40, 2.55108265612582846464000000000e24,          0 },
    { 41, 1.31130704576879886034406250000e25,          0 },
    { 42, 1.07145471557284795514880000000e26,          0 },
    { 43, 5.6386202968058350994794687500e26,           0 },
    { 44, 4.7144007485205310026547200000e27,           0 },
    { 45, 2.53737913356262579476576093750e28,          0 },
    { 46, 2.16862434431944426122117120000e29,          0 },
    { 47, 1.19256819277443412353990764062e30,          0 },
    { 48, 1.04093968527333324538616217600e31,          0 },
    { 49, 5.8435841445947272053455474391e31,           0 },
    { 50, 5.2046984263666662269308108800e32,           0 },
    { 51, 2.98022791374331087472622919392e33,          0 },
    { 52, 2.70644318171066643800402165760e34,          0 },
    { 53, 1.57952079428395476360490147278e35,          0 },
    { 54, 1.46147931812375987652217169510e36,          0 },
    { 55, 8.6873643685617511998269581003e36,           0 },
    { 56, 8.1842841814930553085241614926e37,           0 },
    { 57, 4.9517976900801981839013661172e38,           0 },
    { 58, 4.7468848252659720789440136657e39,           0 },
    { 59, 2.92156063714731692850180600912e40,       0 },
    { 60, 2.84813089515958324736640819942e41,       0 },
    { 61, 1.78215198865986332638610166557e42,       0 },
    { 62, 1.76584115499894161336717308364e43,       0 },
    { 63, 1.12275575285571389562324404931e44,       0 },
    { 64, 1.13013833919932263255499077353e45,       0 },
    { 65, 7.2979123935621403215510863205e45,        0 },
    { 66, 7.4589130387155293748629391053e46,        0 },
    { 67, 4.8896013036866340154392278347e47,        0 },
    { 68, 5.0720608663265599749067985916e48,        0 },
    { 69, 3.3738248995437774706530672060e49,        0 },
    { 70, 3.5504426064285919824347590141e50,        0 },
    { 71, 2.39541567867608200416367771623e51,       0 },
    { 72, 2.55631867662858622735302649017e52,       0 },
    { 73, 1.74865344543353986303948473285e53,       0 },
    { 74, 1.89167582070515380824123960272e54,       0 },
    { 75, 1.31149008407515489727961354964e55,       0 },
    { 76, 1.43767362373591689426334209807e56,       0 },
    { 77, 1.00984736473786927090530243322e57,       0 },
    { 78, 1.12138542651401517752540683649e58,       0 },
    { 79, 7.9777941814291672401518892225e58,        0 },
    { 80, 8.9710834121121214202032546920e59,        0 },
    { 81, 6.4620132869576254645230302702e60,        0 },
    { 82, 7.3562883979319395645666688474e61,        0 },
    { 83, 5.3634710281748291355541151243e62,        0 },
    { 84, 6.1792822542628292342360018318e63,        0 },
    { 85, 4.5589503739486047652209978556e64,        0 },
    { 86, 5.3141827386660331414429615754e65,        0 },
    { 87, 3.9662868253352861457422681344e66,        0 },
    { 88, 4.6764808100261091644698061863e67,        0 },
    { 89, 3.5299952745484046697106186396e68,        0 },
    { 90, 4.2088327290234982480228255677e69,        0 },
    { 91, 3.2122956998390482494366629620e70,        0 },
    { 92, 3.8721261107016183881809995223e71,        0 },
    { 93, 2.98743500085031487197609655470e72,       0 },
    { 94, 3.6397985440595212848901395509e73,        0 },
    { 95, 2.83806325080779912837729172696e74,       0 },
    { 96, 3.4942066022971404334945339689e75,        0 },
    { 97, 2.75292135328356515452597297515e76,       0 },
    { 98, 3.4243224702511976248246432895e77,        0 },
    { 99, 2.72539213975072950298071324540e78,       0 },
    { 100, 3.4243224702511976248246432895e79,       0 },
    { 101, 2.75264606114823679801052037785e80,      0 },
    { 102, 3.4928089196562215773211361553e81,       0 },
    { 103, 2.83522544298268390195083598919e82,      0 },
    { 104, 3.6325212764424704404139816015e83,       0 },
    { 105, 2.97698671513181809704837778865e84,      0 },
    { 106, 3.8504725530290186668388204976e85,       0 },
    { 107, 3.1853757851910453638417642339e86,       0 },
    { 108, 4.1585103572713401601859261374e87,       0 },
    { 109, 3.4720596058582394465875230149e88,       0 },
    { 110, 4.5743613929984741762045187512e89,       0 },
    { 111, 3.8539861625026457857121505465e90,       0 },
    { 112, 5.1232847601582910773490610013e91,       0 },
    { 113, 4.3550043636279897378547301176e92,       0 },
    { 114, 5.8405446265804518281779295415e93,       0 },
    { 115, 5.0082550181721881985329396352e94,       0 },
    { 116, 6.7750317668333241206863982681e95,       0 },
    { 117, 5.8596583712614601922835393732e96,       0 },
    { 118, 7.9945374848633224624099499564e97,       0 },
    { 119, 6.9729934618011376288174118541e98,       0 },
    { 120, 9.5934449818359869548919399477e99,       0 },
    { 121, 8.4373220887793765308690683435e100,      0 },
    { 122, 1.17040028778399040849681667362e102,       0 },
    { 123, 1.03779061691986331329689540625e103,       0 },
    { 124, 1.45129635685214810653605267528e104,       0 },
    { 125, 1.29723827114982914162111925781e105,       0 },
    { 126, 1.82863340963370661423542637086e106,       0 },
    { 127, 1.64749260436028300985882145742e107,       0 },
    { 128, 2.34065076433114446622134575470e108,       0 },
    { 129, 2.12526545962476508271787968008e109,       0 },
    { 130, 3.04284599363048780608774948111e110,       0 },
    { 131, 2.78409775210844225836042238090e111,       0 },
    { 132, 4.0165567115922439040358293151e112,        0 },
    { 133, 3.7028500103042282036193617666e113,        0 },
    { 134, 5.3821859935336068314080112822e114,        0 },
    { 135, 4.9988475139107080748861383849e115,        0 },
    { 136, 7.3197729512057052907148953438e116,        0 },
    { 137, 6.8484210940576700625940095873e117,        0 },
    { 138, 1.01012866726638733011865555744e119,       0 },
    { 139, 9.5193053207401613870056733264e119,        0 },
    { 140, 1.41418013417294226216611778042e121,       0 },
    { 141, 1.34222205022436275556779993902e122,       0 },
    { 142, 2.00813579052557801227588724819e123,       0 },
    { 143, 1.91937753182083874046195391280e124,       0 },
    { 144, 2.89171553835683233767727763739e125,       0 },
    { 145, 2.78309742114021617366983317355e126,       0 },
    { 146, 4.2219046860009752130088253506e127,        0 },
    { 147, 4.0911532090761177752946547651e128,        0 },
    { 148, 6.2484189352814433152530615189e129,        0 },
    { 149, 6.0958182815234154851890356000e130,        0 },
    { 150, 9.3726284029221649728795922783e131,        0 },
    { 151, 9.2046856051003573826354437561e132,        0 },
    { 152, 1.42463951724416907587769802630e134,       0 },
    { 153, 1.40831689758035467954322289468e135,       0 },
    { 154, 2.19394485655602037685165496051e136,       0 },
    { 155, 2.18289119124954975329199548675e137,       0 },
    { 156, 3.4225539762273917878885817384e138,        0 },
    { 157, 3.4271391702617931126684329142e139,        0 },
    { 158, 5.4076352824392790248639591467e140,        0 },
    { 159, 5.4491512807162510491428083336e141,        0 },
    { 160, 8.6522164519028464397823346347e142,        0 },
    { 161, 8.7731335619531641891199214170e143,        0 },
    { 162, 1.40165906520826112324473821082e145,       0 },
    { 163, 1.43002077059836576282654719098e146,       0 },
    { 164, 2.29872086694154824212137066574e147,       0 },
    { 165, 2.35953427148730350866380286512e148,       0 },
    { 166, 3.8158766391229700819214753051e149,        0 },
    { 167, 3.9404222333837968594685507847e150,        0 },
    { 168, 6.4106727537265897376280785126e151,        0 },
    { 169, 6.6593135744186166925018508262e152,        0 },
    { 170, 1.08981436813352025539677334714e154,       0 },
    { 171, 1.13874262122558345441781649128e155,       0 },
    { 172, 1.87448071318965483928245015709e156,       0 },
    { 173, 1.97002473472025937614282252992e157,       0 },
    { 174, 3.2615964409499994203514632733e158,        0 },
    { 175, 3.4475432857604539082499394274e159,        0 },
    { 176, 5.7404097360719989798185753611e160,        0 },
    { 177, 6.1021516157960034176023927864e161,        0 },
    { 178, 1.02179293302081581840770641427e163,       0 },
    { 179, 1.09228513922748461175082830877e164,       0 },
    { 180, 1.83922727943746847313387154568e165,       0 },
    { 181, 1.97703610200174714726899923887e166,       0 },
    { 182, 3.3473936485761926211036462131e167,        0 },
    { 183, 3.6179760666631972795022686071e168,        0 },
    { 184, 6.1592043133801944228307090322e169,        0 },
    { 185, 6.6932557233269149670791969232e170,        0 },
    { 186, 1.14561200228871616264651187999e172,       0 },
    { 187, 1.25163882026213309884380982464e173,       0 },
    { 188, 2.15375056430278638577544233437e174,       0 },
    { 189, 2.36559737029543155681480056857e175,       0 },
    { 190, 4.0921260721752941329733404353e176,        0 },
    { 191, 4.5182909772642742735162690860e177,        0 },
    { 192, 7.8568820585765647353088136358e178,        0 },
    { 193, 8.7203015861200493478863993359e179,        0 },
    { 194, 1.52423511936385355864990984535e181,       0 },
    { 195, 1.70045880929340962283784787050e182,       0 },
    { 196, 2.98750083395315297495382329688e183,       0 },
    { 197, 3.3499038543080169569905603049e184,        0 },
    { 198, 5.9152516512272428904085701278e185,        0 },
    { 199, 6.6663086700729537444112150067e186,        0 },
    { 200, 1.18305033024544857808171402556e188,       0 },
    { 201, 1.33992804268466370262665421635e189,       0 },
    { 202, 2.38976166709580612772506233164e190,       0 },
    { 203, 2.72005392664986731633210805920e191,       0 },
    { 204, 4.8751138008754445005591271565e192,        0 },
    { 205, 5.5761105496322279984808215214e193,        0 },
    { 206, 1.00427344298034156711518019425e195,       0 },
    { 207, 1.15425488377387119568553005492e196,       0 },
    { 208, 2.08888876139911045959957480403e197,       0 },
    { 209, 2.41239270708739079898275781478e198,       0 },
    { 210, 4.3866663989381319651591070885e199,        0 },
    { 211, 5.0901486119543945858536189892e200,        0 },
    { 212, 9.2997327657488397661373070276e201,        0 },
    { 213, 1.08420165434628604678682084470e203,       0 },
    { 214, 1.99014281187025170995338370390e204,       0 },
    { 215, 2.33103355684451500059166481610e205,       0 },
    { 216, 4.2987084736397436934993088004e206,        0 },
    { 217, 5.0583428183525975512839126509e207,        0 },
    { 218, 9.3711844725346412518284931849e208,        0 },
    { 219, 1.10777707721921886373117687056e210,       0 },
    { 220, 2.06166058395762107540226850068e211,       0 },
    { 221, 2.44818734065447368884590088393e212,       0 },
    { 222, 4.5768864963859187873930360715e213,        0 },
    { 223, 5.4594577696594763261263589712e214,        0 },
    { 224, 1.02522257519044580837604008002e216,       0 },
    { 225, 1.22837799817338217337843076851e217,       0 },
    { 226, 2.31700301993040752692985058084e218,       0 },
    { 227, 2.78841805585357753356903784452e219,       0 },
    { 228, 5.2827668854413291614000593243e220,        0 },
    { 229, 6.3854773479046925518730966640e221,        0 },
    { 230, 1.21503638365150570712201364459e223,       0 },
    { 231, 1.47504526736598397948268532937e224,       0 },
    { 232, 2.81888441007149324052307165546e225,       0 },
    { 233, 3.4368554729627426721946568174e226,        0 },
    { 234, 6.5961895195672941828239876738e227,        0 },
    { 235, 8.0766103614624452796574435210e228,        0 },
    { 236, 1.55670072661788142714646109101e230,       0 },
    { 237, 1.91415665566659953127881411447e231,       0 },
    { 238, 3.7049477293505577966085773966e232,        0 },
    { 239, 4.5748344070431728797563657336e233,        0 },
    { 240, 8.8918745504413387118605857518e234,        0 },
    { 241, 1.10253509209740466402128414180e236,       0 },
    { 242, 2.15183364120680396827026175195e237,       0 },
    { 243, 2.67916027379669333357172046456e238,       0 },
    { 244, 5.2504740845446016825794386748e239,        0 },
    { 245, 6.5639426708018986672507151382e240,        0 },
    { 246, 1.29161662479797201391454191399e242,       0 },
    { 247, 1.62129383968806897081092663913e243,       0 },
    { 248, 3.2032092294989705945080639467e244,        0 },
    { 249, 4.0370216608232917373192073314e245,        0 },
    { 250, 8.0080230737474264862701598667e246,        0 },
    { 251, 1.01329243686664622606712104019e248,       0 },
    { 252, 2.01802181458435147454008028642e249,       0 },
    { 253, 2.56362986527261495194981623168e250,       0 },
    { 254, 5.1257754090442527453318039275e251,        0 },
    { 255, 6.5372561564451681274720313908e252,        0 },
    { 256, 1.31219850471532870280494180544e254,       0 },
    { 257, 1.68007483220640820876031206743e255,       0 },
    { 258, 3.3854721421655480532367498580e256,        0 },
    { 259, 4.3513938154145972606892082546e257,        0 },
    { 260, 8.8022275696304249384155496309e258,        0 },
    { 261, 1.13571378582320988503988335446e260,       0 },
    { 262, 2.30618362324317133386487400329e261,       0 },
    { 263, 2.98692725671504199765489322224e262,       0 },
    { 264, 6.0883247653619723214032673687e263,        0 },
    { 265, 7.9153572302948612937854670389e264,        0 },
    { 266, 1.61949438758628463749326912007e266,       0 },
    { 267, 2.11340038048872796544071969939e267,       0 },
    { 268, 4.3402449587312428284819612418e268,        0 },
    { 269, 5.6850470235146782270355359914e269,        0 },
    { 270, 1.17186613885743556369012953528e271,       0 },
    { 271, 1.54064774337247779952663025366e272,       0 },
    { 272, 3.1874758976922247332371523360e273,        0 },
    { 273, 4.2059683394068643927077005925e274,        0 },
    { 274, 8.7336839596766957690697974006e275,        0 },
    { 275, 1.15664129333688770799461766294e277,       0 },
    { 276, 2.41049677287076803226326408256e278,       0 },
    { 277, 3.2038963825431789511450909263e279,        0 },
    { 278, 6.7011810285807351296918741495e280,        0 },
    { 279, 8.9388709072954692736948036845e281,        0 },
    { 280, 1.87633068800260583631372476186e283,       0 },
    { 281, 2.51182272495002686590823983534e284,       0 },
    { 282, 5.2912525401673484584047038284e285,        0 },
    { 283, 7.1084583116085760305203187340e286,        0 },
    { 284, 1.50271572140752696218693588728e288,       0 },
    { 285, 2.02591061880844416869829083919e289,       0 },
    { 286, 4.2977669632255271118546366376e290,        0 },
    { 287, 5.8143634759802347641640947085e291,        0 },
    { 288, 1.23775688540895180821413535163e293,       0 },
    { 289, 1.68035104455828784684342337075e294,       0 },
    { 290, 3.5894949676859602438209925197e295,        0 },
    { 291, 4.8898215396646176343143620089e296,        0 },
    { 292, 1.04813253056430039119572981576e298,       0 },
    { 293, 1.43271771112173296685410806860e299,       0 },
    { 294, 3.08150963985904315011544565835e300,       0 },
    { 295, 4.2265172478091122522196188024e301,        0 },
    { 296, 9.1212685339827677243417191487e302,        0 },
    { 297, 1.25527562259930633890922678431e304,       0 },
    /*
      { 298, 2.71813802312686478185383230631e305,       0 },
      { 299, 3.7532741115719259533385880851e306,        0 },
      { 300, 8.1544140693805943455614969189e307,  }
    */
  };


/* Chebyshev coefficients for Gamma*(3/4(t+1)+1/2), -1<t<1
 */
static double gstar_a_data[30] = {
  2.16786447866463034423060819465,
 -0.05533249018745584258035832802,
  0.01800392431460719960888319748,
 -0.00580919269468937714480019814,
  0.00186523689488400339978881560,
 -0.00059746524113955531852595159,
  0.00019125169907783353925426722,
 -0.00006124996546944685735909697,
  0.00001963889633130842586440945,
 -6.3067741254637180272515795142e-06,
  2.0288698405861392526872789863e-06,
 -6.5384896660838465981983750582e-07,
  2.1108698058908865476480734911e-07,
 -6.8260714912274941677892994580e-08,
  2.2108560875880560555583978510e-08,
 -7.1710331930255456643627187187e-09,
  2.3290892983985406754602564745e-09,
 -7.5740371598505586754890405359e-10,
  2.4658267222594334398525312084e-10,
 -8.0362243171659883803428749516e-11,
  2.6215616826341594653521346229e-11,
 -8.5596155025948750540420068109e-12,
  2.7970831499487963614315315444e-12,
 -9.1471771211886202805502562414e-13,
  2.9934720198063397094916415927e-13,
 -9.8026575909753445931073620469e-14,
  3.2116773667767153777571410671e-14,
 -1.0518035333878147029650507254e-14,
  3.4144405720185253938994854173e-15,
 -1.0115153943081187052322643819e-15
};
static cheb_series gstar_a_cs = {
  gstar_a_data,
  29,
  -1, 1,
  17
};


/* Chebyshev coefficients for
 * x^2(Gamma*(x) - 1 - 1/(12x)), x = 4(t+1)+2, -1 < t < 1
 */
static double gstar_b_data[] = {
  0.0057502277273114339831606096782,
  0.0004496689534965685038254147807,
 -0.0001672763153188717308905047405,
  0.0000615137014913154794776670946,
 -0.0000223726551711525016380862195,
  8.0507405356647954540694800545e-06,
 -2.8671077107583395569766746448e-06,
  1.0106727053742747568362254106e-06,
 -3.5265558477595061262310873482e-07,
  1.2179216046419401193247254591e-07,
 -4.1619640180795366971160162267e-08,
  1.4066283500795206892487241294e-08,
 -4.6982570380537099016106141654e-09,
  1.5491248664620612686423108936e-09,
 -5.0340936319394885789686867772e-10,
  1.6084448673736032249959475006e-10,
 -5.0349733196835456497619787559e-11,
  1.5357154939762136997591808461e-11,
 -4.5233809655775649997667176224e-12,
  1.2664429179254447281068538964e-12,
 -3.2648287937449326771785041692e-13,
  7.1528272726086133795579071407e-14,
 -9.4831735252566034505739531258e-15,
 -2.3124001991413207293120906691e-15,
  2.8406613277170391482590129474e-15,
 -1.7245370321618816421281770927e-15,
  8.6507923128671112154695006592e-16,
 -3.9506563665427555895391869919e-16,
  1.6779342132074761078792361165e-16,
 -6.0483153034414765129837716260e-17
};
static cheb_series gstar_b_cs = {
  gstar_b_data,
  29,
  -1, 1,
  18
};


/* coefficients for gamma=7, kmax=8  Lanczos method */
static double lanczos_7_c[9] = {
  0.99999999999980993227684700473478,
  676.520368121885098567009190444019,
 -1259.13921672240287047156078755283,
  771.3234287776530788486528258894,
 -176.61502916214059906584551354,
  12.507343278686904814458936853,
 -0.13857109526572011689554707,
  9.984369578019570859563e-6,
  1.50563273514931155834e-7
};

/* complex version of Lanczos method; this is not safe for export
 * since it becomes bad in the left half-plane
 */
static
int lngamma_lanczos_complex(double zr, double zi, gsl_sf_result * yr, gsl_sf_result * yi)
{
  int k;
  gsl_sf_result log1_r,    log1_i;
  gsl_sf_result logAg_r,   logAg_i;
  double Ag_r, Ag_i;
  double yi_tmp_val, yi_tmp_err;

  zr -= 1.0; /* Lanczos writes z! instead of Gamma(z) */

  Ag_r = lanczos_7_c[0];
  Ag_i = 0.0;
  for(k=1; k<=8; k++) {
    double R = zr + k;
    double I = zi;
    double a = lanczos_7_c[k] / (R*R + I*I);
    Ag_r +=  a * R;
    Ag_i -=  a * I;
  }

  gsl_sf_complex_log_e(zr + 7.5, zi, &log1_r,  &log1_i);
  gsl_sf_complex_log_e(Ag_r, Ag_i,   &logAg_r, &logAg_i);

  /* (z+0.5)*log(z+7.5) - (z+7.5) + LogRootTwoPi_ + log(Ag(z)) */
  yr->val = (zr+0.5)*log1_r.val - zi*log1_i.val - (zr+7.5) + LogRootTwoPi_ + logAg_r.val;
  yi->val = zi*log1_r.val + (zr+0.5)*log1_i.val - zi + logAg_i.val;
  yr->err = 4.0 * GSL_DBL_EPSILON * fabs(yr->val);
  yi->err = 4.0 * GSL_DBL_EPSILON * fabs(yi->val);
  yi_tmp_val = yi->val;
  yi_tmp_err = yi->err;
  gsl_sf_angle_restrict_symm_err_e(yi_tmp_val, yi);
  yi->err += yi_tmp_err;
  return GSL_SUCCESS;
}


/* Lanczos method for real x > 0;
 * gamma=7, truncated at 1/(z+8) 
 * [J. SIAM Numer. Anal, Ser. B, 1 (1964) 86]
 */
static
int
lngamma_lanczos(double x, gsl_sf_result * result)
{
  int k;
  double Ag;
  double term1, term2;

  x -= 1.0; /* Lanczos writes z! instead of Gamma(z) */

  Ag = lanczos_7_c[0];
  for(k=1; k<=8; k++) { Ag += lanczos_7_c[k]/(x+k); }

  /* (x+0.5)*log(x+7.5) - (x+7.5) + LogRootTwoPi_ + log(Ag(x)) */
  term1 = (x+0.5)*log((x+7.5)/M_E);
  term2 = LogRootTwoPi_ + log(Ag);
  result->val  = term1 + (term2 - 7.0);
  result->err  = 2.0 * GSL_DBL_EPSILON * (fabs(term1) + fabs(term2) + 7.0);
  result->err += GSL_DBL_EPSILON * fabs(result->val);

  return GSL_SUCCESS;
}

/* x = eps near zero
 * gives double-precision for |eps| < 0.02
 */
static
int
lngamma_sgn_0(double eps, gsl_sf_result * lng, double * sgn)
{
  /* calculate series for g(eps) = Gamma(eps) eps - 1/(1+eps) - eps/2 */
  const double c1  = -0.07721566490153286061;
  const double c2  = -0.01094400467202744461;
  const double c3  =  0.09252092391911371098;
  const double c4  = -0.01827191316559981266;
  const double c5  =  0.01800493109685479790;
  const double c6  = -0.00685088537872380685;
  const double c7  =  0.00399823955756846603;
  const double c8  = -0.00189430621687107802;
  const double c9  =  0.00097473237804513221;
  const double c10 = -0.00048434392722255893;
  const double g6  = c6+eps*(c7+eps*(c8 + eps*(c9 + eps*c10)));
  const double g   = eps*(c1+eps*(c2+eps*(c3+eps*(c4+eps*(c5+eps*g6)))));

  /* calculate Gamma(eps) eps, a positive quantity */
  const double gee = g + 1.0/(1.0+eps) + 0.5*eps;

  lng->val = log(gee/fabs(eps));
  lng->err = 4.0 * GSL_DBL_EPSILON * fabs(lng->val);
  *sgn = GSL_SIGN(eps);

  return GSL_SUCCESS;
}


/* x near a negative integer
 * Calculates sign as well as log(|gamma(x)|).
 * x = -N + eps
 * assumes N >= 1
 */
static
int
lngamma_sgn_sing(int N, double eps, gsl_sf_result * lng, double * sgn)
{
  if(eps == 0.0) {
    lng->val = 0.0;
    lng->err = 0.0;
    *sgn = 0.0;
    GSL_ERROR ("error", GSL_EDOM);
  }
  else if(N == 1) {
    /* calculate series for
     * g = eps gamma(-1+eps) + 1 + eps/2 (1+3eps)/(1-eps^2)
     * double-precision for |eps| < 0.02
     */
    const double c0 =  0.07721566490153286061;
    const double c1 =  0.08815966957356030521;
    const double c2 = -0.00436125434555340577;
    const double c3 =  0.01391065882004640689;
    const double c4 = -0.00409427227680839100;
    const double c5 =  0.00275661310191541584;
    const double c6 = -0.00124162645565305019;
    const double c7 =  0.00065267976121802783;
    const double c8 = -0.00032205261682710437;
    const double c9 =  0.00016229131039545456;
    const double g5 = c5 + eps*(c6 + eps*(c7 + eps*(c8 + eps*c9)));
    const double g  = eps*(c0 + eps*(c1 + eps*(c2 + eps*(c3 + eps*(c4 + eps*g5)))));

    /* calculate eps gamma(-1+eps), a negative quantity */
    const double gam_e = g - 1.0 - 0.5*eps*(1.0+3.0*eps)/(1.0 - eps*eps);

    lng->val = log(fabs(gam_e)/fabs(eps));
    lng->err = 2.0 * GSL_DBL_EPSILON * fabs(lng->val);
    *sgn = ( eps > 0.0 ? -1.0 : 1.0 );
    return GSL_SUCCESS;
  }
  else {
    double g;

    /* series for sin(Pi(N+1-eps))/(Pi eps) modulo the sign
     * double-precision for |eps| < 0.02
     */
    const double cs1 = -1.6449340668482264365;
    const double cs2 =  0.8117424252833536436;
    const double cs3 = -0.1907518241220842137;
    const double cs4 =  0.0261478478176548005;
    const double cs5 = -0.0023460810354558236;
    const double e2  = eps*eps;
    const double sin_ser = 1.0 + e2*(cs1+e2*(cs2+e2*(cs3+e2*(cs4+e2*cs5))));

    /* calculate series for ln(gamma(1+N-eps))
     * double-precision for |eps| < 0.02
     */
    double aeps = fabs(eps);
    double c1, c2, c3, c4, c5, c6, c7;
    double lng_ser;
    gsl_sf_result c0;
    gsl_sf_result psi_0;
    gsl_sf_result psi_1;
    gsl_sf_result psi_2;
    gsl_sf_result psi_3;
    gsl_sf_result psi_4;
    gsl_sf_result psi_5;
    gsl_sf_result psi_6;
    psi_2.val = 0.0;
    psi_3.val = 0.0;
    psi_4.val = 0.0;
    psi_5.val = 0.0;
    psi_6.val = 0.0;
    gsl_sf_lnfact_e(N, &c0);
    gsl_sf_psi_int_e(N+1, &psi_0);
    gsl_sf_psi_1_int_e(N+1, &psi_1);
    if(aeps > 0.00001) gsl_sf_psi_n_e(2, N+1.0, &psi_2);
    if(aeps > 0.0002)  gsl_sf_psi_n_e(3, N+1.0, &psi_3);
    if(aeps > 0.001)   gsl_sf_psi_n_e(4, N+1.0, &psi_4);
    if(aeps > 0.005)   gsl_sf_psi_n_e(5, N+1.0, &psi_5);
    if(aeps > 0.01)    gsl_sf_psi_n_e(6, N+1.0, &psi_6);
    c1 = psi_0.val;
    c2 = psi_1.val/2.0;
    c3 = psi_2.val/6.0;
    c4 = psi_3.val/24.0;
    c5 = psi_4.val/120.0;
    c6 = psi_5.val/720.0;
    c7 = psi_6.val/5040.0;
    lng_ser = c0.val-eps*(c1-eps*(c2-eps*(c3-eps*(c4-eps*(c5-eps*(c6-eps*c7))))));

    /* calculate
     * g = ln(|eps gamma(-N+eps)|)
     *   = -ln(gamma(1+N-eps)) + ln(|eps Pi/sin(Pi(N+1+eps))|)
     */
    g = -lng_ser - log(sin_ser);

    lng->val = g - log(fabs(eps));
    lng->err = c0.err + 2.0 * GSL_DBL_EPSILON * (fabs(g) + fabs(lng->val));

    *sgn = ( GSL_IS_ODD(N) ? -1.0 : 1.0 ) * ( eps > 0.0 ? 1.0 : -1.0 );

    return GSL_SUCCESS;
  }
}


/* This gets bad near the negative half axis. However, this
 * region can be avoided by use of the reflection formula, as usual.
 * Only the first two terms of the series are kept.
 */
#if 0
static
int
lngamma_complex_stirling(const double zr, const double zi, double * lg_r, double * arg)
{
  double re_zinv,  im_zinv;
  double re_zinv2, im_zinv2;
  double re_zinv3, im_zinv3;
  double re_zhlnz, im_zhlnz;
  double r, lnr, theta;
  gsl_sf_complex_log_e(zr, zi, &lnr, &theta);  /* z = r e^{i theta} */
  r = exp(lnr);
  re_zinv =  (zr/r)/r;
  im_zinv = -(zi/r)/r;
  re_zinv2 = re_zinv*re_zinv - im_zinv*im_zinv;
  re_zinv2 = 2.0*re_zinv*im_zinv;
  re_zinv3 = re_zinv2*re_zinv - im_zinv2*im_zinv;
  re_zinv3 = re_zinv2*im_zinv + im_zinv2*re_zinv;
  re_zhlnz = (zr - 0.5)*lnr - zi*theta;
  im_zhlnz = zi*lnr + zr*theta;
  *lg_r = re_zhlnz - zr + 0.5*(M_LN2+M_LNPI) + re_zinv/12.0 - re_zinv3/360.0;
  *arg  = im_zhlnz - zi + 1.0/12.0*im_zinv - im_zinv3/360.0;
  return GSL_SUCCESS;
}
#endif /* 0 */


inline
static
int
lngamma_1_pade(const double eps, gsl_sf_result * result)
{
  /* Use (2,2) Pade for Log[Gamma[1+eps]]/eps
   * plus a correction series.
   */
  const double n1 = -1.0017419282349508699871138440;
  const double n2 =  1.7364839209922879823280541733;
  const double d1 =  1.2433006018858751556055436011;
  const double d2 =  5.0456274100274010152489597514;
  const double num = (eps + n1) * (eps + n2);
  const double den = (eps + d1) * (eps + d2);
  const double pade = 2.0816265188662692474880210318 * num / den;
  const double c0 =  0.004785324257581753;
  const double c1 = -0.01192457083645441;
  const double c2 =  0.01931961413960498;
  const double c3 = -0.02594027398725020;
  const double c4 =  0.03141928755021455;
  const double eps5 = eps*eps*eps*eps*eps;
  const double corr = eps5 * (c0 + eps*(c1 + eps*(c2 + eps*(c3 + c4*eps))));
  result->val = eps * (pade + corr);
  result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  return GSL_SUCCESS;
}

inline
static
int
lngamma_2_pade(const double eps, gsl_sf_result * result)
{
  /* Use (2,2) Pade for Log[Gamma[2+eps]]/eps
   * plus a correction series.
   */
  const double n1 = 1.000895834786669227164446568;
  const double n2 = 4.209376735287755081642901277;
  const double d1 = 2.618851904903217274682578255;
  const double d2 = 10.85766559900983515322922936;
  const double num = (eps + n1) * (eps + n2);
  const double den = (eps + d1) * (eps + d2);
  const double pade = 2.85337998765781918463568869 * num/den;
  const double c0 =  0.0001139406357036744;
  const double c1 = -0.0001365435269792533;
  const double c2 =  0.0001067287169183665;
  const double c3 = -0.0000693271800931282;
  const double c4 =  0.0000407220927867950;
  const double eps5 = eps*eps*eps*eps*eps;
  const double corr = eps5 * (c0 + eps*(c1 + eps*(c2 + eps*(c3 + c4*eps))));
  result->val = eps * (pade + corr);
  result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  return GSL_SUCCESS;
}


/* series for gammastar(x)
 * double-precision for x > 10.0
 */
static
int
gammastar_ser(const double x, gsl_sf_result * result)
{
  /* Use the Stirling series for the correction to Log(Gamma(x)),
   * which is better behaved and easier to compute than the
   * regular Stirling series for Gamma(x). 
   */
  const double y = 1.0/(x*x);
  const double c0 =  1.0/12.0;
  const double c1 = -1.0/360.0;
  const double c2 =  1.0/1260.0;
  const double c3 = -1.0/1680.0;
  const double c4 =  1.0/1188.0;
  const double c5 = -691.0/360360.0;
  const double c6 =  1.0/156.0;
  const double c7 = -3617.0/122400.0;
  const double ser = c0 + y*(c1 + y*(c2 + y*(c3 + y*(c4 + y*(c5 + y*(c6 + y*c7))))));
  result->val = exp(ser/x);
  result->err = 2.0 * GSL_DBL_EPSILON * result->val * GSL_MAX_DBL(1.0, ser/x);
  return GSL_SUCCESS;
}


/* Chebyshev expansion for log(gamma(x)/gamma(8))
 * 5 < x < 10
 * -1 < t < 1
 */
static double gamma_5_10_data[24] = {
 -1.5285594096661578881275075214,
  4.8259152300595906319768555035,
  0.2277712320977614992970601978,
 -0.0138867665685617873604917300,
  0.0012704876495201082588139723,
 -0.0001393841240254993658962470,
  0.0000169709242992322702260663,
 -2.2108528820210580075775889168e-06,
  3.0196602854202309805163918716e-07,
 -4.2705675000079118380587357358e-08,
  6.2026423818051402794663551945e-09,
 -9.1993973208880910416311405656e-10,
  1.3875551258028145778301211638e-10,
 -2.1218861491906788718519522978e-11,
  3.2821736040381439555133562600e-12,
 -5.1260001009953791220611135264e-13,
  8.0713532554874636696982146610e-14,
 -1.2798522376569209083811628061e-14,
  2.0417711600852502310258808643e-15,
 -3.2745239502992355776882614137e-16,
  5.2759418422036579482120897453e-17,
 -8.5354147151695233960425725513e-18,
  1.3858639703888078291599886143e-18,
 -2.2574398807738626571560124396e-19
};
static const cheb_series gamma_5_10_cs = {
  gamma_5_10_data,
  23,
  -1, 1,
  11
};


/* gamma(x) for x >= 1/2
 * assumes x >= 1/2
 */
static
int
gamma_xgthalf(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x == 0.5) {
    result->val = 1.77245385090551602729817;
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  } else if (x <= (GSL_SF_FACT_NMAX + 1.0) && x == floor(x)) {
    int n = (int) floor (x);
    result->val = fact_table[n - 1].f;
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }    
  else if(fabs(x - 1.0) < 0.01) {
    /* Use series for Gamma[1+eps] - 1/(1+eps).
     */
    const double eps = x - 1.0;
    const double c1 =  0.4227843350984671394;
    const double c2 = -0.01094400467202744461;
    const double c3 =  0.09252092391911371098;
    const double c4 = -0.018271913165599812664;
    const double c5 =  0.018004931096854797895;
    const double c6 = -0.006850885378723806846;
    const double c7 =  0.003998239557568466030;
    result->val = 1.0/x + eps*(c1+eps*(c2+eps*(c3+eps*(c4+eps*(c5+eps*(c6+eps*c7))))));
    result->err = GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
  else if(fabs(x - 2.0) < 0.01) {
    /* Use series for Gamma[1 + eps].
     */
    const double eps = x - 2.0;
    const double c1 =  0.4227843350984671394;
    const double c2 =  0.4118403304264396948;
    const double c3 =  0.08157691924708626638;
    const double c4 =  0.07424901075351389832;
    const double c5 = -0.00026698206874501476832;
    const double c6 =  0.011154045718130991049;
    const double c7 = -0.002852645821155340816;
    const double c8 =  0.0021039333406973880085;
    result->val = 1.0 + eps*(c1+eps*(c2+eps*(c3+eps*(c4+eps*(c5+eps*(c6+eps*(c7+eps*c8)))))));
    result->err = GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
  else if(x < 5.0) {
    /* Exponentiating the logarithm is fine, as
     * long as the exponential is not so large
     * that it greatly amplifies the error.
     */
    gsl_sf_result lg;
    lngamma_lanczos(x, &lg);
    result->val = exp(lg.val);
    result->err = result->val * (lg.err + 2.0 * GSL_DBL_EPSILON);
    return GSL_SUCCESS;
  }
  else if(x < 10.0) {
    /* This is a sticky area. The logarithm
     * is too large and the gammastar series
     * is not good.
     */
    const double gamma_8 = 5040.0;
    const double t = (2.0*x - 15.0)/5.0;
    gsl_sf_result c;
    cheb_eval_e(&gamma_5_10_cs, t, &c);
    result->val  = exp(c.val) * gamma_8;
    result->err  = result->val * c.err;
    result->err += 2.0 * GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
  else if(x < GSL_SF_GAMMA_XMAX) {
    /* We do not want to exponentiate the logarithm
     * if x is large because of the inevitable
     * inflation of the error. So we carefully
     * use pow() and exp() with exact quantities.
     */
    double p = pow(x, 0.5*x);
    double e = exp(-x);
    double q = (p * e) * p;
    double pre = M_SQRT2 * M_SQRTPI * q/sqrt(x);
    gsl_sf_result gstar;
    int stat_gs = gammastar_ser(x, &gstar);
    result->val = pre * gstar.val;
    result->err = (x + 2.5) * GSL_DBL_EPSILON * result->val;
    return stat_gs;
  }
  else {
    OVERFLOW_ERROR(result);
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/


int gsl_sf_lngamma_e(double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(fabs(x - 1.0) < 0.01) {
    /* Note that we must amplify the errors
     * from the Pade evaluations because of
     * the way we must pass the argument, i.e.
     * writing (1-x) is a loss of precision
     * when x is near 1.
     */
    int stat = lngamma_1_pade(x - 1.0, result);
    result->err *= 1.0/(GSL_DBL_EPSILON + fabs(x - 1.0));
    return stat;
  }
  else if(fabs(x - 2.0) < 0.01) {
    int stat = lngamma_2_pade(x - 2.0, result);
    result->err *= 1.0/(GSL_DBL_EPSILON + fabs(x - 2.0));
    return stat;
  }
  else if(x >= 0.5) {
    return lngamma_lanczos(x, result);
  }
  else if(x == 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(fabs(x) < 0.02) {
    double sgn;
    return lngamma_sgn_0(x, result, &sgn);
  }
  else if(x > -0.5/(GSL_DBL_EPSILON*M_PI)) {
    /* Try to extract a fractional
     * part from x.
     */
    double z  = 1.0 - x;
    double s  = sin(M_PI*z);
    double as = fabs(s);
    if(s == 0.0) {
      DOMAIN_ERROR(result);
    }
    else if(as < M_PI*0.015) {
      /* x is near a negative integer, -N */
      if(x < INT_MIN + 2.0) {
        result->val = 0.0;
        result->err = 0.0;
        GSL_ERROR ("error", GSL_EROUND);
      }
      else {
        int N = -(int)(x - 0.5);
        double eps = x + N;
        double sgn;
        return lngamma_sgn_sing(N, eps, result, &sgn);
      }
    }
    else {
      gsl_sf_result lg_z;
      lngamma_lanczos(z, &lg_z);
      result->val = M_LNPI - (log(as) + lg_z.val);
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val) + lg_z.err;
      return GSL_SUCCESS;
    }
  }
  else {
    /* |x| was too large to extract any fractional part */
    result->val = 0.0;
    result->err = 0.0;
    GSL_ERROR ("error", GSL_EROUND);
  }
}


int gsl_sf_lngamma_sgn_e(double x, gsl_sf_result * result_lg, double * sgn)
{
  if(fabs(x - 1.0) < 0.01) {
    int stat = lngamma_1_pade(x - 1.0, result_lg);
    result_lg->err *= 1.0/(GSL_DBL_EPSILON + fabs(x - 1.0));
    *sgn = 1.0;
    return stat;
  }
  else if(fabs(x - 2.0) < 0.01) {
   int stat = lngamma_2_pade(x - 2.0, result_lg);
    result_lg->err *= 1.0/(GSL_DBL_EPSILON + fabs(x - 2.0));
    *sgn = 1.0;
    return stat;
  }
  else if(x >= 0.5) {
    *sgn = 1.0;
    return lngamma_lanczos(x, result_lg);
  }
  else if(x == 0.0) {
    *sgn = 0.0;
    DOMAIN_ERROR(result_lg);
  }
  else if(fabs(x) < 0.02) {
    return lngamma_sgn_0(x, result_lg, sgn);
  }
  else if(x > -0.5/(GSL_DBL_EPSILON*M_PI)) {
   /* Try to extract a fractional
     * part from x.
     */
    double z = 1.0 - x;
    double s = sin(M_PI*x);
    double as = fabs(s);
    if(s == 0.0) {
      *sgn = 0.0;
      DOMAIN_ERROR(result_lg);
    }
    else if(as < M_PI*0.015) {
      /* x is near a negative integer, -N */
      if(x < INT_MIN + 2.0) {
        result_lg->val = 0.0;
        result_lg->err = 0.0;
        *sgn = 0.0;
        GSL_ERROR ("error", GSL_EROUND);
      }
      else {
        int N = -(int)(x - 0.5);
        double eps = x + N;
        return lngamma_sgn_sing(N, eps, result_lg, sgn);
      }
    }
    else {
      gsl_sf_result lg_z;
      lngamma_lanczos(z, &lg_z);
      *sgn = (s > 0.0 ? 1.0 : -1.0);
      result_lg->val = M_LNPI - (log(as) + lg_z.val);
      result_lg->err = 2.0 * GSL_DBL_EPSILON * fabs(result_lg->val) + lg_z.err;
      return GSL_SUCCESS;
    }
  }
  else {
    /* |x| was too large to extract any fractional part */
    result_lg->val = 0.0;
    result_lg->err = 0.0;
    *sgn = 0.0;
    GSL_ERROR ("error", GSL_EROUND);
  }
}


int
gsl_sf_gamma_e(const double x, gsl_sf_result * result)
{
  if(x < 0.5) {
    int rint_x = (int)floor(x+0.5);
    double f_x = x - rint_x;
    double sgn_gamma = ( GSL_IS_EVEN(rint_x) ? 1.0 : -1.0 );
    double sin_term = sgn_gamma * sin(M_PI * f_x) / M_PI;

    if(sin_term == 0.0) {
      DOMAIN_ERROR(result);
    }
    else if(x > -169.0) {
      gsl_sf_result g;
      gamma_xgthalf(1.0-x, &g);
      if(fabs(sin_term) * g.val * GSL_DBL_MIN < 1.0) {
        result->val  = 1.0/(sin_term * g.val);
        result->err  = fabs(g.err/g.val) * fabs(result->val);
        result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
        return GSL_SUCCESS;
      }
      else {
        UNDERFLOW_ERROR(result);
      }
    }
    else {
      /* It is hard to control it here.
       * We can only exponentiate the
       * logarithm and eat the loss of
       * precision.
       */
      gsl_sf_result lng;
      double sgn;
      int stat_lng = gsl_sf_lngamma_sgn_e(x, &lng, &sgn);
      int stat_e   = gsl_sf_exp_mult_err_e(lng.val, lng.err, sgn, 0.0, result);
      return GSL_ERROR_SELECT_2(stat_e, stat_lng);
    }
  }
  else {
    return gamma_xgthalf(x, result);
  }
}


int
gsl_sf_gammastar_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x < 0.5) {
    gsl_sf_result lg;
    const int stat_lg = gsl_sf_lngamma_e(x, &lg);
    const double lx = log(x);
    const double c  = 0.5*(M_LN2+M_LNPI);
    const double lnr_val = lg.val - (x-0.5)*lx + x - c;
    const double lnr_err = lg.err + 2.0 * GSL_DBL_EPSILON *((x+0.5)*fabs(lx) + c);
    const int stat_e  = gsl_sf_exp_err_e(lnr_val, lnr_err, result);
    return GSL_ERROR_SELECT_2(stat_lg, stat_e);
  }
  else if(x < 2.0) {
    const double t = 4.0/3.0*(x-0.5) - 1.0;
    return cheb_eval_e(&gstar_a_cs, t, result);
  }
  else if(x < 10.0) {
    const double t = 0.25*(x-2.0) - 1.0;
    gsl_sf_result c;
    cheb_eval_e(&gstar_b_cs, t, &c);
    result->val  = c.val/(x*x) + 1.0 + 1.0/(12.0*x);
    result->err  = c.err/(x*x);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < 1.0/GSL_ROOT4_DBL_EPSILON) {
    return gammastar_ser(x, result);
  }
  else if(x < 1.0/GSL_DBL_EPSILON) {
    /* Use Stirling formula for Gamma(x).
     */
    const double xi = 1.0/x;
    result->val = 1.0 + xi/12.0*(1.0 + xi/24.0*(1.0 - xi*(139.0/180.0 + 571.0/8640.0*xi)));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = 1.0;
    result->err = 1.0/x;
    return GSL_SUCCESS;
  }
}


int
gsl_sf_gammainv_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if (x <= 0.0 && x == floor(x)) { /* negative integer */
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  } else if(x < 0.5) {
    gsl_sf_result lng;
    double sgn;
    int stat_lng = gsl_sf_lngamma_sgn_e(x, &lng, &sgn);
    if(stat_lng == GSL_EDOM) {
      result->val = 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
    else if(stat_lng != GSL_SUCCESS) {
      result->val = 0.0;
      result->err = 0.0;
      return stat_lng;
    }
    else {
      return gsl_sf_exp_mult_err_e(-lng.val, lng.err, sgn, 0.0, result);
    }
  }
  else {
    gsl_sf_result g;
    int stat_g = gamma_xgthalf(x, &g);
    if(stat_g == GSL_EOVRFLW) {
      UNDERFLOW_ERROR(result);
    }
    else {
      result->val  = 1.0/g.val;
      result->err  = fabs(g.err/g.val) * fabs(result->val);
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      CHECK_UNDERFLOW(result);
      return GSL_SUCCESS;
    }
  }
}


int
gsl_sf_lngamma_complex_e(double zr, double zi, gsl_sf_result * lnr, gsl_sf_result * arg)
{
  if(zr <= 0.5) {
    /* Transform to right half plane using reflection;
     * in fact we do a little better by stopping at 1/2.
     */
    double x = 1.0-zr;
    double y = -zi;
    gsl_sf_result a, b;
    gsl_sf_result lnsin_r, lnsin_i;

    int stat_l = lngamma_lanczos_complex(x, y, &a, &b);
    int stat_s = gsl_sf_complex_logsin_e(M_PI*zr, M_PI*zi, &lnsin_r, &lnsin_i);

    if(stat_s == GSL_SUCCESS) {
      int stat_r;
      lnr->val = M_LNPI - lnsin_r.val - a.val;
      lnr->err = lnsin_r.err + a.err + 2.0 * GSL_DBL_EPSILON * fabs(lnr->val);
      arg->val = -lnsin_i.val - b.val;
      arg->err = lnsin_i.err + b.err + 2.0 * GSL_DBL_EPSILON * fabs(arg->val);
      stat_r = gsl_sf_angle_restrict_symm_e(&(arg->val));
      return GSL_ERROR_SELECT_2(stat_r, stat_l);
    }
    else {
      DOMAIN_ERROR_2(lnr,arg);
    }
  }
  else {
    /* otherwise plain vanilla Lanczos */
    return lngamma_lanczos_complex(zr, zi, lnr, arg);
  }
}


int gsl_sf_taylorcoeff_e(const int n, const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x < 0.0 || n < 0) {
    DOMAIN_ERROR(result);
  }
  else if(n == 0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(n == 1) {
    result->val = x;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(x == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    const double log2pi = M_LNPI + M_LN2;
    const double ln_test = n*(log(x)+1.0) + 1.0 - (n+0.5)*log(n+1.0) + 0.5*log2pi;

    if(ln_test < GSL_LOG_DBL_MIN+1.0) {
      UNDERFLOW_ERROR(result);
    }
    else if(ln_test > GSL_LOG_DBL_MAX-1.0) {
      OVERFLOW_ERROR(result);
    }
    else {
      double product = 1.0;
      int k;
      for(k=1; k<=n; k++) {
        product *= (x/k);
      }
      result->val = product;
      result->err = n * GSL_DBL_EPSILON * product;
      CHECK_UNDERFLOW(result);
      return GSL_SUCCESS;
    }    
  }
}


int gsl_sf_fact_e(const unsigned int n, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(n < 18) {
    result->val = fact_table[n].f;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(n <= GSL_SF_FACT_NMAX){
    result->val = fact_table[n].f;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    OVERFLOW_ERROR(result);
  }
}


int gsl_sf_doublefact_e(const unsigned int n, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(n < 26) {
    result->val = doub_fact_table[n].f;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(n <= GSL_SF_DOUBLEFACT_NMAX){
    result->val = doub_fact_table[n].f;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    OVERFLOW_ERROR(result);
  }
}


int gsl_sf_lnfact_e(const unsigned int n, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(n <= GSL_SF_FACT_NMAX){
    result->val = log(fact_table[n].f);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_lngamma_e(n+1.0, result);
    return GSL_SUCCESS;
  }
}


int gsl_sf_lndoublefact_e(const unsigned int n, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(n <= GSL_SF_DOUBLEFACT_NMAX){
    result->val = log(doub_fact_table[n].f);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(GSL_IS_ODD(n)) {
    gsl_sf_result lg;
    gsl_sf_lngamma_e(0.5*(n+2.0), &lg);
    result->val = 0.5*(n+1.0) * M_LN2 - 0.5*M_LNPI + lg.val;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val) + lg.err;
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result lg;
    gsl_sf_lngamma_e(0.5*n+1.0, &lg);
    result->val = 0.5*n*M_LN2 + lg.val;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val) + lg.err;
    return GSL_SUCCESS;
  }
}


int gsl_sf_lnchoose_e(unsigned int n, unsigned int m, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(m > n) {
    DOMAIN_ERROR(result);
  }
  else if(m == n || m == 0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result nf;
    gsl_sf_result mf;
    gsl_sf_result nmmf;
    if(m*2 > n) m = n-m;
    gsl_sf_lnfact_e(n, &nf);
    gsl_sf_lnfact_e(m, &mf);
    gsl_sf_lnfact_e(n-m, &nmmf);
    result->val  = nf.val - mf.val - nmmf.val;
    result->err  = nf.err + mf.err + nmmf.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int gsl_sf_choose_e(unsigned int n, unsigned int m, gsl_sf_result * result)
{
  if(m > n) {
    DOMAIN_ERROR(result);
  }
  else if(m == n || m == 0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if (n <= GSL_SF_FACT_NMAX) {
    result->val = (fact_table[n].f / fact_table[m].f) / fact_table[n-m].f;
    result->err = 6.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  } else {
    if(m*2 < n) m = n-m;

    if (n - m < 64)  /* compute product for a manageable number of terms */
      {
        double prod = 1.0;
        unsigned int k;
        
        for(k=n; k>=m+1; k--) {
          double tk = (double)k / (double)(k-m);
          if(tk > GSL_DBL_MAX/prod) {
            OVERFLOW_ERROR(result);
          }
          prod *= tk;
        }
        result->val = prod;
        result->err = 2.0 * GSL_DBL_EPSILON * prod * fabs(n-m);
        return GSL_SUCCESS;
      }
    else
      {
        gsl_sf_result lc;
        const int stat_lc = gsl_sf_lnchoose_e (n, m, &lc);
        const int stat_e  = gsl_sf_exp_err_e(lc.val, lc.err, result);
        return GSL_ERROR_SELECT_2(stat_lc, stat_e);
      }
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_fact(const unsigned int n)
{
  EVAL_RESULT(gsl_sf_fact_e(n, &result));
}

double gsl_sf_lnfact(const unsigned int n)
{
  EVAL_RESULT(gsl_sf_lnfact_e(n, &result));
}

double gsl_sf_doublefact(const unsigned int n)
{
  EVAL_RESULT(gsl_sf_doublefact_e(n, &result));
}

double gsl_sf_lndoublefact(const unsigned int n)
{
  EVAL_RESULT(gsl_sf_lndoublefact_e(n, &result));
}

double gsl_sf_lngamma(const double x)
{
  EVAL_RESULT(gsl_sf_lngamma_e(x, &result));
}

double gsl_sf_gamma(const double x)
{
  EVAL_RESULT(gsl_sf_gamma_e(x, &result));
}

double gsl_sf_gammastar(const double x)
{
  EVAL_RESULT(gsl_sf_gammastar_e(x, &result));
}

double gsl_sf_gammainv(const double x)
{
  EVAL_RESULT(gsl_sf_gammainv_e(x, &result));
}

double gsl_sf_taylorcoeff(const int n, const double x)
{
  EVAL_RESULT(gsl_sf_taylorcoeff_e(n, x, &result));
}

double gsl_sf_choose(unsigned int n, unsigned int m)
{
  EVAL_RESULT(gsl_sf_choose_e(n, m, &result));
}

double gsl_sf_lnchoose(unsigned int n, unsigned int m)
{
  EVAL_RESULT(gsl_sf_lnchoose_e(n, m, &result));
}
} // gamma.c

// gamma_inc.c ///////////////////////////////////////////////////////////////////////
extern "C" { 
/* The dominant part,
 * D(a,x) := x^a e^(-x) / Gamma(a+1)
 */
static
int
gamma_inc_D(const double a, const double x, gsl_sf_result * result)
{
  if(a < 10.0) {
    double lnr;
    gsl_sf_result lg;
    gsl_sf_lngamma_e(a+1.0, &lg);
    lnr = a * log(x) - x - lg.val;
    result->val = exp(lnr);
    result->err = 2.0 * GSL_DBL_EPSILON * (fabs(lnr) + 1.0) * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result gstar;
    gsl_sf_result ln_term;
    double term1;
    if (x < 0.5*a) {
      double u = x/a;   
      double ln_u = log(u);
      ln_term.val = ln_u - u + 1.0;
      ln_term.err = (fabs(ln_u) + fabs(u) + 1.0) * GSL_DBL_EPSILON;
    } else {
      double mu = (x-a)/a;
      gsl_sf_log_1plusx_mx_e(mu, &ln_term);  /* log(1+mu) - mu */
    };
    gsl_sf_gammastar_e(a, &gstar);
    term1 = exp(a*ln_term.val)/sqrt(2.0*M_PI*a);
    result->val  = term1/gstar.val;
    result->err  = 2.0 * GSL_DBL_EPSILON * (fabs(a*ln_term.val) + 1.0) * fabs(result->val);
    result->err += gstar.err/fabs(gstar.val) * fabs(result->val);
    return GSL_SUCCESS;
  }

}


/* P series representation.
 */
static
int
gamma_inc_P_series(const double a, const double x, gsl_sf_result * result)
{
  const int nmax = 5000;

  gsl_sf_result D;
  int stat_D = gamma_inc_D(a, x, &D);

  double sum  = 1.0;
  double term = 1.0;
  int n;
  for(n=1; n<nmax; n++) {
    term *= x/(a+n);
    sum  += term;
    if(fabs(term/sum) < GSL_DBL_EPSILON) break;
  }

  result->val  = D.val * sum;
  result->err  = D.err * fabs(sum);
  result->err += (1.0 + n) * GSL_DBL_EPSILON * fabs(result->val);

  if(n == nmax)
    GSL_ERROR ("error", GSL_EMAXITER);
  else
    return stat_D;
}


/* Q large x asymptotic
 */
static
int
gamma_inc_Q_large_x(const double a, const double x, gsl_sf_result * result)
{
  const int nmax = 5000;

  gsl_sf_result D;
  const int stat_D = gamma_inc_D(a, x, &D);

  double sum  = 1.0;
  double term = 1.0;
  double last = 1.0;
  int n;
  for(n=1; n<nmax; n++) {
    term *= (a-n)/x;
    if(fabs(term/last) > 1.0) break;
    if(fabs(term/sum)  < GSL_DBL_EPSILON) break;
    sum  += term;
    last  = term;
  }

  result->val  = D.val * (a/x) * sum;
  result->err  = D.err * fabs((a/x) * sum);
  result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);

  if(n == nmax)
    GSL_ERROR ("error in large x asymptotic", GSL_EMAXITER);
  else
    return stat_D;
}


/* Uniform asymptotic for x near a, a and x large.
 * See [Temme, p. 285]
 */
static
int
gamma_inc_Q_asymp_unif(const double a, const double x, gsl_sf_result * result)
{
  const double rta = sqrt(a);
  const double eps = (x-a)/a;

  gsl_sf_result ln_term;
  const int stat_ln = gsl_sf_log_1plusx_mx_e(eps, &ln_term);  /* log(1+eps) - eps */
  const double eta  = GSL_SIGN(eps) * sqrt(-2.0*ln_term.val);

  gsl_sf_result erfc;

  double R;
  double c0, c1;

  /* This used to say erfc(eta*M_SQRT2*rta), which is wrong.
   * The sqrt(2) is in the denominator. Oops.
   * Fixed: [GJ] Mon Nov 15 13:25:32 MST 2004
   */
  gsl_sf_erfc_e(eta*rta/M_SQRT2, &erfc);

  if(fabs(eps) < GSL_ROOT5_DBL_EPSILON) {
    c0 = -1.0/3.0 + eps*(1.0/12.0 - eps*(23.0/540.0 - eps*(353.0/12960.0 - eps*589.0/30240.0)));
    c1 = -1.0/540.0 - eps/288.0;
  }
  else {
    const double rt_term = sqrt(-2.0 * ln_term.val/(eps*eps));
    const double lam = x/a;
    c0 = (1.0 - 1.0/rt_term)/eps;
    c1 = -(eta*eta*eta * (lam*lam + 10.0*lam + 1.0) - 12.0 * eps*eps*eps) / (12.0 * eta*eta*eta*eps*eps*eps);
  }

  R = exp(-0.5*a*eta*eta)/(M_SQRT2*M_SQRTPI*rta) * (c0 + c1/a);

  result->val  = 0.5 * erfc.val + R;
  result->err  = GSL_DBL_EPSILON * fabs(R * 0.5 * a*eta*eta) + 0.5 * erfc.err;
  result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);

  return stat_ln;
}


/* Continued fraction which occurs in evaluation
 * of Q(a,x) or Gamma(a,x).
 *
 *              1   (1-a)/x  1/x  (2-a)/x   2/x  (3-a)/x
 *   F(a,x) =  ---- ------- ----- -------- ----- -------- ...
 *             1 +   1 +     1 +   1 +      1 +   1 +
 *
 * Hans E. Plesser, 2002-01-22 (hans dot plesser at itf dot nlh dot no).
 *
 * Split out from gamma_inc_Q_CF() by GJ [Tue Apr  1 13:16:41 MST 2003].
 * See gamma_inc_Q_CF() below.
 *
 */
static int
gamma_inc_F_CF(const double a, const double x, gsl_sf_result * result)
{
  const int    nmax  =  5000;
  const double small =  gsl_pow_3 (GSL_DBL_EPSILON);

  double hn = 1.0;           /* convergent */
  double Cn = 1.0 / small;
  double Dn = 1.0;
  int n;

  /* n == 1 has a_1, b_1, b_0 independent of a,x,
     so that has been done by hand                */
  for ( n = 2 ; n < nmax ; n++ )
  {
    double an;
    double delta;

    if(GSL_IS_ODD(n))
      an = 0.5*(n-1)/x;
    else
      an = (0.5*n-a)/x;

    Dn = 1.0 + an * Dn;
    if ( fabs(Dn) < small )
      Dn = small;
    Cn = 1.0 + an/Cn;
    if ( fabs(Cn) < small )
      Cn = small;
    Dn = 1.0 / Dn;
    delta = Cn * Dn;
    hn *= delta;
    if(fabs(delta-1.0) < GSL_DBL_EPSILON) break;
  }

  result->val = hn;
  result->err = 2.0*GSL_DBL_EPSILON * fabs(hn);
  result->err += GSL_DBL_EPSILON * (2.0 + 0.5*n) * fabs(result->val);

  if(n == nmax)
    GSL_ERROR ("error in CF for F(a,x)", GSL_EMAXITER);
  else
    return GSL_SUCCESS;
}


/* Continued fraction for Q.
 *
 * Q(a,x) = D(a,x) a/x F(a,x)
 *
 * Hans E. Plesser, 2002-01-22 (hans dot plesser at itf dot nlh dot no):
 *
 * Since the Gautschi equivalent series method for CF evaluation may lead
 * to singularities, I have replaced it with the modified Lentz algorithm
 * given in
 *
 * I J Thompson and A R Barnett
 * Coulomb and Bessel Functions of Complex Arguments and Order
 * J Computational Physics 64:490-509 (1986)
 *
 * In consequence, gamma_inc_Q_CF_protected() is now obsolete and has been
 * removed.
 *
 * Identification of terms between the above equation for F(a, x) and
 * the first equation in the appendix of Thompson&Barnett is as follows:
 *
 *    b_0 = 0, b_n = 1 for all n > 0
 *
 *    a_1 = 1
 *    a_n = (n/2-a)/x    for n even
 *    a_n = (n-1)/(2x)   for n odd
 *
 */
static
int
gamma_inc_Q_CF(const double a, const double x, gsl_sf_result * result)
{
  gsl_sf_result D;
  gsl_sf_result F;
  const int stat_D = gamma_inc_D(a, x, &D);
  const int stat_F = gamma_inc_F_CF(a, x, &F);

  result->val  = D.val * (a/x) * F.val;
  result->err  = D.err * fabs((a/x) * F.val) + fabs(D.val * a/x * F.err);

  return GSL_ERROR_SELECT_2(stat_F, stat_D);
}


/* Useful for small a and x. Handles the subtraction analytically.
 */
static
int
gamma_inc_Q_series(const double a, const double x, gsl_sf_result * result)
{
  double term1;  /* 1 - x^a/Gamma(a+1) */
  double sum;    /* 1 + (a+1)/(a+2)(-x)/2! + (a+1)/(a+3)(-x)^2/3! + ... */
  int stat_sum;
  double term2;  /* a temporary variable used at the end */

  {
    /* Evaluate series for 1 - x^a/Gamma(a+1), small a
     */
    const double pg21 = -2.404113806319188570799476;  /* PolyGamma[2,1] */
    const double lnx  = log(x);
    const double el   = M_EULER+lnx;
    const double c1 = -el;
    const double c2 = M_PI*M_PI/12.0 - 0.5*el*el;
    const double c3 = el*(M_PI*M_PI/12.0 - el*el/6.0) + pg21/6.0;
    const double c4 = -0.04166666666666666667
                       * (-1.758243446661483480 + lnx)
                       * (-0.764428657272716373 + lnx)
                       * ( 0.723980571623507657 + lnx)
                       * ( 4.107554191916823640 + lnx);
    const double c5 = -0.0083333333333333333
                       * (-2.06563396085715900 + lnx)
                       * (-1.28459889470864700 + lnx)
                       * (-0.27583535756454143 + lnx)
                       * ( 1.33677371336239618 + lnx)
                       * ( 5.17537282427561550 + lnx);
    const double c6 = -0.0013888888888888889
                       * (-2.30814336454783200 + lnx)
                       * (-1.65846557706987300 + lnx)
                       * (-0.88768082560020400 + lnx)
                       * ( 0.17043847751371778 + lnx)
                       * ( 1.92135970115863890 + lnx)
                       * ( 6.22578557795474900 + lnx);
    const double c7 = -0.00019841269841269841
                       * (-2.5078657901291800 + lnx)
                       * (-1.9478900888958200 + lnx)
                       * (-1.3194837322612730 + lnx)
                       * (-0.5281322700249279 + lnx)
                       * ( 0.5913834939078759 + lnx)
                       * ( 2.4876819633378140 + lnx)
                       * ( 7.2648160783762400 + lnx);
    const double c8 = -0.00002480158730158730
                       * (-2.677341544966400 + lnx)
                       * (-2.182810448271700 + lnx)
                       * (-1.649350342277400 + lnx)
                       * (-1.014099048290790 + lnx)
                       * (-0.191366955370652 + lnx)
                       * ( 0.995403817918724 + lnx)
                       * ( 3.041323283529310 + lnx)
                       * ( 8.295966556941250 + lnx);
    const double c9 = -2.75573192239859e-6
                       * (-2.8243487670469080 + lnx)
                       * (-2.3798494322701120 + lnx)
                       * (-1.9143674728689960 + lnx)
                       * (-1.3814529102920370 + lnx)
                       * (-0.7294312810261694 + lnx)
                       * ( 0.1299079285269565 + lnx)
                       * ( 1.3873333251885240 + lnx)
                       * ( 3.5857258865210760 + lnx)
                       * ( 9.3214237073814600 + lnx);
    const double c10 = -2.75573192239859e-7
                       * (-2.9540329644556910 + lnx)
                       * (-2.5491366926991850 + lnx)
                       * (-2.1348279229279880 + lnx)
                       * (-1.6741881076349450 + lnx)
                       * (-1.1325949616098420 + lnx)
                       * (-0.4590034650618494 + lnx)
                       * ( 0.4399352987435699 + lnx)
                       * ( 1.7702236517651670 + lnx)
                       * ( 4.1231539047474080 + lnx)
                       * ( 10.342627908148680 + lnx);

    term1 = a*(c1+a*(c2+a*(c3+a*(c4+a*(c5+a*(c6+a*(c7+a*(c8+a*(c9+a*c10)))))))));
  }

  {
    /* Evaluate the sum.
     */
    const int nmax = 5000;
    double t = 1.0;
    int n;
    sum = 1.0;

    for(n=1; n<nmax; n++) {
      t *= -x/(n+1.0);
      sum += (a+1.0)/(a+n+1.0)*t;
      if(fabs(t/sum) < GSL_DBL_EPSILON) break;
    }

    if(n == nmax)
      stat_sum = GSL_EMAXITER;
    else
      stat_sum = GSL_SUCCESS;
  }

  term2 = (1.0 - term1) * a/(a+1.0) * x * sum;
  result->val  = term1 + term2;
  result->err  = GSL_DBL_EPSILON * (fabs(term1) + 2.0*fabs(term2));
  result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  return stat_sum;
}


/* series for small a and x, but not defined for a == 0 */
static int
gamma_inc_series(double a, double x, gsl_sf_result * result)
{
  gsl_sf_result Q;
  gsl_sf_result G;
  const int stat_Q = gamma_inc_Q_series(a, x, &Q);
  const int stat_G = gsl_sf_gamma_e(a, &G);
  result->val = Q.val * G.val;
  result->err = fabs(Q.val * G.err) + fabs(Q.err * G.val);
  result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);

  return GSL_ERROR_SELECT_2(stat_Q, stat_G);
}


static int
gamma_inc_a_gt_0(double a, double x, gsl_sf_result * result)
{
  /* x > 0 and a > 0; use result for Q */
  gsl_sf_result Q;
  gsl_sf_result G;
  const int stat_Q = gsl_sf_gamma_inc_Q_e(a, x, &Q);
  const int stat_G = gsl_sf_gamma_e(a, &G);

  result->val = G.val * Q.val;
  result->err = fabs(G.val * Q.err) + fabs(G.err * Q.val);
  result->err += 2.0*GSL_DBL_EPSILON * fabs(result->val);

  return GSL_ERROR_SELECT_2(stat_G, stat_Q);
}


static int
gamma_inc_CF(double a, double x, gsl_sf_result * result)
{
  gsl_sf_result F;
  gsl_sf_result pre;
  const int stat_F = gamma_inc_F_CF(a, x, &F);
  const int stat_E = gsl_sf_exp_e((a-1.0)*log(x) - x, &pre);

  result->val = F.val * pre.val;
  result->err = fabs(F.err * pre.val) + fabs(F.val * pre.err);
  result->err += (2.0 + fabs(a)) * GSL_DBL_EPSILON * fabs(result->val);

  return GSL_ERROR_SELECT_2(stat_F, stat_E);
}


/* evaluate Gamma(0,x), x > 0 */
#define GAMMA_INC_A_0(x, result) gsl_sf_expint_E1_e(x, result)


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_gamma_inc_Q_e(const double a, const double x, gsl_sf_result * result)
{
  if(a < 0.0 || x < 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x == 0.0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(a == 0.0)
  {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(x <= 0.5*a) {
    /* If the series is quick, do that. It is
     * robust and simple.
     */
    gsl_sf_result P;
    int stat_P = gamma_inc_P_series(a, x, &P);
    result->val  = 1.0 - P.val;
    result->err  = P.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_P;
  }
  else if(a >= 1.0e+06 && (x-a)*(x-a) < a) {
    /* Then try the difficult asymptotic regime.
     * This is the only way to do this region.
     */
    return gamma_inc_Q_asymp_unif(a, x, result);
  }
  else if(a < 0.2 && x < 5.0) {
    /* Cancellations at small a must be handled
     * analytically; x should not be too big
     * either since the series terms grow
     * with x and log(x).
     */
    return gamma_inc_Q_series(a, x, result);
  }
  else if(a <= x) {
    if(x <= 1.0e+06) {
      /* Continued fraction is excellent for x >~ a.
       * We do not let x be too large when x > a since
       * it is somewhat pointless to try this there;
       * the function is rapidly decreasing for
       * x large and x > a, and it will just
       * underflow in that region anyway. We
       * catch that case in the standard
       * large-x method.
       */
      return gamma_inc_Q_CF(a, x, result);
    }
    else {
      return gamma_inc_Q_large_x(a, x, result);
    }
  }
  else {
    if(x > a - sqrt(a)) {
      /* Continued fraction again. The convergence
       * is a little slower here, but that is fine.
       * We have to trade that off against the slow
       * convergence of the series, which is the
       * only other option.
       */
      return gamma_inc_Q_CF(a, x, result);
    }
    else {
      gsl_sf_result P;
      int stat_P = gamma_inc_P_series(a, x, &P);
      result->val  = 1.0 - P.val;
      result->err  = P.err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return stat_P;
    }
  }
}


int
gsl_sf_gamma_inc_P_e(const double a, const double x, gsl_sf_result * result)
{
  if(a <= 0.0 || x < 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(x < 20.0 || x < 0.5*a) {
    /* Do the easy series cases. Robust and quick.
     */
    return gamma_inc_P_series(a, x, result);
  }
  else if(a > 1.0e+06 && (x-a)*(x-a) < a) {
    /* Crossover region. Note that Q and P are
     * roughly the same order of magnitude here,
     * so the subtraction is stable.
     */
    gsl_sf_result Q;
    int stat_Q = gamma_inc_Q_asymp_unif(a, x, &Q);
    result->val  = 1.0 - Q.val;
    result->err  = Q.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_Q;
  }
  else if(a <= x) {
    /* Q <~ P in this area, so the
     * subtractions are stable.
     */
    gsl_sf_result Q;
    int stat_Q;
    if(a > 0.2*x) {
      stat_Q = gamma_inc_Q_CF(a, x, &Q);
    }
    else {
      stat_Q = gamma_inc_Q_large_x(a, x, &Q);
    }
    result->val  = 1.0 - Q.val;
    result->err  = Q.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_Q;
  }
  else {
    if((x-a)*(x-a) < a) {
      /* This condition is meant to insure
       * that Q is not very close to 1,
       * so the subtraction is stable.
       */
      gsl_sf_result Q;
      int stat_Q = gamma_inc_Q_CF(a, x, &Q);
      result->val  = 1.0 - Q.val;
      result->err  = Q.err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return stat_Q;
    }
    else {
      return gamma_inc_P_series(a, x, result);
    }
  }
}


int
gsl_sf_gamma_inc_e(const double a, const double x, gsl_sf_result * result)
{
  if(x < 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x == 0.0) {
    return gsl_sf_gamma_e(a, result);
  }
  else if(a == 0.0)
  {
    return GAMMA_INC_A_0(x, result);
  }
  else if(a > 0.0)
  {
    return gamma_inc_a_gt_0(a, x, result);
  }
  else if(x > 0.25)
  {
    /* continued fraction seems to fail for x too small; otherwise
       it is ok, independent of the value of |x/a|, because of the
       non-oscillation in the expansion, i.e. the CF is
       un-conditionally convergent for a < 0 and x > 0
     */
    return gamma_inc_CF(a, x, result);
  }
  else if(fabs(a) < 0.5)
  {
    return gamma_inc_series(a, x, result);
  }
  else
  {
    /* a = fa + da; da >= 0 */
    const double fa = floor(a);
    const double da = a - fa;

    gsl_sf_result g_da;
    const int stat_g_da = ( da > 0.0 ? gamma_inc_a_gt_0(da, x, &g_da)
                                     : GAMMA_INC_A_0(x, &g_da));

    double alpha = da;
    double gax = g_da.val;

    /* Gamma(alpha-1,x) = 1/(alpha-1) (Gamma(a,x) - x^(alpha-1) e^-x) */
    do
    {
      const double shift = exp(-x + (alpha-1.0)*log(x));
      gax = (gax - shift) / (alpha - 1.0);
      alpha -= 1.0;
    } while(alpha > a);

    result->val = gax;
    result->err = 2.0*(1.0 + fabs(a))*GSL_DBL_EPSILON*fabs(gax);
    return stat_g_da;
  }

}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_gamma_inc_P(const double a, const double x)
{
  EVAL_RESULT(gsl_sf_gamma_inc_P_e(a, x, &result));
}

double gsl_sf_gamma_inc_Q(const double a, const double x)
{
  EVAL_RESULT(gsl_sf_gamma_inc_Q_e(a, x, &result));
}

double gsl_sf_gamma_inc(const double a, const double x)
{
   EVAL_RESULT(gsl_sf_gamma_inc_e(a, x, &result));
}
} // gamma_inc.c

// infnan.c ///////////////////////////////////////////////////////////////////////
extern "C" { 
double gsl_nan (void);
double gsl_posinf (void);
double gsl_neginf (void);
double gsl_fdiv (const double x, const double y);

double gsl_nan (void)
{
  return gsl_fdiv (0.0, 0.0);
}

double gsl_posinf (void)
{
  return gsl_fdiv (+1.0, 0.0);
}

double gsl_neginf (void)
{
  return gsl_fdiv (-1.0, 0.0);
}


int gsl_isnan (const double x);
int gsl_isinf (const double x);
int gsl_finite (const double x);

int gsl_isnan (const double x)
{
  int status = (x != x);
  return status;
}

int gsl_isinf (const double x)
{
  double y = x - x;
  int s = (y != y);

  if (s && x > 0)
    return +1;
  else if (s && x < 0)
    return -1;
  else
    return 0;
}

int gsl_finite (const double x)
{
  const double y = x - x;
  int status = (y == y);
  return status;
}

} // infnan.c

// log.c ///////////////////////////////////////////////////////////////////////
extern "C" { 
/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* Chebyshev expansion for log(1 + x(t))/x(t)
 *
 * x(t) = (4t-1)/(2(4-t))
 * t(x) = (8x+1)/(2(x+2))
 * -1/2 < x < 1/2
 * -1 < t < 1
 */
static double lopx_data[21] = {
  2.16647910664395270521272590407,
 -0.28565398551049742084877469679,
  0.01517767255690553732382488171,
 -0.00200215904941415466274422081,
  0.00019211375164056698287947962,
 -0.00002553258886105542567601400,
  2.9004512660400621301999384544e-06,
 -3.8873813517057343800270917900e-07,
  4.7743678729400456026672697926e-08,
 -6.4501969776090319441714445454e-09,
  8.2751976628812389601561347296e-10,
 -1.1260499376492049411710290413e-10,
  1.4844576692270934446023686322e-11,
 -2.0328515972462118942821556033e-12,
  2.7291231220549214896095654769e-13,
 -3.7581977830387938294437434651e-14,
  5.1107345870861673561462339876e-15,
 -7.0722150011433276578323272272e-16,
  9.7089758328248469219003866867e-17,
 -1.3492637457521938883731579510e-17,
  1.8657327910677296608121390705e-18
};
static cheb_series lopx_cs = {
  lopx_data,
  20,
  -1, 1,
  10
};

/* Chebyshev expansion for (log(1 + x(t)) - x(t))/x(t)^2
 *
 * x(t) = (4t-1)/(2(4-t))
 * t(x) = (8x+1)/(2(x+2))
 * -1/2 < x < 1/2
 * -1 < t < 1
 */
static double lopxmx_data[20] = {
 -1.12100231323744103373737274541,
  0.19553462773379386241549597019,
 -0.01467470453808083971825344956,
  0.00166678250474365477643629067,
 -0.00018543356147700369785746902,
  0.00002280154021771635036301071,
 -2.8031253116633521699214134172e-06,
  3.5936568872522162983669541401e-07,
 -4.6241857041062060284381167925e-08,
  6.0822637459403991012451054971e-09,
 -8.0339824424815790302621320732e-10,
  1.0751718277499375044851551587e-10,
 -1.4445310914224613448759230882e-11,
  1.9573912180610336168921438426e-12,
 -2.6614436796793061741564104510e-13,
  3.6402634315269586532158344584e-14,
 -4.9937495922755006545809120531e-15,
  6.8802890218846809524646902703e-16,
 -9.5034129794804273611403251480e-17,
  1.3170135013050997157326965813e-17
};
static cheb_series lopxmx_cs = {
  lopxmx_data,
  19,
  -1, 1,
  9
};


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_log_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else {
    result->val = log(x);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int
gsl_sf_log_abs_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x == 0.0) {
    DOMAIN_ERROR(result);
  }
  else {
    result->val = log(fabs(x));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}

int
gsl_sf_complex_log_e(const double zr, const double zi, gsl_sf_result * lnr, gsl_sf_result * theta)
{
  /* CHECK_POINTER(lnr) */
  /* CHECK_POINTER(theta) */

  if(zr != 0.0 || zi != 0.0) {
    const double ax = fabs(zr);
    const double ay = fabs(zi);
    const double min = GSL_MIN(ax, ay);
    const double max = GSL_MAX(ax, ay);
    lnr->val = log(max) + 0.5 * log(1.0 + (min/max)*(min/max));
    lnr->err = 2.0 * GSL_DBL_EPSILON * fabs(lnr->val);
    theta->val = atan2(zi, zr);
    theta->err = GSL_DBL_EPSILON * fabs(lnr->val);
    return GSL_SUCCESS;
  }
  else {
    DOMAIN_ERROR_2(lnr, theta);
  }
}


int
gsl_sf_log_1plusx_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= -1.0) {
    DOMAIN_ERROR(result);
  }
  else if(fabs(x) < GSL_ROOT6_DBL_EPSILON) {
    const double c1 = -0.5;
    const double c2 =  1.0/3.0;
    const double c3 = -1.0/4.0;
    const double c4 =  1.0/5.0;
    const double c5 = -1.0/6.0;
    const double c6 =  1.0/7.0;
    const double c7 = -1.0/8.0;
    const double c8 =  1.0/9.0;
    const double c9 = -1.0/10.0;
    const double t  =  c5 + x*(c6 + x*(c7 + x*(c8 + x*c9)));
    result->val = x * (1.0 + x*(c1 + x*(c2 + x*(c3 + x*(c4 + x*t)))));
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(fabs(x) < 0.5) {
    double t = 0.5*(8.0*x + 1.0)/(x+2.0);
    gsl_sf_result c;
    cheb_eval_e(&lopx_cs, t, &c);
    result->val = x * c.val;
    result->err = fabs(x * c.err);
    return GSL_SUCCESS;
  }
  else {
    result->val = log(1.0 + x);
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int
gsl_sf_log_1plusx_mx_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= -1.0) {
    DOMAIN_ERROR(result);
  }
  else if(fabs(x) < GSL_ROOT5_DBL_EPSILON) {
    const double c1 = -0.5;
    const double c2 =  1.0/3.0;
    const double c3 = -1.0/4.0;
    const double c4 =  1.0/5.0;
    const double c5 = -1.0/6.0;
    const double c6 =  1.0/7.0;
    const double c7 = -1.0/8.0;
    const double c8 =  1.0/9.0;
    const double c9 = -1.0/10.0;
    const double t  =  c5 + x*(c6 + x*(c7 + x*(c8 + x*c9)));
    result->val = x*x * (c1 + x*(c2 + x*(c3 + x*(c4 + x*t))));
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(fabs(x) < 0.5) {
    double t = 0.5*(8.0*x + 1.0)/(x+2.0);
    gsl_sf_result c;
    cheb_eval_e(&lopxmx_cs, t, &c);
    result->val = x*x * c.val;
    result->err = x*x * c.err;
    return GSL_SUCCESS;
  }
  else {
    const double lterm = log(1.0 + x);
    result->val = lterm - x;
    result->err = GSL_DBL_EPSILON * (fabs(lterm) + fabs(x));
    return GSL_SUCCESS;
  }
}



/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_log(const double x)
{
  EVAL_RESULT(gsl_sf_log_e(x, &result));
}

double gsl_sf_log_abs(const double x)
{
  EVAL_RESULT(gsl_sf_log_abs_e(x, &result));
}

double gsl_sf_log_1plusx(const double x)
{
  EVAL_RESULT(gsl_sf_log_1plusx_e(x, &result));
}

double gsl_sf_log_1plusx_mx(const double x)
{
  EVAL_RESULT(gsl_sf_log_1plusx_mx_e(x, &result));
}
} // log.c

// math.c //////////////////////////////////////////////////////////////////////
extern "C" {
/**********************************************************************
 * Complex numbers
 **********************************************************************/

#ifndef HIDE_INLINE_STATIC
gsl_complex
gsl_complex_rect (double x, double y)
{                               /* return z = x + i y */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, x, y);
  return z;
}
#endif

gsl_complex
gsl_complex_polar (double r, double theta)
{                               /* return z = r exp(i theta) */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, r * cos (theta), r * sin (theta));
  return z;
}

/**********************************************************************
 * Properties of complex numbers
 **********************************************************************/

double
gsl_complex_arg (gsl_complex z)
{                               /* return arg(z),  -pi < arg(z) <= +pi */
  double x = GSL_REAL (z);
  double y = GSL_IMAG (z);

  if (x == 0.0 && y == 0.0)
    {
      return 0;
    }

  return atan2 (y, x);
}

double
gsl_complex_abs (gsl_complex z)
{                               /* return |z| */
  return hypot (GSL_REAL (z), GSL_IMAG (z));
}

double
gsl_complex_abs2 (gsl_complex z)
{                               /* return |z|^2 */
  double x = GSL_REAL (z);
  double y = GSL_IMAG (z);

  return (x * x + y * y);
}

double
gsl_complex_logabs (gsl_complex z)
{                               /* return log|z| */
  double xabs = fabs (GSL_REAL (z));
  double yabs = fabs (GSL_IMAG (z));
  double max, u;

  if (xabs >= yabs)
    {
      max = xabs;
      u = yabs / xabs;
    }
  else
    {
      max = yabs;
      u = xabs / yabs;
    }

  /* Handle underflow when u is close to 0 */

  return log (max) + 0.5 * log1p (u * u);
}


/***********************************************************************
 * Complex arithmetic operators
 ***********************************************************************/

gsl_complex
gsl_complex_add (gsl_complex a, gsl_complex b)
{                               /* z=a+b */
  double ar = GSL_REAL (a), ai = GSL_IMAG (a);
  double br = GSL_REAL (b), bi = GSL_IMAG (b);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, ar + br, ai + bi);
  return z;
}

gsl_complex
gsl_complex_add_real (gsl_complex a, double x)
{                               /* z=a+x */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, GSL_REAL (a) + x, GSL_IMAG (a));
  return z;
}

gsl_complex
gsl_complex_add_imag (gsl_complex a, double y)
{                               /* z=a+iy */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, GSL_REAL (a), GSL_IMAG (a) + y);
  return z;
}


gsl_complex
gsl_complex_sub (gsl_complex a, gsl_complex b)
{                               /* z=a-b */
  double ar = GSL_REAL (a), ai = GSL_IMAG (a);
  double br = GSL_REAL (b), bi = GSL_IMAG (b);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, ar - br, ai - bi);
  return z;
}

gsl_complex
gsl_complex_sub_real (gsl_complex a, double x)
{                               /* z=a-x */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, GSL_REAL (a) - x, GSL_IMAG (a));
  return z;
}

gsl_complex
gsl_complex_sub_imag (gsl_complex a, double y)
{                               /* z=a-iy */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, GSL_REAL (a), GSL_IMAG (a) - y);
  return z;
}

gsl_complex
gsl_complex_mul (gsl_complex a, gsl_complex b)
{                               /* z=a*b */
  double ar = GSL_REAL (a), ai = GSL_IMAG (a);
  double br = GSL_REAL (b), bi = GSL_IMAG (b);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, ar * br - ai * bi, ar * bi + ai * br);
  return z;
}

gsl_complex
gsl_complex_mul_real (gsl_complex a, double x)
{                               /* z=a*x */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, x * GSL_REAL (a), x * GSL_IMAG (a));
  return z;
}

gsl_complex
gsl_complex_mul_imag (gsl_complex a, double y)
{                               /* z=a*iy */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, -y * GSL_IMAG (a), y * GSL_REAL (a));
  return z;
}

gsl_complex
gsl_complex_div (gsl_complex a, gsl_complex b)
{                               /* z=a/b */
  double ar = GSL_REAL (a), ai = GSL_IMAG (a);
  double br = GSL_REAL (b), bi = GSL_IMAG (b);

  double s = 1.0 / gsl_complex_abs (b);

  double sbr = s * br;
  double sbi = s * bi;

  double zr = (ar * sbr + ai * sbi) * s;
  double zi = (ai * sbr - ar * sbi) * s;

  gsl_complex z;
  GSL_SET_COMPLEX (&z, zr, zi);
  return z;
}

gsl_complex
gsl_complex_div_real (gsl_complex a, double x)
{                               /* z=a/x */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, GSL_REAL (a) / x, GSL_IMAG (a) / x);
  return z;
}

gsl_complex
gsl_complex_div_imag (gsl_complex a, double y)
{                               /* z=a/(iy) */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, GSL_IMAG (a) / y,  - GSL_REAL (a) / y);
  return z;
}

gsl_complex
gsl_complex_conjugate (gsl_complex a)
{                               /* z=conj(a) */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, GSL_REAL (a), -GSL_IMAG (a));
  return z;
}

gsl_complex
gsl_complex_negative (gsl_complex a)
{                               /* z=-a */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, -GSL_REAL (a), -GSL_IMAG (a));
  return z;
}

gsl_complex
gsl_complex_inverse (gsl_complex a)
{                               /* z=1/a */
  double s = 1.0 / gsl_complex_abs (a);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, (GSL_REAL (a) * s) * s, -(GSL_IMAG (a) * s) * s);
  return z;
}

/**********************************************************************
 * Elementary complex functions
 **********************************************************************/

gsl_complex
gsl_complex_sqrt (gsl_complex a)
{                               /* z=sqrt(a) */
  gsl_complex z;

  if (GSL_REAL (a) == 0.0 && GSL_IMAG (a) == 0.0)
    {
      GSL_SET_COMPLEX (&z, 0, 0);
    }
  else
    {
      double x = fabs (GSL_REAL (a));
      double y = fabs (GSL_IMAG (a));
      double w;

      if (x >= y)
        {
          double t = y / x;
          w = sqrt (x) * sqrt (0.5 * (1.0 + sqrt (1.0 + t * t)));
        }
      else
        {
          double t = x / y;
          w = sqrt (y) * sqrt (0.5 * (t + sqrt (1.0 + t * t)));
        }

      if (GSL_REAL (a) >= 0.0)
        {
          double ai = GSL_IMAG (a);
          GSL_SET_COMPLEX (&z, w, ai / (2.0 * w));
        }
      else
        {
          double ai = GSL_IMAG (a);
          double vi = (ai >= 0) ? w : -w;
          GSL_SET_COMPLEX (&z, ai / (2.0 * vi), vi);
        }
    }

  return z;
}

gsl_complex
gsl_complex_sqrt_real (double x)
{                               /* z=sqrt(x) */
  gsl_complex z;

  if (x >= 0)
    {
      GSL_SET_COMPLEX (&z, sqrt (x), 0.0);
    }
  else
    {
      GSL_SET_COMPLEX (&z, 0.0, sqrt (-x));
    }

  return z;
}

gsl_complex
gsl_complex_exp (gsl_complex a)
{                               /* z=exp(a) */
  double rho = exp (GSL_REAL (a));
  double theta = GSL_IMAG (a);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, rho * cos (theta), rho * sin (theta));
  return z;
}

gsl_complex
gsl_complex_pow (gsl_complex a, gsl_complex b)
{                               /* z=a^b */
  gsl_complex z;

  if (GSL_REAL (a) == 0 && GSL_IMAG (a) == 0.0)
    {
      GSL_SET_COMPLEX (&z, 0.0, 0.0);
    }
  else
    {
      double logr = gsl_complex_logabs (a);
      double theta = gsl_complex_arg (a);

      double br = GSL_REAL (b), bi = GSL_IMAG (b);

      double rho = exp (logr * br - bi * theta);
      double beta = theta * br + bi * logr;

      GSL_SET_COMPLEX (&z, rho * cos (beta), rho * sin (beta));
    }

  return z;
}

gsl_complex
gsl_complex_pow_real (gsl_complex a, double b)
{                               /* z=a^b */
  gsl_complex z;

  if (GSL_REAL (a) == 0 && GSL_IMAG (a) == 0)
    {
      GSL_SET_COMPLEX (&z, 0, 0);
    }
  else
    {
      double logr = gsl_complex_logabs (a);
      double theta = gsl_complex_arg (a);
      double rho = exp (logr * b);
      double beta = theta * b;
      GSL_SET_COMPLEX (&z, rho * cos (beta), rho * sin (beta));
    }

  return z;
}

gsl_complex
gsl_complex_log (gsl_complex a)
{                               /* z=log(a) */
  double logr = gsl_complex_logabs (a);
  double theta = gsl_complex_arg (a);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, logr, theta);
  return z;
}

gsl_complex
gsl_complex_log10 (gsl_complex a)
{                               /* z = log10(a) */
  return gsl_complex_mul_real (gsl_complex_log (a), 1 / log (10.));
}

gsl_complex
gsl_complex_log_b (gsl_complex a, gsl_complex b)
{
  return gsl_complex_div (gsl_complex_log (a), gsl_complex_log (b));
}

/***********************************************************************
 * Complex trigonometric functions
 ***********************************************************************/

gsl_complex
gsl_complex_sin (gsl_complex a)
{                               /* z = sin(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);

  gsl_complex z;

  if (I == 0.0) 
    {
      /* avoid returing negative zero (-0.0) for the imaginary part  */

      GSL_SET_COMPLEX (&z, sin (R), 0.0);  
    } 
  else 
    {
      GSL_SET_COMPLEX (&z, sin (R) * cosh (I), cos (R) * sinh (I));
    }

  return z;
}

gsl_complex
gsl_complex_cos (gsl_complex a)
{                               /* z = cos(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);

  gsl_complex z;

  if (I == 0.0) 
    {
      /* avoid returing negative zero (-0.0) for the imaginary part  */

      GSL_SET_COMPLEX (&z, cos (R), 0.0);  
    } 
  else 
    {
      GSL_SET_COMPLEX (&z, cos (R) * cosh (I), sin (R) * sinh (-I));
    }

  return z;
}

gsl_complex
gsl_complex_tan (gsl_complex a)
{                               /* z = tan(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);

  gsl_complex z;

  if (fabs (I) < 1)
    {
      double D = pow (cos (R), 2.0) + pow (sinh (I), 2.0);

      GSL_SET_COMPLEX (&z, 0.5 * sin (2 * R) / D, 0.5 * sinh (2 * I) / D);
    }
  else
    {
      double u = exp (-I);
      double C = 2 * u / (1 - pow (u, 2.0));
      double D = 1 + pow (cos (R), 2.0) * pow (C, 2.0);

      double S = pow (C, 2.0);
      double T = 1.0 / tanh (I);

      GSL_SET_COMPLEX (&z, 0.5 * sin (2 * R) * S / D, T / D);
    }

  return z;
}

gsl_complex
gsl_complex_sec (gsl_complex a)
{                               /* z = sec(a) */
  gsl_complex z = gsl_complex_cos (a);
  return gsl_complex_inverse (z);
}

gsl_complex
gsl_complex_csc (gsl_complex a)
{                               /* z = csc(a) */
  gsl_complex z = gsl_complex_sin (a);
  return gsl_complex_inverse(z);
}


gsl_complex
gsl_complex_cot (gsl_complex a)
{                               /* z = cot(a) */
  gsl_complex z = gsl_complex_tan (a);
  return gsl_complex_inverse (z);
}

/**********************************************************************
 * Inverse Complex Trigonometric Functions
 **********************************************************************/

gsl_complex
gsl_complex_arcsin (gsl_complex a)
{                               /* z = arcsin(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);
  gsl_complex z;

  if (I == 0)
    {
      z = gsl_complex_arcsin_real (R);
    }
  else
    {
      double x = fabs (R), y = fabs (I);
      double r = hypot (x + 1, y), s = hypot (x - 1, y);
      double A = 0.5 * (r + s);
      double B = x / A;
      double y2 = y * y;

      double real, imag;

      const double A_crossover = 1.5, B_crossover = 0.6417;

      if (B <= B_crossover)
        {
          real = asin (B);
        }
      else
        {
          if (x <= 1)
            {
              double D = 0.5 * (A + x) * (y2 / (r + x + 1) + (s + (1 - x)));
              real = atan (x / sqrt (D));
            }
          else
            {
              double Apx = A + x;
              double D = 0.5 * (Apx / (r + x + 1) + Apx / (s + (x - 1)));
              real = atan (x / (y * sqrt (D)));
            }
        }

      if (A <= A_crossover)
        {
          double Am1;

          if (x < 1)
            {
              Am1 = 0.5 * (y2 / (r + (x + 1)) + y2 / (s + (1 - x)));
            }
          else
            {
              Am1 = 0.5 * (y2 / (r + (x + 1)) + (s + (x - 1)));
            }

          imag = log1p (Am1 + sqrt (Am1 * (A + 1)));
        }
      else
        {
          imag = log (A + sqrt (A * A - 1));
        }

      GSL_SET_COMPLEX (&z, (R >= 0) ? real : -real, (I >= 0) ? imag : -imag);
    }

  return z;
}

gsl_complex
gsl_complex_arcsin_real (double a)
{                               /* z = arcsin(a) */
  gsl_complex z;

  if (fabs (a) <= 1.0)
    {
      GSL_SET_COMPLEX (&z, asin (a), 0.0);
    }
  else
    {
      if (a < 0.0)
        {
          GSL_SET_COMPLEX (&z, -M_PI_2, acosh (-a));
        }
      else
        {
          GSL_SET_COMPLEX (&z, M_PI_2, -acosh (a));
        }
    }

  return z;
}

gsl_complex
gsl_complex_arccos (gsl_complex a)
{                               /* z = arccos(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);
  gsl_complex z;

  if (I == 0)
    {
      z = gsl_complex_arccos_real (R);
    }
  else
    {
      double x = fabs (R), y = fabs (I);
      double r = hypot (x + 1, y), s = hypot (x - 1, y);
      double A = 0.5 * (r + s);
      double B = x / A;
      double y2 = y * y;

      double real, imag;

      const double A_crossover = 1.5, B_crossover = 0.6417;

      if (B <= B_crossover)
        {
          real = acos (B);
        }
      else
        {
          if (x <= 1)
            {
              double D = 0.5 * (A + x) * (y2 / (r + x + 1) + (s + (1 - x)));
              real = atan (sqrt (D) / x);
            }
          else
            {
              double Apx = A + x;
              double D = 0.5 * (Apx / (r + x + 1) + Apx / (s + (x - 1)));
              real = atan ((y * sqrt (D)) / x);
            }
        }

      if (A <= A_crossover)
        {
          double Am1;

          if (x < 1)
            {
              Am1 = 0.5 * (y2 / (r + (x + 1)) + y2 / (s + (1 - x)));
            }
          else
            {
              Am1 = 0.5 * (y2 / (r + (x + 1)) + (s + (x - 1)));
            }

          imag = log1p (Am1 + sqrt (Am1 * (A + 1)));
        }
      else
        {
          imag = log (A + sqrt (A * A - 1));
        }

      GSL_SET_COMPLEX (&z, (R >= 0) ? real : M_PI - real, (I >= 0) ? -imag : imag);
    }

  return z;
}

gsl_complex
gsl_complex_arccos_real (double a)
{                               /* z = arccos(a) */
  gsl_complex z;

  if (fabs (a) <= 1.0)
    {
      GSL_SET_COMPLEX (&z, acos (a), 0);
    }
  else
    {
      if (a < 0.0)
        {
          GSL_SET_COMPLEX (&z, M_PI, -acosh (-a));
        }
      else
        {
          GSL_SET_COMPLEX (&z, 0, acosh (a));
        }
    }

  return z;
}

gsl_complex
gsl_complex_arctan (gsl_complex a)
{                               /* z = arctan(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);
  gsl_complex z;

  if (I == 0)
    {
      GSL_SET_COMPLEX (&z, atan (R), 0);
    }
  else
    {
      /* FIXME: This is a naive implementation which does not fully
         take into account cancellation errors, overflow, underflow
         etc.  It would benefit from the Hull et al treatment. */

      double r = hypot (R, I);

      double imag;

      double u = 2 * I / (1 + r * r);

      /* FIXME: the following cross-over should be optimized but 0.1
         seems to work ok */

      if (fabs (u) < 0.1)
        {
          imag = 0.25 * (log1p (u) - log1p (-u));
        }
      else
        {
          double A = hypot (R, I + 1);
          double B = hypot (R, I - 1);
          imag = 0.5 * log (A / B);
        }

      if (R == 0)
        {
          if (I > 1)
            {
              GSL_SET_COMPLEX (&z, M_PI_2, imag);
            }
          else if (I < -1)
            {
              GSL_SET_COMPLEX (&z, -M_PI_2, imag);
            }
          else
            {
              GSL_SET_COMPLEX (&z, 0, imag);
            };
        }
      else
        {
          GSL_SET_COMPLEX (&z, 0.5 * atan2 (2 * R, ((1 + r) * (1 - r))), imag);
        }
    }

  return z;
}

gsl_complex
gsl_complex_arcsec (gsl_complex a)
{                               /* z = arcsec(a) */
  gsl_complex z = gsl_complex_inverse (a);
  return gsl_complex_arccos (z);
}

gsl_complex
gsl_complex_arcsec_real (double a)
{                               /* z = arcsec(a) */
  gsl_complex z;

  if (a <= -1.0 || a >= 1.0)
    {
      GSL_SET_COMPLEX (&z, acos (1 / a), 0.0);
    }
  else
    {
      if (a >= 0.0)
        {
          GSL_SET_COMPLEX (&z, 0, acosh (1 / a));
        }
      else
        {
          GSL_SET_COMPLEX (&z, M_PI, -acosh (-1 / a));
        }
    }

  return z;
}

gsl_complex
gsl_complex_arccsc (gsl_complex a)
{                               /* z = arccsc(a) */
  gsl_complex z = gsl_complex_inverse (a);
  return gsl_complex_arcsin (z);
}

gsl_complex
gsl_complex_arccsc_real (double a)
{                               /* z = arccsc(a) */
  gsl_complex z;

  if (a <= -1.0 || a >= 1.0)
    {
      GSL_SET_COMPLEX (&z, asin (1 / a), 0.0);
    }
  else
    {
      if (a >= 0.0)
        {
          GSL_SET_COMPLEX (&z, M_PI_2, -acosh (1 / a));
        }
      else
        {
          GSL_SET_COMPLEX (&z, -M_PI_2, acosh (-1 / a));
        }
    }

  return z;
}

gsl_complex
gsl_complex_arccot (gsl_complex a)
{                               /* z = arccot(a) */
  gsl_complex z;

  if (GSL_REAL (a) == 0.0 && GSL_IMAG (a) == 0.0)
    {
      GSL_SET_COMPLEX (&z, M_PI_2, 0);
    }
  else
    {
      z = gsl_complex_inverse (a);
      z = gsl_complex_arctan (z);
    }

  return z;
}

/**********************************************************************
 * Complex Hyperbolic Functions
 **********************************************************************/

gsl_complex
gsl_complex_sinh (gsl_complex a)
{                               /* z = sinh(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, sinh (R) * cos (I), cosh (R) * sin (I));
  return z;
}

gsl_complex
gsl_complex_cosh (gsl_complex a)
{                               /* z = cosh(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, cosh (R) * cos (I), sinh (R) * sin (I));
  return z;
}

gsl_complex
gsl_complex_tanh (gsl_complex a)
{                               /* z = tanh(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);

  gsl_complex z;

  if (fabs(R) < 1.0) 
    {
      double D = pow (cos (I), 2.0) + pow (sinh (R), 2.0);
      
      GSL_SET_COMPLEX (&z, sinh (R) * cosh (R) / D, 0.5 * sin (2 * I) / D);
    }
  else
    {
      double D = pow (cos (I), 2.0) + pow (sinh (R), 2.0);
      double F = 1 + pow (cos (I) / sinh (R), 2.0);

      GSL_SET_COMPLEX (&z, 1.0 / (tanh (R) * F), 0.5 * sin (2 * I) / D);
    }

  return z;
}

gsl_complex
gsl_complex_sech (gsl_complex a)
{                               /* z = sech(a) */
  gsl_complex z = gsl_complex_cosh (a);
  return gsl_complex_inverse (z);
}

gsl_complex
gsl_complex_csch (gsl_complex a)
{                               /* z = csch(a) */
  gsl_complex z = gsl_complex_sinh (a);
  return gsl_complex_inverse (z);
}

gsl_complex
gsl_complex_coth (gsl_complex a)
{                               /* z = coth(a) */
  gsl_complex z = gsl_complex_tanh (a);
  return gsl_complex_inverse (z);
}

/**********************************************************************
 * Inverse Complex Hyperbolic Functions
 **********************************************************************/

gsl_complex
gsl_complex_arcsinh (gsl_complex a)
{                               /* z = arcsinh(a) */
  gsl_complex z = gsl_complex_mul_imag(a, 1.0);
  z = gsl_complex_arcsin (z);
  z = gsl_complex_mul_imag (z, -1.0);
  return z;
}

gsl_complex
gsl_complex_arccosh (gsl_complex a)
{                               /* z = arccosh(a) */
  gsl_complex z = gsl_complex_arccos (a);
  z = gsl_complex_mul_imag (z, GSL_IMAG(z) > 0 ? -1.0 : 1.0);
  return z;
}

gsl_complex
gsl_complex_arccosh_real (double a)
{                               /* z = arccosh(a) */
  gsl_complex z;

  if (a >= 1)
    {
      GSL_SET_COMPLEX (&z, acosh (a), 0);
    }
  else
    {
      if (a >= -1.0)
        {
          GSL_SET_COMPLEX (&z, 0, acos (a));
        }
      else
        {
          GSL_SET_COMPLEX (&z, acosh (-a), M_PI);
        }
    }

  return z;
}

gsl_complex
gsl_complex_arctanh (gsl_complex a)
{                               /* z = arctanh(a) */
  if (GSL_IMAG (a) == 0.0)
    {
      return gsl_complex_arctanh_real (GSL_REAL (a));
    }
  else
    {
      gsl_complex z = gsl_complex_mul_imag(a, 1.0);
      z = gsl_complex_arctan (z);
      z = gsl_complex_mul_imag (z, -1.0);
      return z;
    }
}

gsl_complex
gsl_complex_arctanh_real (double a)
{                               /* z = arctanh(a) */
  gsl_complex z;

  if (a > -1.0 && a < 1.0)
    {
      GSL_SET_COMPLEX (&z, atanh (a), 0);
    }
  else
    {
      GSL_SET_COMPLEX (&z, atanh (1 / a), (a < 0) ? M_PI_2 : -M_PI_2);
    }

  return z;
}

gsl_complex
gsl_complex_arcsech (gsl_complex a)
{                               /* z = arcsech(a); */
  gsl_complex t = gsl_complex_inverse (a);
  return gsl_complex_arccosh (t);
}

gsl_complex
gsl_complex_arccsch (gsl_complex a)
{                               /* z = arccsch(a) */
  gsl_complex t = gsl_complex_inverse (a);
  return gsl_complex_arcsinh (t);
}

gsl_complex
gsl_complex_arccoth (gsl_complex a)
{                               /* z = arccoth(a) */
  gsl_complex t = gsl_complex_inverse (a);
  return gsl_complex_arctanh (t);
}

} // math.c

// pow_int.c ///////////////////////////////////////////////////////////////////
extern "C" {
#ifndef HIDE_INLINE_STATIC
double gsl_pow_2(const double x) { return x*x;   }
double gsl_pow_3(const double x) { return x*x*x; }
double gsl_pow_4(const double x) { double x2 = x*x;   return x2*x2;    }
double gsl_pow_5(const double x) { double x2 = x*x;   return x2*x2*x;  }
double gsl_pow_6(const double x) { double x2 = x*x;   return x2*x2*x2; }
double gsl_pow_7(const double x) { double x3 = x*x*x; return x3*x3*x;  }
double gsl_pow_8(const double x) { double x2 = x*x;   double x4 = x2*x2; return x4*x4; }
double gsl_pow_9(const double x) { double x3 = x*x*x; return x3*x3*x3; }
#endif

double gsl_pow_int(double x, int n)
{
  double value = 1.0;

  if(n < 0) {
    x = 1.0/x;
    n = -n;
  }

  /* repeated squaring method 
   * returns 0.0^0 = 1.0, so continuous in x
   */
  do {
     if(n & 1) value *= x;  /* for n odd */
     n >>= 1;
     x *= x;
  } while (n);

  return value;
}

} // pow_int.c

// psi.c ///////////////////////////////////////////////////////////////////////
extern "C" { 
/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/


/* Chebyshev fit for f(y) = Re(Psi(1+Iy)) + M_EULER - y^2/(1+y^2) - y^2/(2(4+y^2))
 * 1 < y < 10
 *   ==>
 * y(x) = (9x + 11)/2,  -1 < x < 1
 * x(y) = (2y - 11)/9
 *
 * g(x) := f(y(x))
 */
static double r1py_data[] = {
   1.59888328244976954803168395603,
   0.67905625353213463845115658455,
  -0.068485802980122530009506482524,
  -0.005788184183095866792008831182,
   0.008511258167108615980419855648,
  -0.004042656134699693434334556409,
   0.001352328406159402601778462956,
  -0.000311646563930660566674525382,
   0.000018507563785249135437219139,
   0.000028348705427529850296492146,
  -0.000019487536014574535567541960,
   8.0709788710834469408621587335e-06,
  -2.2983564321340518037060346561e-06,
   3.0506629599604749843855962658e-07,
   1.3042238632418364610774284846e-07,
  -1.2308657181048950589464690208e-07,
   5.7710855710682427240667414345e-08,
  -1.8275559342450963966092636354e-08,
   3.1020471300626589420759518930e-09,
   6.8989327480593812470039430640e-10,
  -8.7182290258923059852334818997e-10,
   4.4069147710243611798213548777e-10,
  -1.4727311099198535963467200277e-10,
   2.7589682523262644748825844248e-11,
   4.1871826756975856411554363568e-12,
  -6.5673460487260087541400767340e-12,
   3.4487900886723214020103638000e-12,
  -1.1807251417448690607973794078e-12,
   2.3798314343969589258709315574e-13,
   2.1663630410818831824259465821e-15
};
static cheb_series r1py_cs = {
  r1py_data,
  29,
  -1,1,
  18
};


/* Chebyshev fits from SLATEC code for psi(x)

 Series for PSI        on the interval  0.         to  1.00000D+00
                                       with weighted error   2.03E-17
                                        log weighted error  16.69
                              significant figures required  16.39
                                   decimal places required  17.37

 Series for APSI       on the interval  0.         to  2.50000D-01
                                       with weighted error   5.54E-17
                                        log weighted error  16.26
                              significant figures required  14.42
                                   decimal places required  16.86

*/

static double psics_data[23] = {
  -.038057080835217922,
   .491415393029387130, 
  -.056815747821244730,
   .008357821225914313,
  -.001333232857994342,
   .000220313287069308,
  -.000037040238178456,
   .000006283793654854,
  -.000001071263908506,
   .000000183128394654,
  -.000000031353509361,
   .000000005372808776,
  -.000000000921168141,
   .000000000157981265,
  -.000000000027098646,
   .000000000004648722,
  -.000000000000797527,
   .000000000000136827,
  -.000000000000023475,
   .000000000000004027,
  -.000000000000000691,
   .000000000000000118,
  -.000000000000000020
};
static cheb_series psi_cs = {
  psics_data,
  22,
  -1, 1,
  17
};

static double apsics_data[16] = {    
  -.0204749044678185,
  -.0101801271534859,
   .0000559718725387,
  -.0000012917176570,
   .0000000572858606,
  -.0000000038213539,
   .0000000003397434,
  -.0000000000374838,
   .0000000000048990,
  -.0000000000007344,
   .0000000000001233,
  -.0000000000000228,
   .0000000000000045,
  -.0000000000000009,
   .0000000000000002,
  -.0000000000000000 
};    
static cheb_series apsi_cs = {
  apsics_data,
  15,
  -1, 1,
  9
};

#define PSI_TABLE_NMAX 100
static double psi_table[PSI_TABLE_NMAX+1] = {
  0.0,  /* Infinity */              /* psi(0) */
 -M_EULER,                          /* psi(1) */
  0.42278433509846713939348790992,  /* ...    */
  0.92278433509846713939348790992,
  1.25611766843180047272682124325,
  1.50611766843180047272682124325,
  1.70611766843180047272682124325,
  1.87278433509846713939348790992,
  2.01564147795560999653634505277,
  2.14064147795560999653634505277,
  2.25175258906672110764745616389,
  2.35175258906672110764745616389,
  2.44266167997581201673836525479,
  2.52599501330914535007169858813,
  2.60291809023222227314862166505,
  2.67434666166079370172005023648,
  2.74101332832746036838671690315,
  2.80351332832746036838671690315,
  2.86233685773922507426906984432,
  2.91789241329478062982462539988,
  2.97052399224214905087725697883,
  3.02052399224214905087725697883,
  3.06814303986119666992487602645,
  3.11359758531574212447033057190,
  3.15707584618530734186163491973,
  3.1987425128519740085283015864,
  3.2387425128519740085283015864,
  3.2772040513135124700667631249,
  3.3142410883505495071038001619,
  3.3499553740648352213895144476,
  3.3844381326855248765619282407,
  3.4177714660188582098952615740,
  3.4500295305349872421533260902,
  3.4812795305349872421533260902,
  3.5115825608380175451836291205,
  3.5409943255438998981248055911,
  3.5695657541153284695533770196,
  3.5973435318931062473311547974,
  3.6243705589201332743581818244,
  3.6506863483938174848844976139,
  3.6763273740348431259101386396,
  3.7013273740348431259101386396,
  3.7257176179372821503003825420,
  3.7495271417468059598241920658,
  3.7727829557002943319172153216,
  3.7955102284275670591899425943,
  3.8177324506497892814121648166,
  3.8394715810845718901078169905,
  3.8607481768292527411716467777,
  3.8815815101625860745049801110,
  3.9019896734278921969539597029,
  3.9219896734278921969539597029,
  3.9415975165651470989147440166,
  3.9608282857959163296839747858,
  3.9796962103242182164764276160,
  3.9982147288427367349949461345,
  4.0163965470245549168131279527,
  4.0342536898816977739559850956,
  4.0517975495308205809735289552,
  4.0690389288411654085597358518,
  4.0859880813835382899156680552,
  4.1026547480502049565823347218,
  4.1190481906731557762544658694,
  4.1351772229312202923834981274,
  4.1510502388042361653993711433,
  4.1666752388042361653993711433,
  4.1820598541888515500147557587,
  4.1972113693403667015299072739,
  4.2121367424746950597388624977,
  4.2268426248276362362094507330,
  4.2413353784508246420065521823,
  4.2556210927365389277208378966,
  4.2697055997787924488475984600,
  4.2835944886676813377364873489,
  4.2972931188046676391063503626,
  4.3108066323181811526198638761,
  4.3241399656515144859531972094,
  4.3372978603883565912163551041,
  4.3502848733753695782293421171,
  4.3631053861958823987421626300,
  4.3757636140439836645649474401,
  4.3882636140439836645649474401,
  4.4006092930563293435772931191,
  4.4128044150075488557724150703,
  4.4248526077786331931218126607,
  4.4367573696833950978837174226,
  4.4485220755657480390601880108,
  4.4601499825424922251066996387,
  4.4716442354160554434975042364,
  4.4830078717796918071338678728,
  4.4942438268358715824147667492,
  4.5053549379469826935258778603,
  4.5163439489359936825368668713,
  4.5272135141533849868846929582,
  4.5379662023254279976373811303,
  4.5486045001977684231692960239,
  4.5591308159872421073798223397,
  4.5695474826539087740464890064,
  4.5798567610044242379640147796,
  4.5900608426370772991885045755,
  4.6001618527380874001986055856
};


#define PSI_1_TABLE_NMAX 100
static double psi_1_table[PSI_1_TABLE_NMAX+1] = {
  0.0,  /* Infinity */              /* psi(1,0) */
  M_PI*M_PI/6.0,                    /* psi(1,1) */
  0.644934066848226436472415,       /* ...      */
  0.394934066848226436472415,
  0.2838229557371153253613041,
  0.2213229557371153253613041,
  0.1813229557371153253613041,
  0.1535451779593375475835263,
  0.1331370146940314251345467,
  0.1175120146940314251345467,
  0.1051663356816857461222010,
  0.0951663356816857461222010,
  0.0869018728717683907503002,
  0.0799574284273239463058557,
  0.0740402686640103368384001,
  0.0689382278476838062261552,
  0.0644937834032393617817108,
  0.0605875334032393617817108,
  0.0571273257907826143768665,
  0.0540409060376961946237801,
  0.0512708229352031198315363,
  0.0487708229352031198315363,
  0.0465032492390579951149830,
  0.0444371335365786562720078,
  0.0425467743683366902984728,
  0.0408106632572255791873617,
  0.0392106632572255791873617,
  0.0377313733163971768204978,
  0.0363596312039143235969038,
  0.0350841209998326909438426,
  0.0338950603577399442137594,
  0.0327839492466288331026483,
  0.0317433665203020901265817,
  0.03076680402030209012658168,
  0.02984853037475571730748159,
  0.02898347847164153045627052,
  0.02816715194102928555831133,
  0.02739554700275768062003973,
  0.02666508681283803124093089,
  0.02597256603721476254286995,
  0.02531510384129102815759710,
  0.02469010384129102815759710,
  0.02409521984367056414807896,
  0.02352832641963428296894063,
  0.02298749353699501850166102,
  0.02247096461137518379091722,
  0.02197713745088135663042339,
  0.02150454765882086513703965,
  0.02105185413233829383780923,
  0.02061782635456051606003145,
  0.02020133322669712580597065,
  0.01980133322669712580597065,
  0.01941686571420193164987683,
  0.01904704322899483105816086,
  0.01869104465298913508094477,
  0.01834810912486842177504628,
  0.01801753061247172756017024,
  0.01769865306145131939690494,
  0.01739086605006319997554452,
  0.01709360088954001329302371,
  0.01680632711763538818529605,
  0.01652854933985761040751827,
  0.01625980437882562975715546,
  0.01599965869724394401313881,
  0.01574770606433893015574400,
  0.01550356543933893015574400,
  0.01526687904880638577704578,
  0.01503731063741979257227076,
  0.01481454387422086185273411,
  0.01459828089844231513993134,
  0.01438824099085987447620523,
  0.01418415935820681325171544,
  0.01398578601958352422176106,
  0.01379288478501562298719316,
  0.01360523231738567365335942,
  0.01342261726990576130858221,
  0.01324483949212798353080444,
  0.01307170929822216635628920,
  0.01290304679189732236910755,
  0.01273868124291638877278934,
  0.01257845051066194236996928,
  0.01242220051066194236996928,
  0.01226978472038606978956995,
  0.01212106372098095378719041,
  0.01197590477193174490346273,
  0.01183418141592267460867815,
  0.01169577311142440471248438,
  0.01156056489076458859566448,
  0.01142844704164317229232189,
  0.01129931481023821361463594,
  0.01117306812421372175754719,
  0.01104961133409026496742374,
  0.01092885297157366069257770,
  0.01081070552355853781923177,
  0.01069508522063334415522437,
  0.01058191183901270133041676,
  0.01047110851491297833872701,
  0.01036260157046853389428257,
  0.01025632035036012704977199,  /* ...        */
  0.01015219706839427948625679,  /* psi(1,99)  */
  0.01005016666333357139524567   /* psi(1,100) */
};


/* digamma for x both positive and negative; we do both
 * cases here because of the way we use even/odd parts
 * of the function
 */
static int
psi_x(const double x, gsl_sf_result * result)
{
  const double y = fabs(x);

  if(x == 0.0 || x == -1.0 || x == -2.0) {
    DOMAIN_ERROR(result);
  }
  else if(y >= 2.0) {
    const double t = 8.0/(y*y)-1.0;
    gsl_sf_result result_c;
    cheb_eval_e(&apsi_cs, t, &result_c);
    if(x < 0.0) {
      const double s = sin(M_PI*x);
      const double c = cos(M_PI*x);
      if(fabs(s) < 2.0*GSL_SQRT_DBL_MIN) {
        DOMAIN_ERROR(result);
      }
      else {
        result->val  = log(y) - 0.5/x + result_c.val - M_PI * c/s;
        result->err  = M_PI*fabs(x)*GSL_DBL_EPSILON/(s*s);
        result->err += result_c.err;
        result->err += GSL_DBL_EPSILON * fabs(result->val);
        return GSL_SUCCESS;
      }
    }
    else {
      result->val  = log(y) - 0.5/x + result_c.val;
      result->err  = result_c.err;
      result->err += GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
  }
  else { /* -2 < x < 2 */
    gsl_sf_result result_c;

    if(x < -1.0) { /* x = -2 + v */
      const double v  = x + 2.0;
      const double t1 = 1.0/x;
      const double t2 = 1.0/(x+1.0);
      const double t3 = 1.0/v;
      cheb_eval_e(&psi_cs, 2.0*v-1.0, &result_c);
      
      result->val  = -(t1 + t2 + t3) + result_c.val;
      result->err  = GSL_DBL_EPSILON * (fabs(t1) + fabs(x/(t2*t2)) + fabs(x/(t3*t3)));
      result->err += result_c.err;
      result->err += GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if(x < 0.0) { /* x = -1 + v */
      const double v  = x + 1.0;
      const double t1 = 1.0/x;
      const double t2 = 1.0/v;
      cheb_eval_e(&psi_cs, 2.0*v-1.0, &result_c);
      
      result->val  = -(t1 + t2) + result_c.val;
      result->err  = GSL_DBL_EPSILON * (fabs(t1) + fabs(x/(t2*t2)));
      result->err += result_c.err;
      result->err += GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if(x < 1.0) { /* x = v */
      const double t1 = 1.0/x;
      cheb_eval_e(&psi_cs, 2.0*x-1.0, &result_c);
      
      result->val  = -t1 + result_c.val;
      result->err  = GSL_DBL_EPSILON * t1;
      result->err += result_c.err;
      result->err += GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else { /* x = 1 + v */
      const double v = x - 1.0;
      return cheb_eval_e(&psi_cs, 2.0*v-1.0, result);
    }
  }
}


/* psi(z) for large |z| in the right half-plane; [Abramowitz + Stegun, 6.3.18] */
static
gsl_complex
psi_complex_asymp(gsl_complex z)
{
  /* coefficients in the asymptotic expansion for large z;
   * let w = z^(-2) and write the expression in the form
   *
   *   ln(z) - 1/(2z) - 1/12 w (1 + c1 w + c2 w + c3 w + ... )
   */
  static const double c1 = -0.1;
  static const double c2 =  1.0/21.0;
  static const double c3 = -0.05;

  gsl_complex zi = gsl_complex_inverse(z);
  gsl_complex w  = gsl_complex_mul(zi, zi);
  gsl_complex cs;

  /* Horner method evaluation of term in parentheses */
  gsl_complex sum;
  sum = gsl_complex_mul_real(w, c3/c2);
  sum = gsl_complex_add_real(sum, 1.0);
  sum = gsl_complex_mul_real(sum, c2/c1);
  sum = gsl_complex_mul(sum, w);
  sum = gsl_complex_add_real(sum, 1.0);
  sum = gsl_complex_mul_real(sum, c1);
  sum = gsl_complex_mul(sum, w);
  sum = gsl_complex_add_real(sum, 1.0);

  /* correction added to log(z) */
  cs = gsl_complex_mul(sum, w);
  cs = gsl_complex_mul_real(cs, -1.0/12.0);
  cs = gsl_complex_add(cs, gsl_complex_mul_real(zi, -0.5));

  return gsl_complex_add(gsl_complex_log(z), cs);
}



/* psi(z) for complex z in the right half-plane */
static int
psi_complex_rhp(
  gsl_complex z,
  gsl_sf_result * result_re,
  gsl_sf_result * result_im
  )
{
  int n_recurse = 0;
  int i;
  gsl_complex a;

  if(GSL_REAL(z) == 0.0 && GSL_IMAG(z) == 0.0)
  {
    result_re->val = 0.0;
    result_im->val = 0.0;
    result_re->err = 0.0;
    result_im->err = 0.0;
    return GSL_EDOM;
  }

  /* compute the number of recurrences to apply */
  if(GSL_REAL(z) < 20.0 && fabs(GSL_IMAG(z)) < 20.0)
  {
    const double sp = sqrt(20.0 + GSL_IMAG(z));
    const double sn = sqrt(20.0 - GSL_IMAG(z));
    const double rhs = sp*sn - GSL_REAL(z);
    if(rhs > 0.0) n_recurse = ceil(rhs);
  }

  /* compute asymptotic at the large value z + n_recurse */
  a = psi_complex_asymp(gsl_complex_add_real(z, n_recurse));

  /* descend recursively, if necessary */
  for(i = n_recurse; i >= 1; --i)
  {
    gsl_complex zn = gsl_complex_add_real(z, i - 1.0);
    gsl_complex zn_inverse = gsl_complex_inverse(zn);
    a = gsl_complex_sub(a, zn_inverse);
  }

  result_re->val = GSL_REAL(a);
  result_im->val = GSL_IMAG(a);
  result_re->err = 2.0 * (1.0 + n_recurse) * GSL_DBL_EPSILON * fabs(result_re->val);
  result_im->err = 2.0 * (1.0 + n_recurse) * GSL_DBL_EPSILON * fabs(result_im->val);

  return GSL_SUCCESS;
}



/* generic polygamma; assumes n >= 0 and x > 0
 */
static int
psi_n_xg0(const int n, const double x, gsl_sf_result * result)
{
  if(n == 0) {
    return gsl_sf_psi_e(x, result);
  }
  else {
    /* Abramowitz + Stegun 6.4.10 */
    gsl_sf_result ln_nf;
    gsl_sf_result hzeta;
    int stat_hz = gsl_sf_hzeta_e(n+1.0, x, &hzeta);
    int stat_nf = gsl_sf_lnfact_e((unsigned int) n, &ln_nf);
    int stat_e  = gsl_sf_exp_mult_err_e(ln_nf.val, ln_nf.err,
                                           hzeta.val, hzeta.err,
                                           result);
    if(GSL_IS_EVEN(n)) result->val = -result->val;
    return GSL_ERROR_SELECT_3(stat_e, stat_nf, stat_hz);
  }
}



/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_psi_int_e(const int n, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(n <= 0) {
    DOMAIN_ERROR(result);
  }
  else if(n <= PSI_TABLE_NMAX) {
    result->val = psi_table[n];
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    /* Abramowitz+Stegun 6.3.18 */
    const double c2 = -1.0/12.0;
    const double c3 =  1.0/120.0;
    const double c4 = -1.0/252.0;
    const double c5 =  1.0/240.0;
    const double ni2 = (1.0/n)*(1.0/n);
    const double ser = ni2 * (c2 + ni2 * (c3 + ni2 * (c4 + ni2*c5)));
    result->val  = log(n) - 0.5/n + ser;
    result->err  = GSL_DBL_EPSILON * (fabs(log(n)) + fabs(0.5/n) + fabs(ser));
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int gsl_sf_psi_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */
  return psi_x(x, result);
}


int
gsl_sf_psi_1piy_e(const double y, gsl_sf_result * result)
{
  const double ay = fabs(y);

  /* CHECK_POINTER(result) */

  if(ay > 1000.0) {
    /* [Abramowitz+Stegun, 6.3.19] */
    const double yi2 = 1.0/(ay*ay);
    const double lny = log(ay);
    const double sum = yi2 * (1.0/12.0 + 1.0/120.0 * yi2 + 1.0/252.0 * yi2*yi2);
    result->val = lny + sum;
    result->err = 2.0 * GSL_DBL_EPSILON * (fabs(lny) + fabs(sum));
    return GSL_SUCCESS;
  }
  else if(ay > 10.0) {
    /* [Abramowitz+Stegun, 6.3.19] */
    const double yi2 = 1.0/(ay*ay);
    const double lny = log(ay);
    const double sum = yi2 * (1.0/12.0 +
                         yi2 * (1.0/120.0 +
                           yi2 * (1.0/252.0 +
                             yi2 * (1.0/240.0 +
                               yi2 * (1.0/132.0 + 691.0/32760.0 * yi2)))));
    result->val = lny + sum;
    result->err = 2.0 * GSL_DBL_EPSILON * (fabs(lny) + fabs(sum));
    return GSL_SUCCESS;
  }
  else if(ay > 1.0){
    const double y2 = ay*ay;
    const double x  = (2.0*ay - 11.0)/9.0;
    const double v  = y2*(1.0/(1.0+y2) + 0.5/(4.0+y2));
    gsl_sf_result result_c;
    cheb_eval_e(&r1py_cs, x, &result_c);
    result->val  = result_c.val - M_EULER + v;
    result->err  = result_c.err;
    result->err += 2.0 * GSL_DBL_EPSILON * (fabs(v) + M_EULER + fabs(result_c.val));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    result->err *= 5.0; /* FIXME: losing a digit somewhere... maybe at x=... ? */
    return GSL_SUCCESS;
  }
  else {
    /* [Abramowitz+Stegun, 6.3.17]
     *
     * -M_EULER + y^2 Sum[1/n 1/(n^2 + y^2), {n,1,M}]
     *   +     Sum[1/n^3, {n,M+1,Infinity}]
     *   - y^2 Sum[1/n^5, {n,M+1,Infinity}]
     *   + y^4 Sum[1/n^7, {n,M+1,Infinity}]
     *   - y^6 Sum[1/n^9, {n,M+1,Infinity}]
     *   + O(y^8)
     *
     * We take M=50 for at least 15 digit precision.
     */
    const int M = 50;
    const double y2 = y*y;
    const double c0 = 0.00019603999466879846570;
    const double c2 = 3.8426659205114376860e-08;
    const double c4 = 1.0041592839497643554e-11;
    const double c6 = 2.9516743763500191289e-15;
    const double p  = c0 + y2 *(-c2 + y2*(c4 - y2*c6));
    double sum = 0.0;
    double v;
    
    int n;
    for(n=1; n<=M; n++) {
      sum += 1.0/(n * (n*n + y*y));
    }

    v = y2 * (sum + p);
    result->val  = -M_EULER + v;
    result->err  = GSL_DBL_EPSILON * (M_EULER + fabs(v));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int gsl_sf_psi_1_int_e(const int n, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */
  if(n <= 0) {
    DOMAIN_ERROR(result);
  }
  else if(n <= PSI_1_TABLE_NMAX) {
    result->val = psi_1_table[n];
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
  else {
    /* Abramowitz+Stegun 6.4.12
     * double-precision for n > 100
     */
    const double c0 = -1.0/30.0;
    const double c1 =  1.0/42.0;
    const double c2 = -1.0/30.0;
    const double ni2 = (1.0/n)*(1.0/n);
    const double ser =  ni2*ni2 * (c0 + ni2*(c1 + c2*ni2));
    result->val = (1.0 + 0.5/n + 1.0/(6.0*n*n) + ser) / n;
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
}


int gsl_sf_psi_1_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x == 0.0 || x == -1.0 || x == -2.0) {
    DOMAIN_ERROR(result);
  }
  else if(x > 0.0)
  {
    return psi_n_xg0(1, x, result);
  }
  else if(x > -5.0)
  {
    /* Abramowitz + Stegun 6.4.6 */
    int M = -floor(x);
    double fx = x + M;
    double sum = 0.0;
    int m;

    if(fx == 0.0)
      DOMAIN_ERROR(result);

    for(m = 0; m < M; ++m)
      sum += 1.0/((x+m)*(x+m));

    {
      int stat_psi = psi_n_xg0(1, fx, result);
      result->val += sum;
      result->err += M * GSL_DBL_EPSILON * sum;
      return stat_psi;
    }
  }
  else
  {
    /* Abramowitz + Stegun 6.4.7 */
    const double sin_px = sin(M_PI * x);
    const double d = M_PI*M_PI/(sin_px*sin_px);
    gsl_sf_result r;
    int stat_psi = psi_n_xg0(1, 1.0-x, &r);
    result->val = d - r.val;
    result->err = r.err + 2.0*GSL_DBL_EPSILON*d;
    return stat_psi;
  }
}


int gsl_sf_psi_n_e(const int n, const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(n == 0)
  {
    return gsl_sf_psi_e(x, result);
  }
  else if(n == 1)
  {
    return gsl_sf_psi_1_e(x, result);
  }
  else if(n < 0 || x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else {
    gsl_sf_result ln_nf;
    gsl_sf_result hzeta;
    int stat_hz = gsl_sf_hzeta_e(n+1.0, x, &hzeta);
    int stat_nf = gsl_sf_lnfact_e((unsigned int) n, &ln_nf);
    int stat_e  = gsl_sf_exp_mult_err_e(ln_nf.val, ln_nf.err,
                                           hzeta.val, hzeta.err,
                                           result);
    if(GSL_IS_EVEN(n)) result->val = -result->val;
    return GSL_ERROR_SELECT_3(stat_e, stat_nf, stat_hz);
  }
}


int
gsl_sf_complex_psi_e(
  const double x,
  const double y,
  gsl_sf_result * result_re,
  gsl_sf_result * result_im
  )
{
  if(x >= 0.0)
  {
    gsl_complex z = gsl_complex_rect(x, y);
    return psi_complex_rhp(z, result_re, result_im);
  }
  else
  {
    /* reflection formula [Abramowitz+Stegun, 6.3.7] */
    gsl_complex z = gsl_complex_rect(x, y);
    gsl_complex omz = gsl_complex_rect(1.0 - x, -y);
    gsl_complex zpi = gsl_complex_mul_real(z, M_PI);
    gsl_complex cotzpi = gsl_complex_cot(zpi);
    int ret_val = psi_complex_rhp(omz, result_re, result_im);

    if(GSL_IS_REAL(GSL_REAL(cotzpi)) && GSL_IS_REAL(GSL_IMAG(cotzpi)))
    {
      result_re->val -= M_PI * GSL_REAL(cotzpi);
      result_im->val -= M_PI * GSL_IMAG(cotzpi);
      return ret_val;
    }
    else
    {
      GSL_ERROR("singularity", GSL_EDOM);
    }
  }
}



/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_psi_int(const int n)
{
  EVAL_RESULT(gsl_sf_psi_int_e(n, &result));
}

double gsl_sf_psi(const double x)
{
  EVAL_RESULT(gsl_sf_psi_e(x, &result));
}

double gsl_sf_psi_1piy(const double x)
{
  EVAL_RESULT(gsl_sf_psi_1piy_e(x, &result));
}

double gsl_sf_psi_1_int(const int n)
{
  EVAL_RESULT(gsl_sf_psi_1_int_e(n, &result));
}

double gsl_sf_psi_1(const double x)
{
  EVAL_RESULT(gsl_sf_psi_1_e(x, &result));
}

double gsl_sf_psi_n(const int n, const double x)
{
  EVAL_RESULT(gsl_sf_psi_n_e(n, x, &result));
}
} // psi.c

// stream.c ///////////////////////////////////////////////////////////////////////
extern "C" { 
FILE * gsl_stream = NULL ;
gsl_stream_handler_t * gsl_stream_handler = NULL;

void
gsl_stream_printf (const char *label, const char *file, int line, 
                   const char *reason)
{
  if (gsl_stream == NULL)
    {
      gsl_stream = stderr;
    }
  if (gsl_stream_handler)
    {
      (*gsl_stream_handler) (label, file, line, reason);
      return;
    }
  fprintf (gsl_stream, "gsl: %s:%d: %s: %s\n", file, line, label, reason);

}

gsl_stream_handler_t *
gsl_set_stream_handler (gsl_stream_handler_t * new_handler)
{
  gsl_stream_handler_t * previous_handler = gsl_stream_handler;
  gsl_stream_handler = new_handler;
  return previous_handler;
}

FILE *
gsl_set_stream (FILE * new_stream)
{
  FILE * previous_stream;
  if (gsl_stream == NULL) {
    gsl_stream = stderr;
  }
  previous_stream = gsl_stream;
  gsl_stream = new_stream;
  return previous_stream;
}
} // stream.c

// trig.c ///////////////////////////////////////////////////////////////////////
extern "C" { 
/* sinh(x) series
 * double-precision for |x| < 1.0
 */
inline
static
int
sinh_series(const double x, double * result)
{
  const double y = x*x;
  const double c0 = 1.0/6.0;
  const double c1 = 1.0/120.0;
  const double c2 = 1.0/5040.0;
  const double c3 = 1.0/362880.0;
  const double c4 = 1.0/39916800.0;
  const double c5 = 1.0/6227020800.0;
  const double c6 = 1.0/1307674368000.0;
  const double c7 = 1.0/355687428096000.0;
  *result = x*(1.0 + y*(c0+y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*c7))))))));
  return GSL_SUCCESS;
}


/* cosh(x)-1 series
 * double-precision for |x| < 1.0
 */
inline
static
int
cosh_m1_series(const double x, double * result)
{
  const double y = x*x;
  const double c0 = 0.5;
  const double c1 = 1.0/24.0;
  const double c2 = 1.0/720.0;
  const double c3 = 1.0/40320.0;
  const double c4 = 1.0/3628800.0;
  const double c5 = 1.0/479001600.0;
  const double c6 = 1.0/87178291200.0;
  const double c7 = 1.0/20922789888000.0;
  const double c8 = 1.0/6402373705728000.0;
  *result = y*(c0+y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*(c7+y*c8))))))));
  return GSL_SUCCESS;
}


/* Chebyshev expansion for f(t) = sinc((t+1)/2), -1 < t < 1
 */
static double sinc_data[17] = {
  1.133648177811747875422,
 -0.532677564732557348781,
 -0.068293048346633177859,
  0.033403684226353715020,
  0.001485679893925747818,
 -0.000734421305768455295,
 -0.000016837282388837229,
  0.000008359950146618018,
  0.000000117382095601192,
 -0.000000058413665922724,
 -0.000000000554763755743,
  0.000000000276434190426,
  0.000000000001895374892,
 -0.000000000000945237101,
 -0.000000000000004900690,
  0.000000000000002445383,
  0.000000000000000009925
};
static cheb_series sinc_cs = {
  sinc_data,
  16,
  -1, 1,
  10
};


/* Chebyshev expansion for f(t) = g((t+1)Pi/8), -1<t<1
 * g(x) = (sin(x)/x - 1)/(x*x)
 */
static double sin_data[12] = {
  -0.3295190160663511504173,
   0.0025374284671667991990,
   0.0006261928782647355874,
  -4.6495547521854042157541e-06,
  -5.6917531549379706526677e-07,
   3.7283335140973803627866e-09,
   3.0267376484747473727186e-10,
  -1.7400875016436622322022e-12,
  -1.0554678305790849834462e-13,
   5.3701981409132410797062e-16,
   2.5984137983099020336115e-17,
  -1.1821555255364833468288e-19
};
static cheb_series sin_cs = {
  sin_data,
  11,
  -1, 1,
  11
};

/* Chebyshev expansion for f(t) = g((t+1)Pi/8), -1<t<1
 * g(x) = (2(cos(x) - 1)/(x^2) + 1) / x^2
 */
static double cos_data[11] = {
  0.165391825637921473505668118136,
 -0.00084852883845000173671196530195,
 -0.000210086507222940730213625768083,
  1.16582269619760204299639757584e-6,
  1.43319375856259870334412701165e-7,
 -7.4770883429007141617951330184e-10,
 -6.0969994944584252706997438007e-11,
  2.90748249201909353949854872638e-13,
  1.77126739876261435667156490461e-14,
 -7.6896421502815579078577263149e-17,
 -3.7363121133079412079201377318e-18
};
static cheb_series cos_cs = {
  cos_data,
  10,
  -1, 1,
  10
};


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

/* I would have prefered just using the library sin() function.
 * But after some experimentation I decided that there was
 * no good way to understand the error; library sin() is just a black box.
 * So we have to roll our own.
 */
int
gsl_sf_sin_e(double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  {
    const double P1 = 7.85398125648498535156e-1;
    const double P2 = 3.77489470793079817668e-8;
    const double P3 = 2.69515142907905952645e-15;

    const double sgn_x = GSL_SIGN(x);
    const double abs_x = fabs(x);

    if(abs_x < GSL_ROOT4_DBL_EPSILON) {
      const double x2 = x*x;
      result->val = x * (1.0 - x2/6.0);
      result->err = fabs(x*x2*x2 / 100.0);
      return GSL_SUCCESS;
    }
    else {
      double sgn_result = sgn_x;
      double y = floor(abs_x/(0.25*M_PI));
      int octant = y - ldexp(floor(ldexp(y,-3)),3);
      int stat_cs;
      double z;

      if(GSL_IS_ODD(octant)) {
        octant += 1;
        octant &= 07;
        y += 1.0;
      }

      if(octant > 3) {
        octant -= 4;
        sgn_result = -sgn_result;
      }
      
      z = ((abs_x - y * P1) - y * P2) - y * P3;

      if(octant == 0) {
        gsl_sf_result sin_cs_result;
        const double t = 8.0*fabs(z)/M_PI - 1.0;
        stat_cs = cheb_eval_e(&sin_cs, t, &sin_cs_result);
        result->val = z * (1.0 + z*z * sin_cs_result.val);
      }
      else { /* octant == 2 */
        gsl_sf_result cos_cs_result;
        const double t = 8.0*fabs(z)/M_PI - 1.0;
        stat_cs = cheb_eval_e(&cos_cs, t, &cos_cs_result);
        result->val = 1.0 - 0.5*z*z * (1.0 - z*z * cos_cs_result.val);
      }

      result->val *= sgn_result;

      if(abs_x > 1.0/GSL_DBL_EPSILON) {
        result->err = fabs(result->val);
      }
      else if(abs_x > 100.0/GSL_SQRT_DBL_EPSILON) {
        result->err = 2.0 * abs_x * GSL_DBL_EPSILON * fabs(result->val);
      }
      else if(abs_x > 0.1/GSL_SQRT_DBL_EPSILON) {
        result->err = 2.0 * GSL_SQRT_DBL_EPSILON * fabs(result->val);
      }
      else {
        result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      }

      return stat_cs;
    }
  }
}


int
gsl_sf_cos_e(double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  {
    const double P1 = 7.85398125648498535156e-1;
    const double P2 = 3.77489470793079817668e-8;
    const double P3 = 2.69515142907905952645e-15;

    const double abs_x = fabs(x);

    if(abs_x < GSL_ROOT4_DBL_EPSILON) {
      const double x2 = x*x;
      result->val = 1.0 - 0.5*x2;
      result->err = fabs(x2*x2/12.0);
      return GSL_SUCCESS;
    }
    else {
      double sgn_result = 1.0;
      double y = floor(abs_x/(0.25*M_PI));
      int octant = y - ldexp(floor(ldexp(y,-3)),3);
      int stat_cs;
      double z;

      if(GSL_IS_ODD(octant)) {
        octant += 1;
        octant &= 07;
        y += 1.0;
      }

      if(octant > 3) {
        octant -= 4;
        sgn_result = -sgn_result;
      }

      if(octant > 1) {
        sgn_result = -sgn_result;
      }

      z = ((abs_x - y * P1) - y * P2) - y * P3;

      if(octant == 0) {
        gsl_sf_result cos_cs_result;
        const double t = 8.0*fabs(z)/M_PI - 1.0;
        stat_cs = cheb_eval_e(&cos_cs, t, &cos_cs_result);
        result->val = 1.0 - 0.5*z*z * (1.0 - z*z * cos_cs_result.val);
      }
      else { /* octant == 2 */
        gsl_sf_result sin_cs_result;
        const double t = 8.0*fabs(z)/M_PI - 1.0;
        stat_cs = cheb_eval_e(&sin_cs, t, &sin_cs_result);
        result->val = z * (1.0 + z*z * sin_cs_result.val);
      }

      result->val *= sgn_result;

      if(abs_x > 1.0/GSL_DBL_EPSILON) {
        result->err = fabs(result->val);
      }
      else if(abs_x > 100.0/GSL_SQRT_DBL_EPSILON) {
        result->err = 2.0 * abs_x * GSL_DBL_EPSILON * fabs(result->val);
      }
      else if(abs_x > 0.1/GSL_SQRT_DBL_EPSILON) {
        result->err = 2.0 * GSL_SQRT_DBL_EPSILON * fabs(result->val);
      }
      else {
        result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      }

      return stat_cs;
    }
  }
}


int
gsl_sf_hypot_e(const double x, const double y, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x == 0.0 && y == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    const double a = fabs(x);
    const double b = fabs(y);
    const double min = GSL_MIN_DBL(a,b);
    const double max = GSL_MAX_DBL(a,b);
    const double rat = min/max;
    const double root_term = sqrt(1.0 + rat*rat);

    if(max < GSL_DBL_MAX/root_term) {
      result->val = max * root_term;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      OVERFLOW_ERROR(result);
    }
  }
}


int
gsl_sf_complex_sin_e(const double zr, const double zi,
                        gsl_sf_result * szr, gsl_sf_result * szi)
{
  /* CHECK_POINTER(szr) */
  /* CHECK_POINTER(szi) */

  if(fabs(zi) < 1.0) {
    double ch_m1, sh;
    sinh_series(zi, &sh);
    cosh_m1_series(zi, &ch_m1);
    szr->val = sin(zr)*(ch_m1 + 1.0);
    szi->val = cos(zr)*sh;
    szr->err = 2.0 * GSL_DBL_EPSILON * fabs(szr->val);
    szi->err = 2.0 * GSL_DBL_EPSILON * fabs(szi->val);
    return GSL_SUCCESS;
  }
  else if(fabs(zi) < GSL_LOG_DBL_MAX) {
    double ex = exp(zi);
    double ch = 0.5*(ex+1.0/ex);
    double sh = 0.5*(ex-1.0/ex);
    szr->val = sin(zr)*ch;
    szi->val = cos(zr)*sh;
    szr->err = 2.0 * GSL_DBL_EPSILON * fabs(szr->val);
    szi->err = 2.0 * GSL_DBL_EPSILON * fabs(szi->val);
    return GSL_SUCCESS;
  }
  else {
    OVERFLOW_ERROR_2(szr, szi);
  }
}


int
gsl_sf_complex_cos_e(const double zr, const double zi,
                        gsl_sf_result * czr, gsl_sf_result * czi)
{
  /* CHECK_POINTER(czr) */
  /* CHECK_POINTER(czi) */

  if(fabs(zi) < 1.0) {
    double ch_m1, sh;
    sinh_series(zi, &sh);
    cosh_m1_series(zi, &ch_m1);
    czr->val =  cos(zr)*(ch_m1 + 1.0);
    czi->val = -sin(zr)*sh;
    czr->err = 2.0 * GSL_DBL_EPSILON * fabs(czr->val);
    czi->err = 2.0 * GSL_DBL_EPSILON * fabs(czi->val);
    return GSL_SUCCESS;
  }
  else if(fabs(zi) < GSL_LOG_DBL_MAX) {
    double ex = exp(zi);
    double ch = 0.5*(ex+1.0/ex);
    double sh = 0.5*(ex-1.0/ex);
    czr->val =  cos(zr)*ch;
    czi->val = -sin(zr)*sh;
    czr->err = 2.0 * GSL_DBL_EPSILON * fabs(czr->val);
    czi->err = 2.0 * GSL_DBL_EPSILON * fabs(czi->val);
    return GSL_SUCCESS;
  }
  else {
    OVERFLOW_ERROR_2(czr,czi);
  }
}


int
gsl_sf_complex_logsin_e(const double zr, const double zi,
                           gsl_sf_result * lszr, gsl_sf_result * lszi)
{
  /* CHECK_POINTER(lszr) */
  /* CHECK_POINTER(lszi) */

  if(zi > 60.0) {
    lszr->val = -M_LN2 + zi;
    lszi->val =  0.5*M_PI - zr;
    lszr->err = 2.0 * GSL_DBL_EPSILON * fabs(lszr->val);
    lszi->err = 2.0 * GSL_DBL_EPSILON * fabs(lszi->val);
  }
  else if(zi < -60.0) {
    lszr->val = -M_LN2 - zi;
    lszi->val = -0.5*M_PI + zr;
    lszr->err = 2.0 * GSL_DBL_EPSILON * fabs(lszr->val);
    lszi->err = 2.0 * GSL_DBL_EPSILON * fabs(lszi->val);
  }
  else {
    gsl_sf_result sin_r, sin_i;
    int status;
    gsl_sf_complex_sin_e(zr, zi, &sin_r, &sin_i); /* ok by construction */
    status = gsl_sf_complex_log_e(sin_r.val, sin_i.val, lszr, lszi);
    if(status == GSL_EDOM) {
      DOMAIN_ERROR_2(lszr, lszi);
    }
  }
  return gsl_sf_angle_restrict_symm_e(&(lszi->val));
}


int
gsl_sf_lnsinh_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(fabs(x) < 1.0) {
    double eps;
    sinh_series(x, &eps);
    result->val = log(eps);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < -0.5*GSL_LOG_DBL_EPSILON) {
    result->val = x + log(0.5*(1.0 - exp(-2.0*x)));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = -M_LN2 + x;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int gsl_sf_lncosh_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(fabs(x) < 1.0) {
    double eps;
    cosh_m1_series(x, &eps);
    return gsl_sf_log_1plusx_e(eps, result);
  }
  else if(x < -0.5*GSL_LOG_DBL_EPSILON) {
    result->val = x + log(0.5*(1.0 + exp(-2.0*x)));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = -M_LN2 + x;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


/*
inline int gsl_sf_sincos_e(const double theta, double * s, double * c)
{
  double tan_half = tan(0.5 * theta);
  double den = 1. + tan_half*tan_half;
  double cos_theta = (1.0 - tan_half*tan_half) / den;
  double sin_theta = 2.0 * tan_half / den;
}
*/

int
gsl_sf_polar_to_rect(const double r, const double theta,
                          gsl_sf_result * x, gsl_sf_result * y)
{
  double t   = theta;
  int status = gsl_sf_angle_restrict_symm_e(&t);
  double c = cos(t);
  double s = sin(t);
  x->val = r * cos(t);
  y->val = r * sin(t);
  x->err  = r * fabs(s * GSL_DBL_EPSILON * t);
  x->err += 2.0 * GSL_DBL_EPSILON * fabs(x->val);
  y->err  = r * fabs(c * GSL_DBL_EPSILON * t);
  y->err += 2.0 * GSL_DBL_EPSILON * fabs(y->val);
  return status;
}


int
gsl_sf_rect_to_polar(const double x, const double y,
                          gsl_sf_result * r, gsl_sf_result * theta)
{
  int stat_h = gsl_sf_hypot_e(x, y, r);
  if(r->val > 0.0) {
    theta->val = atan2(y, x);
    theta->err = 2.0 * GSL_DBL_EPSILON * fabs(theta->val);
    return stat_h;
  }
  else {
    DOMAIN_ERROR(theta);
  }
}


int gsl_sf_angle_restrict_symm_err_e(const double theta, gsl_sf_result * result)
{
  /* synthetic extended precision constants */
  const double P1 = 4 * 7.8539812564849853515625e-01;
  const double P2 = 4 * 3.7748947079307981766760e-08;
  const double P3 = 4 * 2.6951514290790594840552e-15;
  const double TwoPi = 2*(P1 + P2 + P3);

  const double y = GSL_SIGN(theta) * 2 * floor(fabs(theta)/TwoPi);
  double r = ((theta - y*P1) - y*P2) - y*P3;

  if(r >  M_PI) { r = (((r-2*P1)-2*P2)-2*P3); }  /* r-TwoPi */
  else if (r < -M_PI) r = (((r+2*P1)+2*P2)+2*P3); /* r+TwoPi */

  result->val = r;

  if(fabs(theta) > 0.0625/GSL_DBL_EPSILON) {
    result->val = GSL_NAN;
    result->err = GSL_NAN;
    GSL_ERROR ("error", GSL_ELOSS);
  }
  else if(fabs(theta) > 0.0625/GSL_SQRT_DBL_EPSILON) {
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val - theta);
    return GSL_SUCCESS;
  }
  else {
    double delta = fabs(result->val - theta);
    result->err = 2.0 * GSL_DBL_EPSILON * ((delta < M_PI) ? delta : M_PI);
    return GSL_SUCCESS;
  }
}


int gsl_sf_angle_restrict_pos_err_e(const double theta, gsl_sf_result * result)
{
  /* synthetic extended precision constants */
  const double P1 = 4 * 7.85398125648498535156e-01;
  const double P2 = 4 * 3.77489470793079817668e-08;
  const double P3 = 4 * 2.69515142907905952645e-15;
  const double TwoPi = 2*(P1 + P2 + P3);

  const double y = 2*floor(theta/TwoPi);

  double r = ((theta - y*P1) - y*P2) - y*P3;

  if(r > TwoPi) {r = (((r-2*P1)-2*P2)-2*P3); }  /* r-TwoPi */
  else if (r < 0) { /* may happen due to FP rounding */
    r = (((r+2*P1)+2*P2)+2*P3); /* r+TwoPi */
  }

  result->val = r;

  if(fabs(theta) > 0.0625/GSL_DBL_EPSILON) {
    result->val = GSL_NAN;
    result->err = fabs(result->val);
    GSL_ERROR ("error", GSL_ELOSS);
  }
  else if(fabs(theta) > 0.0625/GSL_SQRT_DBL_EPSILON) {
    result->err = GSL_DBL_EPSILON * fabs(result->val - theta);
    return GSL_SUCCESS;
  }
  else {
    double delta = fabs(result->val - theta);
    result->err = 2.0 * GSL_DBL_EPSILON * ((delta < M_PI) ? delta : M_PI);
    return GSL_SUCCESS;
  }
}


int gsl_sf_angle_restrict_symm_e(double * theta)
{
  gsl_sf_result r;
  int stat = gsl_sf_angle_restrict_symm_err_e(*theta, &r);
  *theta = r.val;
  return stat;
}


int gsl_sf_angle_restrict_pos_e(double * theta)
{
  gsl_sf_result r;
  int stat = gsl_sf_angle_restrict_pos_err_e(*theta, &r);
  *theta = r.val;
  return stat;
}


int gsl_sf_sin_err_e(const double x, const double dx, gsl_sf_result * result)
{
  int stat_s = gsl_sf_sin_e(x, result);
  result->err += fabs(cos(x) * dx);
  result->err += GSL_DBL_EPSILON * fabs(result->val);
  return stat_s;
}


int gsl_sf_cos_err_e(const double x, const double dx, gsl_sf_result * result)
{
  int stat_c = gsl_sf_cos_e(x, result);
  result->err += fabs(sin(x) * dx);
  result->err += GSL_DBL_EPSILON * fabs(result->val);
  return stat_c;
}


#if 0
int
gsl_sf_sin_pi_x_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(-100.0 < x && x < 100.0) {
    result->val = sin(M_PI * x) / (M_PI * x);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const double N = floor(x + 0.5);
    const double f = x - N;

    if(N < INT_MAX && N > INT_MIN) {
      /* Make it an integer if we can. Saves another
       * call to floor().
       */
      const int intN    = (int)N;
      const double sign = ( GSL_IS_ODD(intN) ? -1.0 : 1.0 );
      result->val = sign * sin(M_PI * f);
      result->err = GSL_DBL_EPSILON * fabs(result->val);
    }
    else if(N > 2.0/GSL_DBL_EPSILON || N < -2.0/GSL_DBL_EPSILON) {
      /* All integer-valued floating point numbers
       * bigger than 2/eps=2^53 are actually even.
       */
      result->val = 0.0;
      result->err = 0.0;
    }
    else {
      const double resN = N - 2.0*floor(0.5*N); /* 0 for even N, 1 for odd N */
      const double sign = ( fabs(resN) > 0.5 ? -1.0 : 1.0 );
      result->val = sign * sin(M_PI*f);
      result->err = GSL_DBL_EPSILON * fabs(result->val);
    }

    return GSL_SUCCESS;
  }
}
#endif


int gsl_sf_sinc_e(double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  {
    const double ax = fabs(x);

    if(ax < 0.8) {
      /* Do not go to the limit of the fit since
       * there is a zero there and the Chebyshev
       * accuracy will go to zero.
       */
      return cheb_eval_e(&sinc_cs, 2.0*ax-1.0, result);
    }
    else if(ax < 100.0) {
      /* Small arguments are no problem.
       * We trust the library sin() to
       * roughly machine precision.
       */
      result->val = sin(M_PI * ax)/(M_PI * ax);
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      /* Large arguments must be handled separately.
       */
      const double r = M_PI*ax;
      gsl_sf_result s;
      int stat_s = gsl_sf_sin_e(r, &s);
      result->val = s.val/r;
      result->err = s.err/r + 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return stat_s;
    }
  }
}



/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_sin(const double x)
{
  EVAL_RESULT(gsl_sf_sin_e(x, &result));
}

double gsl_sf_cos(const double x)
{
  EVAL_RESULT(gsl_sf_cos_e(x, &result));
}

double gsl_sf_hypot(const double x, const double y)
{
  EVAL_RESULT(gsl_sf_hypot_e(x, y, &result));
}

double gsl_sf_lnsinh(const double x)
{
  EVAL_RESULT(gsl_sf_lnsinh_e(x, &result));
}

double gsl_sf_lncosh(const double x)
{
  EVAL_RESULT(gsl_sf_lncosh_e(x, &result));
}

double gsl_sf_angle_restrict_symm(const double theta)
{
  double result = theta;
  EVAL_DOUBLE(gsl_sf_angle_restrict_symm_e(&result));
}

double gsl_sf_angle_restrict_pos(const double theta)
{
  double result = theta;
  EVAL_DOUBLE(gsl_sf_angle_restrict_pos_e(&result));
}

#if 0
double gsl_sf_sin_pi_x(const double x)
{
  EVAL_RESULT(gsl_sf_sin_pi_x_e(x, &result));
}
#endif

double gsl_sf_sinc(const double x)
{
  EVAL_RESULT(gsl_sf_sinc_e(x, &result));
}
} // trig.c

// zeta.c ///////////////////////////////////////////////////////////////////////
extern "C" { 
#define LogTwoPi_  1.8378770664093454835606594728111235279723


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* chebyshev fit for (s(t)-1)Zeta[s(t)]
 * s(t)= (t+1)/2
 * -1 <= t <= 1
 */
static double zeta_xlt1_data[14] = {
  1.48018677156931561235192914649,
  0.25012062539889426471999938167,
  0.00991137502135360774243761467,
 -0.00012084759656676410329833091,
 -4.7585866367662556504652535281e-06,
  2.2229946694466391855561441361e-07,
 -2.2237496498030257121309056582e-09,
 -1.0173226513229028319420799028e-10,
  4.3756643450424558284466248449e-12,
 -6.2229632593100551465504090814e-14,
 -6.6116201003272207115277520305e-16,
  4.9477279533373912324518463830e-17,
 -1.0429819093456189719660003522e-18,
  6.9925216166580021051464412040e-21,
};
static cheb_series zeta_xlt1_cs = {
  zeta_xlt1_data,
  13,
  -1, 1,
  8
};

/* chebyshev fit for (s(t)-1)Zeta[s(t)]
 * s(t)= (19t+21)/2
 * -1 <= t <= 1
 */
static double zeta_xgt1_data[30] = {
  19.3918515726724119415911269006,
   9.1525329692510756181581271500,
   0.2427897658867379985365270155,
  -0.1339000688262027338316641329,
   0.0577827064065028595578410202,
  -0.0187625983754002298566409700,
   0.0039403014258320354840823803,
  -0.0000581508273158127963598882,
  -0.0003756148907214820704594549,
   0.0001892530548109214349092999,
  -0.0000549032199695513496115090,
   8.7086484008939038610413331863e-6,
   6.4609477924811889068410083425e-7,
  -9.6749773915059089205835337136e-7,
   3.6585400766767257736982342461e-7,
  -8.4592516427275164351876072573e-8,
   9.9956786144497936572288988883e-9,
   1.4260036420951118112457144842e-9,
  -1.1761968823382879195380320948e-9,
   3.7114575899785204664648987295e-10,
  -7.4756855194210961661210215325e-11,
   7.8536934209183700456512982968e-12,
   9.9827182259685539619810406271e-13,
  -7.5276687030192221587850302453e-13,
   2.1955026393964279988917878654e-13,
  -4.1934859852834647427576319246e-14,
   4.6341149635933550715779074274e-15,
   2.3742488509048340106830309402e-16,
  -2.7276516388124786119323824391e-16,
   7.8473570134636044722154797225e-17
};
static cheb_series zeta_xgt1_cs = {
  zeta_xgt1_data,
  29,
  -1, 1,
  17
};


/* chebyshev fit for Ln[Zeta[s(t)] - 1 - 2^(-s(t))]
 * s(t)= 10 + 5t
 * -1 <= t <= 1; 5 <= s <= 15
 */
static double zetam1_inter_data[24] = {
  -21.7509435653088483422022339374,
  -5.63036877698121782876372020472,
   0.0528041358684229425504861579635,
  -0.0156381809179670789342700883562,
   0.00408218474372355881195080781927,
  -0.0010264867349474874045036628282,
   0.000260469880409886900143834962387,
  -0.0000676175847209968878098566819447,
   0.0000179284472587833525426660171124,
  -4.83238651318556188834107605116e-6,
   1.31913788964999288471371329447e-6,
  -3.63760500656329972578222188542e-7,
   1.01146847513194744989748396574e-7,
  -2.83215225141806501619105289509e-8,
   7.97733710252021423361012829496e-9,
  -2.25850168553956886676250696891e-9,
   6.42269392950164306086395744145e-10,
  -1.83363861846127284505060843614e-10,
   5.25309763895283179960368072104e-11,
  -1.50958687042589821074710575446e-11,
   4.34997545516049244697776942981e-12,
  -1.25597782748190416118082322061e-12,
   3.61280740072222650030134104162e-13,
  -9.66437239205745207188920348801e-14
}; 
static cheb_series zetam1_inter_cs = {
  zetam1_inter_data,
  22,
  -1, 1,
  12
};



/* assumes s >= 0 and s != 1.0 */
inline
static int
riemann_zeta_sgt0(double s, gsl_sf_result * result)
{
  if(s < 1.0) {
    gsl_sf_result c;
    cheb_eval_e(&zeta_xlt1_cs, 2.0*s - 1.0, &c);
    result->val = c.val / (s - 1.0);
    result->err = c.err / fabs(s-1.0) + GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(s <= 20.0) {
    double x = (2.0*s - 21.0)/19.0;
    gsl_sf_result c;
    cheb_eval_e(&zeta_xgt1_cs, x, &c);
    result->val = c.val / (s - 1.0);
    result->err = c.err / (s - 1.0) + GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    double f2 = 1.0 - pow(2.0,-s);
    double f3 = 1.0 - pow(3.0,-s);
    double f5 = 1.0 - pow(5.0,-s);
    double f7 = 1.0 - pow(7.0,-s);
    result->val = 1.0/(f2*f3*f5*f7);
    result->err = 3.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


inline
static int
riemann_zeta1ms_slt0(double s, gsl_sf_result * result)
{
  if(s > -19.0) {
    double x = (-19 - 2.0*s)/19.0;
    gsl_sf_result c;
    cheb_eval_e(&zeta_xgt1_cs, x, &c);
    result->val = c.val / (-s);
    result->err = c.err / (-s) + GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    double f2 = 1.0 - pow(2.0,-(1.0-s));
    double f3 = 1.0 - pow(3.0,-(1.0-s));
    double f5 = 1.0 - pow(5.0,-(1.0-s));
    double f7 = 1.0 - pow(7.0,-(1.0-s));
    result->val = 1.0/(f2*f3*f5*f7);
    result->err = 3.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


/* works for 5 < s < 15*/
static int
riemann_zeta_minus_1_intermediate_s(double s, gsl_sf_result * result)
{
  double t = (s - 10.0)/5.0;
  gsl_sf_result c;
  cheb_eval_e(&zetam1_inter_cs, t, &c);
  result->val = exp(c.val) + pow(2.0, -s);
  result->err = (c.err + 2.0*GSL_DBL_EPSILON)*result->val;
  return GSL_SUCCESS;
}


/* assumes s is large and positive
 * write: zeta(s) - 1 = zeta(s) * (1 - 1/zeta(s))
 * and expand a few terms of the product formula to evaluate 1 - 1/zeta(s)
 *
 * works well for s > 15
 */
static int
riemann_zeta_minus1_large_s(double s, gsl_sf_result * result)
{
  double a = pow( 2.0,-s);
  double b = pow( 3.0,-s);
  double c = pow( 5.0,-s);
  double d = pow( 7.0,-s);
  double e = pow(11.0,-s);
  double f = pow(13.0,-s);
  double t1 = a + b + c + d + e + f;
  double t2 = a*(b+c+d+e+f) + b*(c+d+e+f) + c*(d+e+f) + d*(e+f) + e*f;
  /*
  double t3 = a*(b*(c+d+e+f) + c*(d+e+f) + d*(e+f) + e*f) +
              b*(c*(d+e+f) + d*(e+f) + e*f) +
              c*(d*(e+f) + e*f) +
              d*e*f;
  double t4 = a*(b*(c*(d + e + f) + d*(e + f) + e*f) + c*(d*(e+f) + e*f) + d*e*f) +
              b*(c*(d*(e+f) + e*f) + d*e*f) +
              c*d*e*f;
  double t5 = b*c*d*e*f + a*c*d*e*f+ a*b*d*e*f+ a*b*c*e*f+ a*b*c*d*f+ a*b*c*d*e;
  double t6 = a*b*c*d*e*f;
  */
  double numt = t1 - t2 /* + t3 - t4 + t5 - t6 */;
  double zeta = 1.0/((1.0-a)*(1.0-b)*(1.0-c)*(1.0-d)*(1.0-e)*(1.0-f));
  result->val = numt*zeta;
  result->err = (15.0/s + 1.0) * 6.0*GSL_DBL_EPSILON*result->val;
  return GSL_SUCCESS;
}


#if 0
/* zeta(n) */
#define ZETA_POS_TABLE_NMAX   100
static double zeta_pos_int_table_OLD[ZETA_POS_TABLE_NMAX+1] = {
 -0.50000000000000000000000000000,       /* zeta(0) */
  0.0 /* FIXME: DirectedInfinity() */,   /* zeta(1) */
  1.64493406684822643647241516665,       /* ...     */
  1.20205690315959428539973816151,
  1.08232323371113819151600369654,
  1.03692775514336992633136548646,
  1.01734306198444913971451792979,
  1.00834927738192282683979754985,
  1.00407735619794433937868523851,
  1.00200839282608221441785276923,
  1.00099457512781808533714595890,
  1.00049418860411946455870228253,
  1.00024608655330804829863799805,
  1.00012271334757848914675183653,
  1.00006124813505870482925854511,
  1.00003058823630702049355172851,
  1.00001528225940865187173257149,
  1.00000763719763789976227360029,
  1.00000381729326499983985646164,
  1.00000190821271655393892565696,
  1.00000095396203387279611315204,
  1.00000047693298678780646311672,
  1.00000023845050272773299000365,
  1.00000011921992596531107306779,
  1.00000005960818905125947961244,
  1.00000002980350351465228018606,
  1.00000001490155482836504123466,
  1.00000000745071178983542949198,
  1.00000000372533402478845705482,
  1.00000000186265972351304900640,
  1.00000000093132743241966818287,
  1.00000000046566290650337840730,
  1.00000000023283118336765054920,
  1.00000000011641550172700519776,
  1.00000000005820772087902700889,
  1.00000000002910385044497099687,
  1.00000000001455192189104198424,
  1.00000000000727595983505748101,
  1.00000000000363797954737865119,
  1.00000000000181898965030706595,
  1.00000000000090949478402638893,
  1.00000000000045474737830421540,
  1.00000000000022737368458246525,
  1.00000000000011368684076802278,
  1.00000000000005684341987627586,
  1.00000000000002842170976889302,
  1.00000000000001421085482803161,
  1.00000000000000710542739521085,
  1.00000000000000355271369133711,
  1.00000000000000177635684357912,
  1.00000000000000088817842109308,
  1.00000000000000044408921031438,
  1.00000000000000022204460507980,
  1.00000000000000011102230251411,
  1.00000000000000005551115124845,
  1.00000000000000002775557562136,
  1.00000000000000001387778780973,
  1.00000000000000000693889390454,
  1.00000000000000000346944695217,
  1.00000000000000000173472347605,
  1.00000000000000000086736173801,
  1.00000000000000000043368086900,
  1.00000000000000000021684043450,
  1.00000000000000000010842021725,
  1.00000000000000000005421010862,
  1.00000000000000000002710505431,
  1.00000000000000000001355252716,
  1.00000000000000000000677626358,
  1.00000000000000000000338813179,
  1.00000000000000000000169406589,
  1.00000000000000000000084703295,
  1.00000000000000000000042351647,
  1.00000000000000000000021175824,
  1.00000000000000000000010587912,
  1.00000000000000000000005293956,
  1.00000000000000000000002646978,
  1.00000000000000000000001323489,
  1.00000000000000000000000661744,
  1.00000000000000000000000330872,
  1.00000000000000000000000165436,
  1.00000000000000000000000082718,
  1.00000000000000000000000041359,
  1.00000000000000000000000020680,
  1.00000000000000000000000010340,
  1.00000000000000000000000005170,
  1.00000000000000000000000002585,
  1.00000000000000000000000001292,
  1.00000000000000000000000000646,
  1.00000000000000000000000000323,
  1.00000000000000000000000000162,
  1.00000000000000000000000000081,
  1.00000000000000000000000000040,
  1.00000000000000000000000000020,
  1.00000000000000000000000000010,
  1.00000000000000000000000000005,
  1.00000000000000000000000000003,
  1.00000000000000000000000000001,
  1.00000000000000000000000000001,
  1.00000000000000000000000000000,
  1.00000000000000000000000000000,
  1.00000000000000000000000000000
};
#endif /* 0 */


/* zeta(n) - 1 */
#define ZETA_POS_TABLE_NMAX   100
static double zetam1_pos_int_table[ZETA_POS_TABLE_NMAX+1] = {
 -1.5,                               /* zeta(0) */
  0.0,       /* FIXME: Infinity */   /* zeta(1) - 1 */
  0.644934066848226436472415166646,  /* zeta(2) - 1 */
  0.202056903159594285399738161511,
  0.082323233711138191516003696541,
  0.036927755143369926331365486457,
  0.017343061984449139714517929790,
  0.008349277381922826839797549849,
  0.004077356197944339378685238508,
  0.002008392826082214417852769232,
  0.000994575127818085337145958900,
  0.000494188604119464558702282526,
  0.000246086553308048298637998047,
  0.000122713347578489146751836526,
  0.000061248135058704829258545105,
  0.000030588236307020493551728510,
  0.000015282259408651871732571487,
  7.6371976378997622736002935630e-6,
  3.8172932649998398564616446219e-6,
  1.9082127165539389256569577951e-6,
  9.5396203387279611315203868344e-7,
  4.7693298678780646311671960437e-7,
  2.3845050272773299000364818675e-7,
  1.1921992596531107306778871888e-7,
  5.9608189051259479612440207935e-8,
  2.9803503514652280186063705069e-8,
  1.4901554828365041234658506630e-8,
  7.4507117898354294919810041706e-9,
  3.7253340247884570548192040184e-9,
  1.8626597235130490064039099454e-9,
  9.3132743241966818287176473502e-10,
  4.6566290650337840729892332512e-10,
  2.3283118336765054920014559759e-10,
  1.1641550172700519775929738354e-10,
  5.8207720879027008892436859891e-11,
  2.9103850444970996869294252278e-11,
  1.4551921891041984235929632245e-11,
  7.2759598350574810145208690123e-12,
  3.6379795473786511902372363558e-12,
  1.8189896503070659475848321007e-12,
  9.0949478402638892825331183869e-13,
  4.5474737830421540267991120294e-13,
  2.2737368458246525152268215779e-13,
  1.1368684076802278493491048380e-13,
  5.6843419876275856092771829675e-14,
  2.8421709768893018554550737049e-14,
  1.4210854828031606769834307141e-14,
  7.1054273952108527128773544799e-15,
  3.5527136913371136732984695340e-15,
  1.7763568435791203274733490144e-15,
  8.8817842109308159030960913863e-16,
  4.4408921031438133641977709402e-16,
  2.2204460507980419839993200942e-16,
  1.1102230251410661337205445699e-16,
  5.5511151248454812437237365905e-17,
  2.7755575621361241725816324538e-17,
  1.3877787809725232762839094906e-17,
  6.9388939045441536974460853262e-18,
  3.4694469521659226247442714961e-18,
  1.7347234760475765720489729699e-18,
  8.6736173801199337283420550673e-19,
  4.3368086900206504874970235659e-19,
  2.1684043449972197850139101683e-19,
  1.0842021724942414063012711165e-19,
  5.4210108624566454109187004043e-20,
  2.7105054312234688319546213119e-20,
  1.3552527156101164581485233996e-20,
  6.7762635780451890979952987415e-21,
  3.3881317890207968180857031004e-21,
  1.6940658945097991654064927471e-21,
  8.4703294725469983482469926091e-22,
  4.2351647362728333478622704833e-22,
  2.1175823681361947318442094398e-22,
  1.0587911840680233852265001539e-22,
  5.2939559203398703238139123029e-23,
  2.6469779601698529611341166842e-23,
  1.3234889800848990803094510250e-23,
  6.6174449004244040673552453323e-24,
  3.3087224502121715889469563843e-24,
  1.6543612251060756462299236771e-24,
  8.2718061255303444036711056167e-25,
  4.1359030627651609260093824555e-25,
  2.0679515313825767043959679193e-25,
  1.0339757656912870993284095591e-25,
  5.1698788284564313204101332166e-26,
  2.5849394142282142681277617708e-26,
  1.2924697071141066700381126118e-26,
  6.4623485355705318034380021611e-27,
  3.2311742677852653861348141180e-27,
  1.6155871338926325212060114057e-27,
  8.0779356694631620331587381863e-28,
  4.0389678347315808256222628129e-28,
  2.0194839173657903491587626465e-28,
  1.0097419586828951533619250700e-28,
  5.0487097934144756960847711725e-29,
  2.5243548967072378244674341938e-29,
  1.2621774483536189043753999660e-29,
  6.3108872417680944956826093943e-30,
  3.1554436208840472391098412184e-30,
  1.5777218104420236166444327830e-30,
  7.8886090522101180735205378276e-31
};


#define ZETA_NEG_TABLE_NMAX  99
#define ZETA_NEG_TABLE_SIZE  50
static double zeta_neg_int_table[ZETA_NEG_TABLE_SIZE] = {
 -0.083333333333333333333333333333,     /* zeta(-1) */
  0.008333333333333333333333333333,     /* zeta(-3) */
 -0.003968253968253968253968253968,     /* ...      */
  0.004166666666666666666666666667,
 -0.007575757575757575757575757576,
  0.021092796092796092796092796093,
 -0.083333333333333333333333333333,
  0.44325980392156862745098039216,
 -3.05395433027011974380395433027,
  26.4562121212121212121212121212,
 -281.460144927536231884057971014,
  3607.5105463980463980463980464,
 -54827.583333333333333333333333,
  974936.82385057471264367816092,
 -2.0052695796688078946143462272e+07,
  4.7238486772162990196078431373e+08,
 -1.2635724795916666666666666667e+10,
  3.8087931125245368811553022079e+11,
 -1.2850850499305083333333333333e+13,
  4.8241448354850170371581670362e+14,
 -2.0040310656516252738108421663e+16,
  9.1677436031953307756992753623e+17,
 -4.5979888343656503490437943262e+19,
  2.5180471921451095697089023320e+21,
 -1.5001733492153928733711440151e+23,
  9.6899578874635940656497942895e+24,
 -6.7645882379292820990945242302e+26,
  5.0890659468662289689766332916e+28,
 -4.1147288792557978697665486068e+30,
  3.5666582095375556109684574609e+32,
 -3.3066089876577576725680214670e+34,
  3.2715634236478716264211227016e+36,
 -3.4473782558278053878256455080e+38,
  3.8614279832705258893092720200e+40,
 -4.5892974432454332168863989006e+42,
  5.7775386342770431824884825688e+44,
 -7.6919858759507135167410075972e+46,
  1.0813635449971654696354033351e+49,
 -1.6029364522008965406067102346e+51,
  2.5019479041560462843656661499e+53,
 -4.1067052335810212479752045004e+55,
  7.0798774408494580617452972433e+57,
 -1.2804546887939508790190849756e+60,
  2.4267340392333524078020892067e+62,
 -4.8143218874045769355129570066e+64,
  9.9875574175727530680652777408e+66,
 -2.1645634868435185631335136160e+69,
  4.8962327039620553206849224516e+71,    /* ...        */
 -1.1549023923963519663954271692e+74,    /* zeta(-97)  */
  2.8382249570693706959264156336e+76     /* zeta(-99)  */
};


/* coefficients for Maclaurin summation in hzeta()
 * B_{2j}/(2j)!
 */
static double hzeta_c[15] = {
  1.00000000000000000000000000000,
  0.083333333333333333333333333333,
 -0.00138888888888888888888888888889,
  0.000033068783068783068783068783069,
 -8.2671957671957671957671957672e-07,
  2.0876756987868098979210090321e-08,
 -5.2841901386874931848476822022e-10,
  1.3382536530684678832826980975e-11,
 -3.3896802963225828668301953912e-13,
  8.5860620562778445641359054504e-15,
 -2.1748686985580618730415164239e-16,
  5.5090028283602295152026526089e-18,
 -1.3954464685812523340707686264e-19,
  3.5347070396294674716932299778e-21,
 -8.9535174270375468504026113181e-23
};

#define ETA_POS_TABLE_NMAX  100
static double eta_pos_int_table[ETA_POS_TABLE_NMAX+1] = {
0.50000000000000000000000000000,  /* eta(0) */
M_LN2,                            /* eta(1) */
0.82246703342411321823620758332,  /* ...    */
0.90154267736969571404980362113,
0.94703282949724591757650323447,
0.97211977044690930593565514355,
0.98555109129743510409843924448,
0.99259381992283028267042571313,
0.99623300185264789922728926008,
0.99809429754160533076778303185,
0.99903950759827156563922184570,
0.99951714349806075414409417483,
0.99975768514385819085317967871,
0.99987854276326511549217499282,
0.99993917034597971817095419226,
0.99996955121309923808263293263,
0.99998476421490610644168277496,
0.99999237829204101197693787224,
0.99999618786961011347968922641,
0.99999809350817167510685649297,
0.99999904661158152211505084256,
0.99999952325821554281631666433,
0.99999976161323082254789720494,
0.99999988080131843950322382485,
0.99999994039889239462836140314,
0.99999997019885696283441513311,
0.99999998509923199656878766181,
0.99999999254955048496351585274,
0.99999999627475340010872752767,
0.99999999813736941811218674656,
0.99999999906868228145397862728,
0.99999999953434033145421751469,
0.99999999976716989595149082282,
0.99999999988358485804603047265,
0.99999999994179239904531592388,
0.99999999997089618952980952258,
0.99999999998544809143388476396,
0.99999999999272404460658475006,
0.99999999999636202193316875550,
0.99999999999818101084320873555,
0.99999999999909050538047887809,
0.99999999999954525267653087357,
0.99999999999977262633369589773,
0.99999999999988631316532476488,
0.99999999999994315658215465336,
0.99999999999997157829090808339,
0.99999999999998578914539762720,
0.99999999999999289457268000875,
0.99999999999999644728633373609,
0.99999999999999822364316477861,
0.99999999999999911182158169283,
0.99999999999999955591079061426,
0.99999999999999977795539522974,
0.99999999999999988897769758908,
0.99999999999999994448884878594,
0.99999999999999997224442439010,
0.99999999999999998612221219410,
0.99999999999999999306110609673,
0.99999999999999999653055304826,
0.99999999999999999826527652409,
0.99999999999999999913263826204,
0.99999999999999999956631913101,
0.99999999999999999978315956551,
0.99999999999999999989157978275,
0.99999999999999999994578989138,
0.99999999999999999997289494569,
0.99999999999999999998644747284,
0.99999999999999999999322373642,
0.99999999999999999999661186821,
0.99999999999999999999830593411,
0.99999999999999999999915296705,
0.99999999999999999999957648353,
0.99999999999999999999978824176,
0.99999999999999999999989412088,
0.99999999999999999999994706044,
0.99999999999999999999997353022,
0.99999999999999999999998676511,
0.99999999999999999999999338256,
0.99999999999999999999999669128,
0.99999999999999999999999834564,
0.99999999999999999999999917282,
0.99999999999999999999999958641,
0.99999999999999999999999979320,
0.99999999999999999999999989660,
0.99999999999999999999999994830,
0.99999999999999999999999997415,
0.99999999999999999999999998708,
0.99999999999999999999999999354,
0.99999999999999999999999999677,
0.99999999999999999999999999838,
0.99999999999999999999999999919,
0.99999999999999999999999999960,
0.99999999999999999999999999980,
0.99999999999999999999999999990,
0.99999999999999999999999999995,
0.99999999999999999999999999997,
0.99999999999999999999999999999,
0.99999999999999999999999999999,
1.00000000000000000000000000000,
1.00000000000000000000000000000,
1.00000000000000000000000000000,
};


#define ETA_NEG_TABLE_NMAX  99
#define ETA_NEG_TABLE_SIZE  50
static double eta_neg_int_table[ETA_NEG_TABLE_SIZE] = {
 0.25000000000000000000000000000,   /* eta(-1) */
-0.12500000000000000000000000000,   /* eta(-3) */
 0.25000000000000000000000000000,   /* ...      */
-1.06250000000000000000000000000,
 7.75000000000000000000000000000,
-86.3750000000000000000000000000,
 1365.25000000000000000000000000,
-29049.0312500000000000000000000,
 800572.750000000000000000000000,
-2.7741322625000000000000000000e+7,
 1.1805291302500000000000000000e+9,
-6.0523980051687500000000000000e+10,
 3.6794167785377500000000000000e+12,
-2.6170760990658387500000000000e+14,
 2.1531418140800295250000000000e+16,
-2.0288775575173015930156250000e+18,
 2.1708009902623770590275000000e+20,
-2.6173826968455814932120125000e+22,
 3.5324148876863877826668602500e+24,
-5.3042033406864906641493838981e+26,
 8.8138218364311576767253114668e+28,
-1.6128065107490778547354654864e+31,
 3.2355470001722734208527794569e+33,
-7.0876727476537493198506645215e+35,
 1.6890450341293965779175629389e+38,
-4.3639690731216831157655651358e+40,
 1.2185998827061261322605065672e+43,
-3.6670584803153006180101262324e+45,
 1.1859898526302099104271449748e+48,
-4.1120769493584015047981746438e+50,
 1.5249042436787620309090168687e+53,
-6.0349693196941307074572991901e+55,
 2.5437161764210695823197691519e+58,
-1.1396923802632287851130360170e+61,
 5.4180861064753979196802726455e+63,
-2.7283654799994373847287197104e+66,
 1.4529750514918543238511171663e+69,
-8.1705519371067450079777183386e+71,
 4.8445781606678367790247757259e+74,
-3.0246694206649519336179448018e+77,
 1.9858807961690493054169047970e+80,
-1.3694474620720086994386818232e+83,
 9.9070382984295807826303785989e+85,
-7.5103780796592645925968460677e+88,
 5.9598418264260880840077992227e+91,
-4.9455988887500020399263196307e+94,
 4.2873596927020241277675775935e+97,
-3.8791952037716162900707994047e+100,
 3.6600317773156342245401829308e+103,
-3.5978775704117283875784869570e+106    /* eta(-99)  */
};


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/


int gsl_sf_hzeta_e(const double s, const double q, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(s <= 1.0 || q <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else {
    const double max_bits = 54.0;
    const double ln_term0 = -s * log(q);  

    if(ln_term0 < GSL_LOG_DBL_MIN + 1.0) {
      UNDERFLOW_ERROR(result);
    }
    else if(ln_term0 > GSL_LOG_DBL_MAX - 1.0) {
      OVERFLOW_ERROR (result);
    }
    else if((s > max_bits && q < 1.0) || (s > 0.5*max_bits && q < 0.25)) {
      result->val = pow(q, -s);
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if(s > 0.5*max_bits && q < 1.0) {
      const double p1 = pow(q, -s);
      const double p2 = pow(q/(1.0+q), s);
      const double p3 = pow(q/(2.0+q), s);
      result->val = p1 * (1.0 + p2 + p3);
      result->err = GSL_DBL_EPSILON * (0.5*s + 2.0) * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      /* Euler-Maclaurin summation formula 
       * [Moshier, p. 400, with several typo corrections]
       */
      const int jmax = 12;
      const int kmax = 10;
      int j, k;
      const double pmax  = pow(kmax + q, -s);
      double scp = s;
      double pcp = pmax / (kmax + q);
      double ans = pmax*((kmax+q)/(s-1.0) + 0.5);

      for(k=0; k<kmax; k++) {
        ans += pow(k + q, -s);
      }

      for(j=0; j<=jmax; j++) {
        double delta = hzeta_c[j+1] * scp * pcp;
        ans += delta;
        if(fabs(delta/ans) < 0.5*GSL_DBL_EPSILON) break;
        scp *= (s+2*j+1)*(s+2*j+2);
        pcp /= (kmax + q)*(kmax + q);
      }

      result->val = ans;
      result->err = 2.0 * (jmax + 1.0) * GSL_DBL_EPSILON * fabs(ans);
      return GSL_SUCCESS;
    }
  }
}


int gsl_sf_zeta_e(const double s, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(s == 1.0) {
    DOMAIN_ERROR(result);
  }
  else if(s >= 0.0) {
    return riemann_zeta_sgt0(s, result);
  }
  else {
    /* reflection formula, [Abramowitz+Stegun, 23.2.5] */

    gsl_sf_result zeta_one_minus_s;
    const int stat_zoms = riemann_zeta1ms_slt0(s, &zeta_one_minus_s);
    const double sin_term = (fmod(s,2.0) == 0.0) ? 0.0 : sin(0.5*M_PI*fmod(s,4.0))/M_PI;

    if(sin_term == 0.0) {
      result->val = 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
    else if(s > -170) {
      /* We have to be careful about losing digits
       * in calculating pow(2 Pi, s). The gamma
       * function is fine because we were careful
       * with that implementation.
       * We keep an array of (2 Pi)^(10 n).
       */
      const double twopi_pow[18] = { 1.0,
                                     9.589560061550901348e+007,
                                     9.195966217409212684e+015,
                                     8.818527036583869903e+023,
                                     8.456579467173150313e+031,
                                     8.109487671573504384e+039,
                                     7.776641909496069036e+047,
                                     7.457457466828644277e+055,
                                     7.151373628461452286e+063,
                                     6.857852693272229709e+071,
                                     6.576379029540265771e+079,
                                     6.306458169130020789e+087,
                                     6.047615938853066678e+095,
                                     5.799397627482402614e+103,
                                     5.561367186955830005e+111,
                                     5.333106466365131227e+119,
                                     5.114214477385391780e+127,
                                     4.904306689854036836e+135
                                    };
      const int n = floor((-s)/10.0);
      const double fs = s + 10.0*n;
      const double p = pow(2.0*M_PI, fs) / twopi_pow[n];

      gsl_sf_result g;
      const int stat_g = gsl_sf_gamma_e(1.0-s, &g);
      result->val  = p * g.val * sin_term * zeta_one_minus_s.val;
      result->err  = fabs(p * g.val * sin_term) * zeta_one_minus_s.err;
      result->err += fabs(p * sin_term * zeta_one_minus_s.val) * g.err;
      result->err += GSL_DBL_EPSILON * (fabs(s)+2.0) * fabs(result->val);
      return GSL_ERROR_SELECT_2(stat_g, stat_zoms);
    }
    else {
      /* The actual zeta function may or may not
       * overflow here. But we have no easy way
       * to calculate it when the prefactor(s)
       * overflow. Trying to use log's and exp
       * is no good because we loose a couple
       * digits to the exp error amplification.
       * When we gather a little more patience,
       * we can implement something here. Until
       * then just give up.
       */
      OVERFLOW_ERROR(result);
    }
  }
}


int gsl_sf_zeta_int_e(const int n, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(n < 0) {
    if(!GSL_IS_ODD(n)) {
      result->val = 0.0; /* exactly zero at even negative integers */
      result->err = 0.0;
      return GSL_SUCCESS;
    }
    else if(n > -ZETA_NEG_TABLE_NMAX) {
      result->val = zeta_neg_int_table[-(n+1)/2];
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      return gsl_sf_zeta_e((double)n, result);
    }
  }
  else if(n == 1){
    DOMAIN_ERROR(result);
  }
  else if(n <= ZETA_POS_TABLE_NMAX){
    result->val = 1.0 + zetam1_pos_int_table[n];
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = 1.0;
    result->err = GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
}


int gsl_sf_zetam1_e(const double s, gsl_sf_result * result)
{
  if(s <= 5.0)
  {
    int stat = gsl_sf_zeta_e(s, result);
    result->val = result->val - 1.0;
    return stat;
  }
  else if(s < 15.0)
  {
    return riemann_zeta_minus_1_intermediate_s(s, result);
  }
  else
  {
    return riemann_zeta_minus1_large_s(s, result);
  }
}


int gsl_sf_zetam1_int_e(const int n, gsl_sf_result * result)
{
  if(n < 0) {
    if(!GSL_IS_ODD(n)) {
      result->val = -1.0; /* at even negative integers zetam1 == -1 since zeta is exactly zero */
      result->err = 0.0;
      return GSL_SUCCESS;
    }
    else if(n > -ZETA_NEG_TABLE_NMAX) {
      result->val = zeta_neg_int_table[-(n+1)/2] - 1.0;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      /* could use gsl_sf_zetam1_e here but subtracting 1 makes no difference
         for such large values, so go straight to the result */
      return gsl_sf_zeta_e((double)n, result);  
    }
  }
  else if(n == 1){
    DOMAIN_ERROR(result);
  }
  else if(n <= ZETA_POS_TABLE_NMAX){
    result->val = zetam1_pos_int_table[n];
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    return gsl_sf_zetam1_e(n, result);
  }
}


int gsl_sf_eta_int_e(int n, gsl_sf_result * result)
{
  if(n > ETA_POS_TABLE_NMAX) {
    result->val = 1.0;
    result->err = GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
  else if(n >= 0) {
    result->val = eta_pos_int_table[n];
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    /* n < 0 */

    if(!GSL_IS_ODD(n)) {
      /* exactly zero at even negative integers */
      result->val = 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
    else if(n > -ETA_NEG_TABLE_NMAX) {
      result->val = eta_neg_int_table[-(n+1)/2];
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      gsl_sf_result z;
      gsl_sf_result p;
      int stat_z = gsl_sf_zeta_int_e(n, &z);
      int stat_p = gsl_sf_exp_e((1.0-n)*M_LN2, &p);
      int stat_m = gsl_sf_multiply_e(-p.val, z.val, result);
      result->err  = fabs(p.err * (M_LN2*(1.0-n)) * z.val) + z.err * fabs(p.val);
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_ERROR_SELECT_3(stat_m, stat_p, stat_z);
    }
  }
}


int gsl_sf_eta_e(const double s, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(s > 100.0) {
    result->val = 1.0;
    result->err = GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
  else if(fabs(s-1.0) < 10.0*GSL_ROOT5_DBL_EPSILON) {
    double del = s-1.0;
    double c0  = M_LN2;
    double c1  = M_LN2 * (M_EULER - 0.5*M_LN2);
    double c2  = -0.0326862962794492996;
    double c3  =  0.0015689917054155150;
    double c4  =  0.00074987242112047532;
    result->val = c0 + del * (c1 + del * (c2 + del * (c3 + del * c4)));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result z;
    gsl_sf_result p;
    int stat_z = gsl_sf_zeta_e(s, &z);
    int stat_p = gsl_sf_exp_e((1.0-s)*M_LN2, &p);
    int stat_m = gsl_sf_multiply_e(1.0-p.val, z.val, result);
    result->err  = fabs(p.err * (M_LN2*(1.0-s)) * z.val) + z.err * fabs(p.val);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_ERROR_SELECT_3(stat_m, stat_p, stat_z);
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_zeta(const double s)
{
  EVAL_RESULT(gsl_sf_zeta_e(s, &result));
}

double gsl_sf_hzeta(const double s, const double a)
{
  EVAL_RESULT(gsl_sf_hzeta_e(s, a, &result));
}

double gsl_sf_zeta_int(const int s)
{
  EVAL_RESULT(gsl_sf_zeta_int_e(s, &result));
}

double gsl_sf_zetam1(const double s)
{
  EVAL_RESULT(gsl_sf_zetam1_e(s, &result));
}

double gsl_sf_zetam1_int(const int s)
{
  EVAL_RESULT(gsl_sf_zetam1_int_e(s, &result));
}

double gsl_sf_eta_int(const int s)
{
  EVAL_RESULT(gsl_sf_eta_int_e(s, &result));
}

double gsl_sf_eta(const double s)
{
  EVAL_RESULT(gsl_sf_eta_e(s, &result));
}
} // zeta.c
