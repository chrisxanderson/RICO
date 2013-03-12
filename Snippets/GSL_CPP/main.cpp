#include "gsl.h"
#include "stdio.h"

void test_gamma() {
  int N = 14;
  double a[]  = {1e-100, 0.001, 0.001,  0.001,   1.0,  1.0,  1.0,  10.0, 10.0, 1000.0, 1000.0,\
                 34.0, 37.0, 10.0};
  double x[]  = {0.001,  0.001,   1.0,   10.0, 0.001, 1.01, 10.0, 10.01, 20.0, 1000.1, 2000.0,\
                 32.0, 3.499999999999999289e+01, 1e-16};
  double y[]  = {1.0,                            0.9936876467088602902,\
                 0.9997803916424144436,          0.9999999958306921828,\
                 0.0009995001666250083319,       0.6357810204284766802, \
                 0.9999546000702375151,          0.5433207586693410570, \
                 0.9950045876916924128,          0.5054666401440661753, \
                 1.0,                            0.3849626436463866776322932129,\
                 0.3898035054195570860969333039, 2.755731922398588814734648067e-167};
  int i;
  for(i = 0; i < N; i++) {
    double ig = gsl_sf_gamma_inc_P (a[i], x[i]);
    printf("Incomplete Gamma[a=%f](%f) = %.15f, error =  %e\n", a[i],x[i],ig,(ig-y[i]));
  }
}

void test_beta() {
  int N = 19;
  double a[]  = {1.0, 1.0, 0.1, 1.0, 0.1,\
		 10.0, 50.0, 1.0, 1.0, 1.0,\
                 1.0, 1.0, 1.0, 50.0, 90.0,\
                 500.0, 5000.0, 5000.0, 5000.0};
  double b[]  = {1.0, 1.0, 0.1, 1.0, 1.0,\
		 1.0,  1.0, 0.1, 10.0, 50.0,\
                 1.0, 2.0,  2.0, 60.0, 90.0,\
                 500.0, 5000.0, 5000.0, 2000.0};
  double x[]  = {0.0, 1.0, 1.0, 0.5, 0.5,\
		 0.5, 0.5, 0.5, 0.5, 0.5,\
                 0.1, 0.1, 0.9, 0.5, 0.5,\
                 0.6, 0.4, 0.6, 0.6};
  double y[]  = {0.0, 1.0, 1.0, 0.5, 0.9330329915368074160,\
		 0.0009765625000000000000, 8.881784197001252323e-16,\
                 0.06696700846319258402, 0.99902343750000000000,\
                 0.99999999999999911180, 0.10, 0.19, 0.99, \
		 0.8309072939016694143, 0.5, 0.9999999999157549630, \
                 4.518543727260666383e-91, 1.0, 8.445388773903332659e-89};
  int i;
  for(i = 0; i < N; i++) {
    double ib = gsl_sf_beta_inc (a[i], b[i], x[i]);
    printf("Incomplete Beta[a=%f. b=%f](%f) = %.15f, error =  %e\n", a[i],b[i],x[i],ib,(ib-y[i]));
  }
}

void test_chi_squared(double k) {
  int N = 100;
  double L = 0;
  double R = 8;
  for(int i = 0; i < N; i++) {
    double x = ((double)i)/N*(R-L) + L;
    double c = gsl_sf_gamma_inc_P (k/2, x/2);
    printf("%f\t%f\n", x, c);
  }
}

double get_student_t(double df, double t) {
   if (t == 0) return 0.5;
   //
   // Calculate probability of Student's t using the incomplete beta function.
   // probability = ibeta(df / 2, 1/2, df / (df + t*t))
   //
   // However when t is small compared to the degrees of freedom, that formula
   // suffers from rounding error, use the identity formula to work around
   // the problem:
   //
   // I[x](a,b) = 1 - I[1-x](b,a)
   //
   // and:
   //
   //     x = df / (df + t^2)
   //
   // so:
   //
   // 1 - x = t^2 / (df + t^2)
   //
   double t2 = t * t;
   double probability;
   if(df > 2 * t2 && false)  // This doesn't work
   {
      double z = t2 / (df + t2);
      probability = gsl_sf_beta_inc(0.5, df / 2, z) / 2;  // TRF: want incomplete beta complement!
   } else {
      double z = df / (df + t2);
      probability = gsl_sf_beta_inc(df / 2, 0.5, z) / 2;
   }
   return (t > 0 ? 1 - probability : probability);
}

int main() {

  //test_gamma();
  //test_beta();

  //test_chi_squared(6);
  //printf("ChiSq[6] = %f\n", gsl_sf_gamma_inc_P(6/2,14.45/2));

  printf("Student-t: %f\n",get_student_t(4, 0.741));

  return 0;
}
