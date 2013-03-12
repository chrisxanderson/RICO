#include <list>
#include <vector>
#include <cmath>

#include <iostream>
#include <iterator>
#include <string>
#include <sstream>

#include "Util.h"

#include <stdio.h>

#include "UnitTest++/UnitTest++.h"

////////////////////////////////////////////////////////////////////////////////

vector<double> Util::Slice(vector<double> x, int start, int end) {
  vector<double> y;
  while(end < 0) end = x.size() + end;
  for(int i = start; i < end; i++) y.push_back(x[i]);
  return y;
}

void Util::print(vector<double> x) {
  int N = x.size();
  printf("(");
  for(int i = 0; i < N-1; i++) printf("%f, ", x[i]);
  if(N > 0) printf("%f", x[N-1]);
  printf(")\n");
}

void Util::print(vector<int> x) {
  int N = x.size();
  printf("(");
  for(int i = 0; i < N-1; i++) printf("%d, ", x[i]);
  if(N > 0) printf("%d", x[N-1]);
  printf(")\n");
}

vector<bool> Util::False(int n) {
  vector<bool> v;
  for(int i = 0; i < n; i++) v.push_back(0);
  return v;
}


vector<double> Util::divide(vector<double> x, vector<double> y) {
  vector<double> z;
  for(int i = 0; i < x.size(); i++) z.push_back(x[i] / y[i]);
  return z;
}

vector<double> Util::diff(vector<double> x) {
  vector<double> d;
  for(int i = 0; i < x.size() - 1; i++) d.push_back(x[i+1] - x[i]);
  return d;
}

string Util::join(const list<string>& L, const char* SEP) {
  stringstream sstr;

  ostream_iterator<string> it (sstr, SEP);
  copy(L.begin(), --L.end(), it);
  sstr << L.back();

  return sstr.str();
}

TEST(TestStringJoin) {
  list<string> x;

  x.push_back("fee");
  x.push_back("fie");
  x.push_back("foe");
  x.push_back("fum");

  // cout << Util::join(x, ", ") << endl; // DEBUG

  CHECK_EQUAL(Util::join(x, ", "), "fee, fie, foe, fum");
}

//////////////////////////////////////////////////////////////////

namespace rico {
  namespace recipe {

    // log Gamma(x), x > 0
    double lgamma(double x) {
      static double ser    = 1.000000000190015;
      static double cof[6] = {76.18009172947146, -86.50532032941677,
			      24.01409824083091, -1.231739572450155,
			      0.1208650973955169e-2, -0.5395239384953e-5};
      double y    = x;
      double tmp  = x + 5.5;
      tmp        -= (x + 0.5) * log(tmp);
      
      for(int j = 0; j < 5; j++) ser += cof[j]/++y;

      return -tmp+log(2.5066282746310005 * ser / x);
    }

  } // namespace recipe
} // namespace rico

