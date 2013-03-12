#ifndef __UTIL_H__
#define __UTIL_H__

#include <string>
#include <list>
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////

class Util {

 public:
  static vector<double> Slice(vector<double> x, int start, int end);
  static void print(vector<double> x);
  static void print(vector<int> x);

  static vector<bool> False(int n);

  static vector<double> divide(vector<double> x, vector<double> y);
  static vector<double> diff(vector<double> x);
  static string join(const list<string>& L, const char* SEP);
};

////////////////////////////////////////////////////////////////////////////////////////

namespace rico {
  namespace recipe {
    double lgamma(double x);
  }
}

#endif // __UTIL_H__
