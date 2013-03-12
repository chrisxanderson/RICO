#include "MyVector.h"
#include "UnitTest++/UnitTest++.h"

////////////////////////////////////////////////////////////////////////

TEST(TestMyVector) {
  MyVector<double> x(4,7,2);
  MyVector<bool>   m(1,0,1);

  //cout << "(4, 7, 2) == " << x << endl;

  const float xm_true[] = {4,2};
  CHECK_ARRAY_EQUAL(xm_true, x[m].get(), 2);

  MyVector<double> y(x);
  const float y_true[] = {4,7,2};
  CHECK_ARRAY_EQUAL(y_true, y.get(), 3);

  x += y;
  const float x_true[] = {4,7,2,4,7,2};
  CHECK_ARRAY_EQUAL(x_true, x.get(), 6);

  m += m;
  const float xm2_true[] = {4,2,4,2};
  CHECK_ARRAY_EQUAL(xm2_true, x[m].get(), 4);
}

////////////////////////////////////////////////////////////////////////
