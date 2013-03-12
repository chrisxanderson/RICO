#include <R.h>
#include "X.h"

// static Y y;    // This works, but it's just a nuisance now that we have a "real" example.

X:: X()  { REprintf("TRF: constructor X\n"); }
X::~X()  { REprintf("TRF: destructor X\n"); }
Y:: Y()  { REprintf("TRF: constructor Y\n"); }
Y::~Y()  { REprintf("TRF: destructor Y\n"); }
