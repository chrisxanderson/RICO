/*
g++ inherit.cpp -o inherit

  What we want is a way to cast B to A, but remember the symbolic name of the class and still
  recover the correct "update" when asked. 
  How about we create a vector of each virtual function allocated by enum size
  What's the basic trick? We need a C-function that does the casting
  We declare: B b; but we really have ba = (A)b;
  Wish that ba.update() would just call the B-version, but it calls the A-version
  I'm not really seeing how I can get around the swtiched version of each virtual function.
 */

#include <stdio.h>

/////////////// A ///////////////////

class A {
protected:
  int m_value;
  int m_type;
  
public:
  enum SYMBOLIC_CLASS {A_CLASS, B_CLASS, C_CLASS};

  A(int type) : m_type(type), m_value(10) {}

  virtual void update();

  int value() const {return m_value;}
  int type()  const {return m_type;}
};

/////////////// B ///////////////////

class B: public A {
public:
  B() : A(B_CLASS) {}
  void update() {printf("B::update\n"); m_value++;}
};

/////////////// C ///////////////////

class C: public A {
public:
  C() : A(C_CLASS) {}
  void update() {printf("C::update\n"); m_value--;}
};

/////////////// interface ////////////////

// This works. It lives in C space
void virtual_update(A* a) {
  switch(a->type()) {
  case A::A_CLASS: return a->update();
  case A::B_CLASS: return ((B*)a)->update();
  case A::C_CLASS: return ((C*)a)->update();
  }
}

// This works and lives in place of the pure virtual base function that is a hole anyway!
void A::update() {
  switch(type()) {
  case A::B_CLASS: return ((B*)this)->update();
  case A::C_CLASS: return ((C*)this)->update();
  default:
    printf("A::update failed to find appropraite member function.\n");
  }
}


/////////////// main ////////////////

int main() {
  B b;
  C c;

  b.update();
  c.update();
  printf("b = %d, b.type = %d\n", b.value(), b.type());  // b = 11, type = 1
  printf("c = %d, c.type = %d\n", c.value(), c.type());  // c = 9,  type = 2

  A* pb = (A*)(&b);
  A* pc = (A*)(&c);
  printf("pb = %d, type = %d\n", pb->value(), pb->type());  // b = 11, type = 1
  printf("pc = %d, type = %d\n", pc->value(), pc->type());  // c = 9,  type = 2

  virtual_update(pb);
  virtual_update(pc);
  printf("pb = %d, type = %d\n", pb->value(), pb->type());  // b = 12, type = 1
  printf("pc = %d, type = %d\n", pc->value(), pc->type());  // c = 8,  type = 2

  pb->update();
  pc->update();
  printf("pb = %d, type = %d\n", pb->value(), pb->type());  // b = 13, type = 1
  printf("pc = %d, type = %d\n", pc->value(), pc->type());  // c = 7,  type = 2
}

