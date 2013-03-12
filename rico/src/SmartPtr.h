//
// Generic Refcounted pointer template
//
// This version uses the "indirect container" approach, which is non-
// intrusive, but gives a bit more of a performance hit than having
// an intrusive form which keeps the refcnt in itself.
//
// Example of how to use it:
//   SmartPtr<foo> a = new foo;
// or
//   SmartPtr<foo> a(new foo);
// or
//   SmartPtr<foo> a;
//   foo * f = new foo;
//   a = f;  // Adds to refcnt
//
// NOTE: You must be careful how much you use the last construct, since
// multiple assignments to different SmartPtr pointers of the same
// underlying object pointer can create multiple reference counts on
// the same object, potentially leading to early and incorrect deletion.
//

#ifndef REFCNT_H
#define REFCNT_H

#include <assert.h>

template <class T> class SmartPtr {

 private:
  // Indirect Container, so the refcounts are shared
  struct Cont {
    Cont(T * t) : p(t), cnt(1) { }
    ~Cont()                    { if (p) delete p; }
    T *p;
    unsigned int cnt;
  };
  Cont * c;
  void unassign_cont() {
    if(c) {
      if(!--(c->cnt)) delete c;
      c = 0;
    }
  }

 public:
  SmartPtr() : c(0)         {}
  SmartPtr(T * t)           { c = (t ? (new Cont(t)) : 0); }
  SmartPtr(const SmartPtr& r) { c = r.c; if(c) ++(c->cnt); }
  ~SmartPtr()               { unassign_cont(); }
  
  SmartPtr& operator=(const SmartPtr& r) {
    if(r.c) ++(r.c->cnt);
    if(c)   unassign_cont();
    c = r.c;
    return *this;
  }

  SmartPtr& operator=(const T* t) {
    // For some strict behaviors, if desired
    //assert(!c && "Can only set pointer when not yet assigned");
    if(!t)      
      unassign_cont();
    else if(!c) 
      c = new Cont(t);
    else if(c->p != t) {
      unassign_cont();
      c = new Cont(t);
    }
    return *this;
  }

  T* operator->() const { return (c ? c->p : 0); }
  operator T*()   const { return (c ? c->p : 0); }
};

#endif // REFCNT_H
