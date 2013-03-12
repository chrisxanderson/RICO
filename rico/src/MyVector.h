#ifndef __MYVECTOR_H__
#define __MYVECTOR_H__

#include <iostream>
#include <sstream>
#include <string>

#include <vector>

using namespace std;
////////////////////////////////////////////////////////////////////////////////////

template <typename T> class MyVector {
 private:
  vector<T> m_vector;

 public:
  MyVector()                                   {}
  MyVector(const T& x)                         {push_back(x);}
  MyVector(const T& x, const T& y)             {push_back(x); push_back(y);}
  MyVector(const T& x, const T& y, const T& z) {push_back(x); push_back(y); push_back(z);}
  MyVector(const MyVector& v);

 public:
  string      show() const;
  int         size() const {return m_vector.size();}
  vector<T>   get()        {return m_vector;}

  void        push_back (const T& x)           {m_vector.push_back(x);}
  void        operator+=(const T& x)           {m_vector.push_back(x);}
  void        operator+=(const MyVector<T>& v);

  T           operator[](int i)                      const {return m_vector[i];}
  MyVector<T> operator[](const MyVector<bool>& mask) const;

  template<class U> friend ostream& operator<<(ostream& os, const MyVector<U>& rhs);

};

////////////////////////////////////////////////////////////////////////////////////

template<typename T> string MyVector<T>::show() const {
  stringstream sstr;
  sstr << "(";
  for(typename vector<T>::const_iterator ii = m_vector.begin(); ii != m_vector.end(); ii++) {
    sstr << *ii;
    if(ii+1 != m_vector.end()) sstr << ", ";
  }
  sstr << ")";
  return sstr.str();
}

template <typename T> MyVector<T>::MyVector(const MyVector& v) {
  int N = v.size();
  for(int i = 0; i < N; i++) push_back(v[i]);
}

template <typename T> void MyVector<T>::operator+=(const MyVector<T>& v) {
  int N = v.size();
  for(int i = 0; i < N; i++) push_back(v[i]);
}

template <typename T> MyVector<T> MyVector<T>::operator[](const MyVector<bool>& mask) const {
  MyVector<T> result;

  int N = (*this).size();
  for(int i = 0; i < N; i++) if(mask[i]) result += m_vector[i];
  
  return result;
}

template<class U> ostream& operator<<(ostream& os, const MyVector<U>& r) {
  os << r.show();
  return os;
}

////////////////////////////////////////////////////////////////////////////////////
#endif // __MYVECTOR_H__
