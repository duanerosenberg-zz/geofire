//==================================================================================
// Module       : gtvector.hpp
// Date         : 1/1/18 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a template vector object, with contiguous memory
//                for basic types.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#include "gtypes.h"
#include <cstdlib>
#include <iostream>

template <class T> class GTVector;
template<class T> std::ostream &operator<<(std::ostream &, const GTVector<T> &);


#if !defined(_GTVECTOR_HPP)
#define _GTVECTOR_HPP


template <class T> class GTVector
{
  public:
    GTVector<T>();
    GTVector<T>(GSIZET n);
    GTVector<T>(GTVector<T> &obj);
    GTVector<T>(const GTVector<T> &obj);
    GTVector<T>(T *, GSIZET n);
//  GTVector<T>(GSIZET n, GSIZET ib, GSIZET ie);
   ~GTVector<T>();
    
    T *data() const;
    GSIZET size() const;
    void resize(GSIZET n);
//  void resize(GSIZET n, GSIZET ib, GSIZET ie);

    void updatehost() {
      #pragma acc update self( data_[0:n_-1] )
    };
    void updatedev() {
      #pragma acc update device(data_[0:n_-1])
    };

inline    T &operator[](GSIZET i) {
#if defined(_G_BOUNDS_CHK)
      const char serr[] = "GTVector<T>::operator[]: ";
      if ( i < ibeg_ || i > iend_ )
      {
        std::cout << serr << "Access error: " << i << std::endl;
        exit(1);
      }
#endif
      return data_[i];
    };

inline    T &operator()(GSIZET i) {
#if defined(_G_BOUNDS_CHK)
      const char serr[] = "GTVector<T>::operator(): ";
      if ( i < ibeg_ || i > iend_ )
      {
        std::cout << serr << "Access error: " << i << std::endl;
        exit(1);
      }
#endif
      return data_[i];
    };

inline    T operator[](GSIZET i) const {
#if defined(_G_BOUNDS_CHK)
      const char serr[] = "GTVector<T>::operator[] const: ";
      if ( i < ibeg_ || i > iend_ )
      {
        std::cout << serr << "Access error: " << i << std::endl;
        exit(1);
      }
#endif
      return data_[i];
    };

inline    T operator()(GSIZET i) const {
#if defined(_G_BOUNDS_CHK)
      const char serr[] = "GTVector<T>::operator() const: ";
      if ( i < ibeg_ || i > iend_ )
      {
        std::cout << serr << "Access error: " << i << std::endl;
        exit(1);
      }
#endif
      return data_[i];
    };

    GTVector<T> &operator=(const GTVector<T> &b);
    GTVector<T> &operator=(T b);
    GTVector<T> &operator+(GTVector<T> &b);
    GTVector<T> &operator+(T b);
    void        &operator+=(T b);
    void        &operator+=(GTVector<T> &b);
    GTVector<T> &operator-(GTVector<T> &b);
    GTVector<T> &operator*(GTVector<T> &b);
    GTVector<T> &operator*(T b);
    void        &operator*=(T b);
    void        &operator*=(GTVector<T> &b);

    friend std::ostream &operator<< <> (std::ostream &, const GTVector<T> &);


  private:
    T     *data_;
    GSIZET ibeg_;
    GSIZET iend_;
    GSIZET n_;

};

#include "gtvector.cpp"
#endif
