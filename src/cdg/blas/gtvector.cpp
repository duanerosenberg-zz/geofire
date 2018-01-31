//==================================================================================
// Module      : gtvector
// Date        : 1/1/2018 (DLR)
// Description : Basic template vector class, provides
//               access to contiguous memory 
// Copyright   : Copyright 2018. Colorado State University. All rights reserved
// Derived from: none.
//==================================================================================
#include <cstddef>
#include <cstring>
#include "gtvector.hpp"

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method
// DESC   : Basic
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<class T> GTVector<T>::GTVector():
data_   (NULLPTR),
n_      (0),
ibeg_   (0),
iend_   (0)
{
}

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method
// DESC   : Instantitate with array size
// ARGS   : GSIZET n: array size
// RETURNS: none
//**********************************************************************************
template<class T> GTVector<T>::GTVector(GSIZET n):
data_   (NULLPTR),
n_      (n),
ibeg_   (0),
iend_   (n-1)
{
  data_ = new T [n_];
  #if defined(_G_AUTO_CREATE_DEV)
    #pragma acc enter data copyin( this[0:1] ) create( data_[0:n_-1] )
  #endif
}

/*
//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method
// DESC   : Instantitate with array size, offsets
// ARGS   : 
// RETURNS: none
//**********************************************************************************
template<class T> GTVector<T>::GTVector(GSIZET n, GSIZET ib, GSIZET ie):
data_   (NULLPTR),
n_      (n),
ibeg_   (ib),
iend_   (ie)
{
  data_ = new T [n_+ibeg_+iend_];
  #if defined(_G_AUTO_CREATE_DEV)
    #pragma acc enter data copyin( this[0:1] ) create( data_[0:n_-1] )
  #endif
}
*/


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method: GTVector argument
// DESC   : Instantitate with GTVector argument
// ARGS   : GTVector<T> argument
// RETURNS: 
//**********************************************************************************
template<class T> GTVector<T>::GTVector(GTVector<T> &obj):
data_   (NULLPTR),
n_      (obj.size()),
ibeg_   (0),
iend_   (obj.size()-1)
{
  data_ = new T [n_];
  std::memcpy(data_, obj.data(), sizeof(T)*n_);
  #if defined(_G_AUTO_CREATE_DEV)
    #pragma acc enter data copyin( this[0:1] ) create( data_[0:n_-1] )
  #endif
}


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method: raw pointer and size args
// DESC   : Instantitate with data block and size
// ARGS   : T *indata: pointer to external data block
//          GSIZET n : size of external data block
// RETURNS: 
//**********************************************************************************
template<class T> GTVector<T>::GTVector(T *indata, GSIZET n):
data_   (NULLPTR),
n_      (n),
ibeg_   (0),
iend_   (n-1)
{
  data_ = new T [n_];
  std::memcpy(data_, indata, sizeof(T)*n_);
  #if defined(_G_AUTO_CREATE_DEV)
    #pragma acc enter data copyin( this[0:1] ) create( data_[0:n_-1] )
  #endif
}

//**********************************************************************************
//**********************************************************************************
// METHOD : Copy constructor method
// DESC   : Override degault copy constructor
// ARGS   : 
// RETURNS: 
//**********************************************************************************
template<class T> GTVector<T>::GTVector(const GTVector<T> &obj):
data_   (NULLPTR),
n_      (obj.size()),
ibeg_   (0),
iend_   (obj.size()-1)
{
  data_ = new T [n_];
  std::memcpy(data_, obj.data(), sizeof(T)*n_);
  #if defined(_G_AUTO_CREATE_DEV)
    #pragma acc enter data copyin( this[0:1] ) create( data_[0:n_-1] )
  #endif
}

//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor
// DESC   :
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
template<class T> GTVector<T>::~GTVector()
{
  #pragma acc exit data delete( data_[0:n_=1], this[0:1] )
  if ( data_ != NULLPTR ) delete [] data_;
}

//**********************************************************************************
//**********************************************************************************
// METHOD : size
// DESC   : Get total size of data block
// ARGS   : 
// RETURNS: GSIZET total size of data block
//**********************************************************************************
template<class T> GSIZET GTVector<T>::size() const
{
  return n_;
} // end, method size

//**********************************************************************************
// METHOD : resize 
// DESC   : Resize vector/datablock
// ARGS   : GSIZET n: new vector size
// RETURNS: none.
//**********************************************************************************
template<class T> void GTVector<T>::resize(GSIZET nnew)
{
  if ( nnew != n_ )
  {
    delete [] data_;
    data_ = new T [nnew];
    n_ = nnew;
  }

  updatedev(); // Update data on device if necessary
 
} // end, method resize

//**********************************************************************************
//**********************************************************************************
// METHOD : data
// DESC   : Get pointer to data block
// ARGS   : 
// RETURNS: T* pointer to data
//**********************************************************************************
template<class T> T *GTVector<T>::data() const
{
  return data_;
} // end, method data 

/*
//**********************************************************************************
//**********************************************************************************
// METHOD : [] access method
// DESC   : Get dereferenced data at specified location
// ARGS   : GSIZET index
// RETURNS: T& data item
//**********************************************************************************
template<class T> T &GTVector<T>::operator[](GSIZET i)
{

#if defined(_G_BOUNDS_CHK)
  const char serr[] = "GTVector<T>::operator[]: ";
  if ( i < ibeg_ || i > iend_ )
  {
    std::cout << serr << "Invalid index" << std::endl;
    exit(1);
  }
#endif
  return data_[i];
} // end, operator []
*/

//**********************************************************************************
//**********************************************************************************
// METHOD : assignment operator= GTVector
// DESC   : Equate to another GTVector
// ARGS   : GTVector<T> & right-hand side arg 
// RETURNS: GTVector & 
//**********************************************************************************
template<class T> GTVector<T> &GTVector<T>::operator=(const GTVector<T> &obj)
{
  if ( this == &obj)
  {
    return *this;
  }
  
  if ( data_ != NULLPTR &&  n_ != obj.size() ) 
  {
    delete [] data_;
    n_ = obj.size();
    data_ = new T[n_];
  }
  std::memcpy(data_, obj.data(), sizeof(T)*n_);
  
  return *this;
} // end, operator=(GTVector &)


//**********************************************************************************
//**********************************************************************************
// METHOD : operator = T
// DESC   : Equate to constant
// ARGS   : T arg
// RETURNS: GTVector & 
//**********************************************************************************
template<class T> GTVector<T> &GTVector<T>::operator=(T a)
{
  for ( GSIZET j=0; j<n_; j++ ) 
  {
    data_[j] = a;
  }

  return *this;
} // end, operator=


//**********************************************************************************
//**********************************************************************************
// METHOD : operator * (1)
// DESC   : point-by-point multiplication
// ARGS   : GTVector
// RETURNS: GTVector & 
//**********************************************************************************
template<class T> GTVector<T> &GTVector<T>::operator*(GTVector<T> &obj)
{
  const char *serr = "GTVector<T>::operator*(1): ";

  if ( obj.size() != n_ )
  {
    std::cout << serr << "incompatible size" << std::endl;
    exit(1);
  }

  T *dobj=obj.data();
  for ( GSIZET j=0; j<n_; j++ ) 
  {
    data_[j] *= dobj[j];
  }

  return *this;
} // end, operator* (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : operator * (2)
// DESC   : multiplication by constant
// ARGS   : T argument
// RETURNS: GTVector & 
//**********************************************************************************
template<class T> GTVector<T> &GTVector<T>::operator*(T b)
{
  const char *serr = "GTVector<T>::operator*(2): ";

  for ( GSIZET j=0; j<n_; j++ ) 
  {
    data_[j] *= b;
  }

  return *this;
} // end, operator* (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : operator + (1)
// DESC   : point-by-point addition
// ARGS   : GTVector
// RETURNS: GTVector & 
//**********************************************************************************
template<class T> GTVector<T> &GTVector<T>::operator+(GTVector<T> &obj)
{
  const char *serr = "GTVector<T>::operator+(1): ";

  if ( obj.size() != n_ )
  {
    std::cout << serr << "incompatible size" << std::endl;
    exit(1);
  }

  GTVector<T> ret(obj.size());
  T *dobj=obj.data();
  T *dret=ret.data();
  for ( GSIZET j=0; j<n_; j++ ) 
  {
    dret[j] = data_[j] + dobj[j];
  }

  return ret;
} // end, operator* (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : operator + (2)
// DESC   : constant addition
// ARGS   : T arg
// RETURNS: GTVector & 
//**********************************************************************************
template<class T> GTVector<T> &GTVector<T>::operator+(T b)
{
  const char *serr = "GTVector<T>::operator+(2): ";

  GTVector<T> ret(n_);
  T *dret=ret.data();
  for ( GSIZET j=0; j<n_; j++ ) 
  {
    ret[j] = data_[j] + b;
  }

  return ret;
} // end, operator+ (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : operator += (1)
// DESC   :
// ARGS   : GTVector &
// RETURNS: void
//**********************************************************************************
template<class T> void  GTVector<T>::operator+=(GTVector<T> &b)
{
  const char *serr = "GTVector<T>::operator+=(1): ";

  T *p = b.data();
  for ( GSIZET j=0; j<n_; j++ ) 
  {
    data_[j] += p[j];
  }

} // end, operator+= (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : operator += (2)
// DESC   :
// ARGS   : T arg
// RETURNS: void
//**********************************************************************************
template<class T> void GTVector<T>::operator+=(T b)
{
  const char *serr = "GTVector<T>::operator+=(2): ";

  for ( GSIZET j=0; j<n_; j++ ) 
  {
    data_[j] += b;
  }

} // end, operator+= (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : operator -= (1)
// DESC   :
// ARGS   : T arg
// RETURNS: void
//**********************************************************************************
template<class T> void GTVector<T>::operator-=(GTVector<T> &b)
{
  const char *serr = "GTVector<T>::operator-=(1): ";

  T *p = b.data();
  for ( GSIZET j=0; j<n_; j++ ) 
  {
    data_[j] -= p[j];
  }

} // end, operatori= (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : operator -= (2)
// DESC   :
// ARGS   : T arg
// RETURNS: void
//**********************************************************************************
template<class T> void GTVector<T>::operator-=(T b)
{
  const char *serr = "GTVector<T>::operator-=(2): ";

  for ( GSIZET j=0; j<n_; j++ ) 
  {
    data_[j] -= b;
  }

} // end, operatori= (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : operator *= (1)
// DESC   : point-by-point product
// ARGS   : GTVector &
// RETURNS: void
//**********************************************************************************
template<class T> void GTVector<T>::operator*=(GTVector<T> &b)
{
  const char *serr = "GTVector<T>::operator*=(1): ";

  T *p = b.data();
  for ( GSIZET j=0; j<n_; j++ ) 
  {
    data_[j] *= p[j];
  }

} // end, operator*= (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : operator *= (2)
// DESC   : product of vector and constant
// ARGS   : GTVector &
// RETURNS: void
//**********************************************************************************
template<class T> void GTVector<T>::operator*=(T b)
{
  const char *serr = "GTVector<T>::operator*=(2): ";

  for ( GSIZET j=0; j<n_; j++ ) 
  {
    data_[j] *= b;
  }

} // end, operator*= (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : operator - (1)
// DESC   : point-by-point subtraction
// ARGS   : GTVector
// RETURNS: GTVector & 
//**********************************************************************************
template<class T> GTVector<T> &GTVector<T>::operator-(GTVector<T> &obj)
{
  const char *serr = "GTVector<T>::operator-(1): ";

  if ( obj.size() != n_ )
  {
    std::cout << serr << "incompatible size" << std::endl;
    exit(1);
  }

  GTVector<T> ret(n_);
  T *dobj=obj.data();
  T *dret=ret.data();
  for ( GSIZET j=0; j<n_; j++ ) 
  {
    dret[j] = data_[j] - dobj[j];
  }

  return ret;
} // end, operator- (1)



//**********************************************************************************
//**********************************************************************************
// METHOD : operator<<
// DESC   : output stream operator
// ARGS   : 
// RETURNS: ostream & 
//**********************************************************************************
template<class T> std::ostream &operator<< (std::ostream &os, const GTVector<T> &obj)
{
  T *d=obj.data_;
  for ( GSIZET j=obj.ibeg_; j<obj.n_-1; j++ )
  {
    os << d[j] << " ";
  } 
  os << d[obj.n_-1] << std::endl;
  return os;
} // end, operator<<

template class GTVector <GFLOAT>;
template class GTVector<GDOUBLE>;
template class GTVector   <GINT>;
template class GTVector <GSIZET>;

