//==================================================================================
// Module       : gtmatrix.cpp
// Date         : 1/1/18 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a template matrix object, whose data is composed of
//                a regular array of contiguous data, ordered like a Fortran
//                matrix in row-major order (row data changing fastest over column data).
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived from : none.
//==================================================================================
#include <typeinfo>
#include "gtmatrix.hpp"
#include "gcomm.hpp"
#include "mtk.hpp"
#include "lapack_dwrappers_ext.h"

#define GTMATRIX_ROTATE(a,i,j,k,l) g=a(i,j);h=a(k,l);a(i,j)=g-s*(h+g*tau);\
	a(k,l)=h+s*(g-h*tau);

//************************************************************************************
//************************************************************************************
// constructor (1)
//************************************************************************************
template<class T> GTMatrix<T>::GTMatrix()
:
n1_                    (0),
n2_                    (0),
singzero_             (1e-12),
dtype_                (G_GDOUBLE),
cdtype_               (GC_GDOUBLE)
{
  

  data_.resize(1);
  if      ( typeid(T) == typeid(GUSHORT ) ) {
   dtype_  = G_GUSHORT ;
   cdtype_ = GC_GUSHORT ; }
  else if ( typeid(T) == typeid(GSHORT ) ) {
    dtype_  = G_GSHORT ;
    cdtype_ = GC_GSHORT ;}
  else if ( typeid(T) == typeid(GINT ) ) {
    dtype_  = G_GINT ;
    cdtype_ = GC_GINT ; }
  else if ( typeid(T) == typeid(GDOUBLE) ) {
    dtype_  = G_GDOUBLE;
    cdtype_ = GC_GDOUBLE;}
  else if ( typeid(T) == typeid(GQUAD) ) {
    dtype_  = G_GQUAD;
    cdtype_ = GC_GQUAD; }
} // end of constructor 1 



//************************************************************************************
//************************************************************************************
// constructor (2)
//************************************************************************************
template<class T>  GTMatrix<T>::GTMatrix(const GINT   size1, const GINT   size2)
:
n1_                    (size1),
n2_                    (size2),
singzero_             (1e-12),
dtype_                (G_GDOUBLE),
cdtype_               (GC_GDOUBLE)
{

  data_.resize(n1_*n2_);

  if      ( typeid(T) == typeid(GUSHORT ) ) {
   dtype_  = G_GUSHORT ;
   cdtype_ = GC_GUSHORT ; }
  else if ( typeid(T) == typeid(GSHORT ) ) {
    dtype_  = G_GSHORT ;
    cdtype_ = GC_GSHORT ;}
  else if ( typeid(T) == typeid(GINT ) ) {
    dtype_  = G_GINT ;
    cdtype_ = GC_GINT ; }
  else if ( typeid(T) == typeid(GDOUBLE) ) {
    dtype_  = G_GDOUBLE;
    cdtype_ = GC_GDOUBLE;}
  else if ( typeid(T) == typeid(GQUAD) ) {
    dtype_  = G_GQUAD;
    cdtype_ = GC_GQUAD; }

  Zero();
} // end of constructor 2


//************************************************************************************
//************************************************************************************
// constructor (3): square matrix constructor
//************************************************************************************
template<class T>  GTMatrix<T>::GTMatrix(const GINT size1)
:
n1_                    (size1),
n2_                    (size1),
singzero_             (1e-12),
dtype_                (G_GDOUBLE),
cdtype_               (GC_GDOUBLE)
{

  data_.resize(n1_*n2_);

  if      ( typeid(T) == typeid(GUSHORT ) ) {
   dtype_  = G_GUSHORT ;
   cdtype_ = GC_GUSHORT ; }
  else if ( typeid(T) == typeid(GSHORT ) ) {
    dtype_  = G_GSHORT ;
    cdtype_ = GC_GSHORT ;}
  else if ( typeid(T) == typeid(GINT ) ) {
    dtype_  = G_GINT ;
    cdtype_ = GC_GINT ; }
  else if ( typeid(T) == typeid(GDOUBLE) ) {
    dtype_  = G_GDOUBLE;
    cdtype_ = GC_GDOUBLE;}
  else if ( typeid(T) == typeid(GQUAD) ) {
    dtype_  = G_GQUAD;
    cdtype_ = GC_GQUAD; }

  Zero();
} // end of constructor 3


//************************************************************************************
//************************************************************************************
// constructor (4)
//************************************************************************************
template<class T>  GTMatrix<T>::GTMatrix(T *array, GINT   m1, GINT   m2)
:
n1_                    (m1),
n2_                    (m2),
singzero_             (1e-12),
dtype_                (G_GDOUBLE),
cdtype_               (GC_GDOUBLE)
{



  if      ( typeid(T) == typeid(GUSHORT ) ) {
   dtype_  = G_GUSHORT ;
   cdtype_ = GC_GUSHORT ; }
  else if ( typeid(T) == typeid(GSHORT ) ) {
    dtype_  = G_GSHORT ;
    cdtype_ = GC_GSHORT ;}
  else if ( typeid(T) == typeid(GINT ) ) {
    dtype_  = G_GINT ;
    cdtype_ = GC_GINT ; }
  else if ( typeid(T) == typeid(GDOUBLE) ) {
    dtype_  = G_GDOUBLE;
    cdtype_ = GC_GDOUBLE;}
  else if ( typeid(T) == typeid(GQUAD) ) {
    dtype_  = G_GQUAD;
    cdtype_ = GC_GQUAD; }

  // build matrix data_ structure:
  data_.resize(n1_*n2_);
  memcpy(data_.data(), array, n1_*n2_*G_TYPESZ[dtype_]);

  Zero();

} // end of constructor 4


//************************************************************************************
//************************************************************************************
// constructor (5)
//************************************************************************************
template<class T>  GTMatrix<T>::GTMatrix(GTMatrix<T> *m, GINT *ind, GINT  nn)
:
n1_                    (0),
n2_                    (0),
singzero_             (1e-12),
dtype_                (G_GDOUBLE),
cdtype_               (GC_GDOUBLE)
{



  if      ( typeid(T) == typeid(GUSHORT ) ) {
   dtype_  = G_GUSHORT ;
   cdtype_ = GC_GUSHORT ; }
  else if ( typeid(T) == typeid(GSHORT ) ) {
    dtype_  = G_GSHORT ;
    cdtype_ = GC_GSHORT ;}
  else if ( typeid(T) == typeid(GINT ) ) {
    dtype_  = G_GINT ;
    cdtype_ = GC_GINT ; }
  else if ( typeid(T) == typeid(GDOUBLE) ) {
    dtype_  = G_GDOUBLE;
    cdtype_ = GC_GDOUBLE;}
  else if ( typeid(T) == typeid(GQUAD) ) {
    dtype_  = G_GQUAD;
    cdtype_ = GC_GQUAD; }

  // build matrix data structure:
  if ( ind == NULL ) {
    cout << "GTMatrix<T>::GTMatrix (4): NULL reference index set" << endl;
    exit(1);
  }
  n1_ = m->dim(1);
  n2_ = nn;
  data_.resize(n1_*n2_);
  for ( GINT j=0; j<n2_; j++ ) {
    for ( GINT i=0; i<n1_; i++ ) {
      data_[i+j*n1_] = (*m)(i,ind[j]);
    }
  }
  Zero();

} // end of constructor 5


//************************************************************************************
//************************************************************************************
// Copy constructor:
//************************************************************************************
template<class T> GTMatrix<T>::GTMatrix(const GTMatrix<T> &m)
{



  // copy member data:
  n1_        = m.n1_;
  n2_        = m.n2_;
  singzero_ = m.singzero_;
  dtype_    = m.dtype_;
  cdtype_   = m.cdtype_;
  
  data_.resize(n1_*n2_);

  data_ = m.data_;
} // end of copy constructor


//************************************************************************************
//************************************************************************************
// Destructor
//************************************************************************************
template<class T> GTMatrix<T>::~GTMatrix<T>()
{
   DeleteDynamic();
  #pragma acc exit data delete( data_[0:1], this[0:1] )

}

//************************************************************************************
//************************************************************************************
// Assignment operator method
//************************************************************************************
template<class T> GTMatrix<T> &GTMatrix<T>::operator=(const GTMatrix<T> &m)
{
  

  if ( &m != this ) 
  {
    if ( m.n1_ != n1_ || m.n2_ != n2_ )
    {
      cout << "GTMatrix<T>::=: incompatible matrices" << endl;
      while(1);
      exit(1);
    }
    // copy member data:
    n1_        = m.n1_;
    n2_        = m.n2_;
    singzero_ = m.singzero_;
    dtype_    = m.dtype_;
    cdtype_   = m.cdtype_;
    data_      = m.data_;
  }

  return *this;

} // end = operator


//************************************************************************************
//************************************************************************************
// METHOD : operator =
// DESC   : Assignment operator
// ARGS   : typename T
// RETURNS: assigns this to m
//************************************************************************************
template<class T> void  GTMatrix<T>::operator=(T m)
{
  GINT   i;

  for ( i=0; i<n1_*n2_; i++ ) { 
      data_[i] = m; 
  }

} // end = operator


//************************************************************************************
//************************************************************************************
// METHOD : operator *
// DESC   : multiplies this by constant, and returns
//          result, without destroying *this data
// ARGS   :
// RETURNS: product matrix
//************************************************************************************
template<class T> GTMatrix<T> &GTMatrix<T>::operator*(T a) 
{
  GINT          i;
  T             *adata;
  GTMatrix<T>   aprod(n1_,n2_);

  adata = aprod.data();
  memcpy(adata, data_.data(), n1_*n2_*G_TYPESZ[dtype_]);
  for ( i=0; i<n1_*n2_; i++ ) {
    adata[i] *= a; 
   }

  return aprod;

} // end of * operator (for constant ref.)

//************************************************************************************
//************************************************************************************
// METHOD : operator *=
// DESC   : multiplies this by constant, and returns
//          result, modifying *this data
// ARGS   :
// RETURNS: product matrix
//************************************************************************************
template<class T> void GTMatrix<T>::operator*=(T a) 
{
  GINT          i;
  T             *adata;
  GTMatrix<T>   aprod(n1_,n2_);

  adata = aprod.data();
  memcpy(adata, data_.data(), n1_*n2_*G_TYPESZ[dtype_]);
  for ( i=0; i<n1_*n2_; i++ ) {
    adata[i] *= a; 
   }

  return aprod;

} // end of *= operator (for constant ref.)


//************************************************************************************
//************************************************************************************
// METHOD : operator *
// DESC   : matrix-vector product, returns product
//          without destroying *this data
// ARGS   :
// RETURNS: GTMatrix product matrix
//************************************************************************************
template<class T> GTVector<T> GTMatrix<T>::operator*(GTVector<T> &local_a)
{

  GTVector<T>   aret(n1_);

  switch (dtype_) {
    case G_GDOUBLE:
      MTK::fmatvec_prod((GTMatrix<GDOUBLE>&)*this,(GTVector<GDOUBLE>&)local_a,(GTVector<GDOUBLE>&)*aret);
      break;
    case G_GQUAD:
      MTK::qmatvec_prod((GTMatrix<GQUAD>&)*this,(GTVector<GQUAD>&)local_a,(GTVector<GQUAD>&)*aret);
      break;
    default:
      cout << " GTMatrix<T>::operator*(vector): invalid data type" << endl;
      exit(1);
  }
  return aret;

} // end of operator *


//************************************************************************************
//************************************************************************************
// METHOD : operator *
// DESC   : Multiplies this by matrix m, and returns
//          result
// ARGS   :
// RETURNS: GTMatrix product matrix
//************************************************************************************
template<class T> GTMatrix<T> GTMatrix<T>::operator*(GTMatrix<T> m) 
{
  GTMatrix<T> *mret;

  if ( this->n2_ != m.dim(1) ) {
    cout << "GTMatrix<T>::*: (Matrix) incompatible matrix"<< endl;
    exit(1);
  }

  mret = new GTMatrix<T>(n1_,m.dim(2));

  switch (dtype_) {
    case G_GDOUBLE:
      MTK::fmatmat_prod(*((GTMatrix<GDOUBLE>*)this),(GTMatrix<GDOUBLE>&)m,(GTMatrix<GDOUBLE>&)*mret);
      break;
    case G_GQUAD:
//    MTK::qmatmat_prod(*((GTMatrix<GQUAD>*)this),(GTMatrix<GQUAD>&)m,(GTMatrix<GQUAD>&)*mret);
      cout << " GTMatrix<T>::operator*(matrix): Use MTK!" << endl;
      exit(1);
    default:
      cout << " GTMatrix<T>::operator*(matrix): invalid data type" << endl;
      exit(1);
  }
  return *mret;

} // end of operator * (GTMatrix<T>)


//************************************************************************************
//************************************************************************************
// METHOD : operator +=
// DESC   : Matrix addition: this += a, modifying member data
// ARGS   :
// RETURNS: GTMatrix 
//************************************************************************************
template<class T> void GTMatrix<T>::operator+=(GTMatrix<T> &a) 
{

  if ( this->n1_ != a.dim(1) || this->n2_ !=a.dim(2) ) {
    cout << "GTMatrix<T>::+: incompatible matrices"<< endl;
    exit(1);
  }

  switch (dtype_) {
    case G_GDOUBLE:
      MTK::fvec_add((GTVector<GDOUBLE>&)data_,*((GTVector<GDOUBLE>*)a.Vdata()),*((GTVector<GDOUBLE>*)asum->Vdata()));
      break;
    case G_GQUAD:
      MTK::qvec_add((GTVector<GQUAD>&)data,*((GTVector<GQUAD>*)a.Vdata()),*((GTVector<GQUAD>*)(asum->Vdata())));
      exit(1);
    default:
      cout << " GTMatrix<T>::operator+(matrix): invalid data type" << endl;
      exit(1);
  }

}

//************************************************************************************
//************************************************************************************
// METHOD : operator +
// DESC   : Matrix addition: this + a
// ARGS   :
// RETURNS: GTMatrix 
//************************************************************************************
template<class T> GTMatrix<T> &GTMatrix<T>::operator+(GTMatrix<T> &a) 
{

  GTMatrix<T>  *asum;

  if ( this->n1_ != a.dim(1) || this->n2_ !=a.dim(2) ) {
    cout << "GTMatrix<T>::+: incompatible matrices"<< endl;
    exit(1);
  }

  asum  = new GTMatrix<T>(n1_,n2_);

  switch (dtype_) {
    case G_GDOUBLE:
      MTK::fvec_add((GTVector<GDOUBLE>&)data,*((GTVector<GDOUBLE>*)a.Vdata()),*((GTVector<GDOUBLE>*)asum->Vdata()));
      break;
    case G_GQUAD:
      MTK::qvec_add((GTVector<GQUAD>&)data,*((GTVector<GQUAD>*)a.Vdata()),*((GTVector<GQUAD>*)(asum->Vdata())));
      exit(1);
    default:
      cout << " GTMatrix<T>::operator+(matrix): invalid data type" << endl;
      exit(1);
  }

  return *asum;

}

//************************************************************************************
//************************************************************************************
// METHOD : operator -
// DESC   : Matrix subtraction: this - a
// ARGS   :
// RETURNS: GTMatrix 
//************************************************************************************
template<class T> GTMatrix<T> &GTMatrix<T>::operator-(GTMatrix<T> &a) 
{
  GTMatrix<T>  *asum;
  
  if ( this->n1_ != a.dim(1) || this->n2_ !=a.dim(2) ) {
    cout << "GTMatrix<T>::-: incompatible matrices"<< endl;
    exit(1);
  } 

  asum  = new GTMatrix<T>(n1_,n2_);

  switch (dtype_) {
    case G_GDOUBLE:
      MTK::fvec_sub((GTVector<GDOUBLE>&)data,*((GTVector<GDOUBLE>*)a.Vdata()),*((GTVector<GDOUBLE>*)asum->Vdata()));
      break;
    case G_GQUAD:
      MTK::qvec_sub((GTVector<GQUAD>&)data,*((GTVector<GQUAD>*)a.Vdata()),*((GTVector<GQUAD>*)(asum->Vdata())));
      exit(1);
    default:
      cout << " GTMatrix<T>::operator+(matrix): invalid data type" << endl;
      exit(1);
  }

  return *asum;

}


//************************************************************************************
//************************************************************************************
// METHOD : << operator method 
// DESC   : 
// ARGS   :
// RETURNS:  ostream &
//************************************************************************************
template<class T> std::ostream &operator<<(std::ostream &os, GTMatrix<T> &obj)
{

  GINT i, j, N1, N2;
  T    *adata;
   
  adata = obj.data();
  N1    = obj.dim(1);
  N2    = obj.dim(2);
#if 1
  str << endl << "[ ";
  for ( i=0; i<N1; i++ )
  {
    os << "[ ";
    for ( j=0; j<N2-1; j++ )
      os
//        << setiosflags(ios::scientific)
          << setw(16)
          << setprecision(14)
          << setiosflags(ios::fixed)
          << a(i,j)
//        << setw(1)
          << ", ";

      os
//        << setiosflags(ios::scientific)
          << setw(16)
          << setprecision(14)
          << setiosflags(ios::fixed)
          << a(i,j) ;
//        << setw(1) ;
    if ( i < N1-1 ) os << " ]; " << endl;
    else os  << " ]  ]" ;
  }
#else
  os << endl;
  for ( i=0; i<N1; i++ ) {
    for ( j=0; j<N2 ; j++ ) {
      os
#if 1
          << setiosflags(ios::scientific)
//        << setw(18)
//        << setprecision(15)
          << setiosflags(ios::fixed)
#endif
          << (fabs((GQUAD)adata[i+j*N1]) < GTINY*10.0 ? 0.0 : adata[i+j*N1])
          << " ";
    }
    os << endl;
  }
#endif

  return str;
} // end of << operator 


//************************************************************************************
//************************************************************************************
// METHOD : Vdata
// DESC   : 
// ARGS   :
// RETURNS:  reerence to the data container 
//************************************************************************************
template<class T> GTVector<T> &GTMatrix<T>::Vdata() 
{
  return data_;
} // end of method Vdata


//************************************************************************************
//************************************************************************************
// METHOD : DeleteDynamic
// DESC   : Deletes dynamically allocated quantities
// ARGS   :
// RETURNS:  none
//************************************************************************************
template<class T> void GTMatrix<T>::DeleteDynamic()
{
} // end of method DeleteDynamic


//************************************************************************************
//************************************************************************************
// METHOD : resize
// DESC   : resizes dynamically allocated quantities
//          if required
// ARGS   :
// RETURNS:  TRUE on success, else FALSE
//************************************************************************************
template<class T> GBOOL GTMatrix<T>::resize(GINT   new1, GINT   new2)
{
  

  if ( n1_*n2_ == new1*new2 ) return TRUE;

  DeleteDynamic();

  n1_ = new1;
  n2_ = new2;

  data.resize(n1_*n2_);

  updatedev(); // Update data on device if necessary

  return TRUE;

} // end of method resize


//************************************************************************************
//************************************************************************************
// METHOD : resizeM
// DESC   : resizes dynamically allocated quantities
//          if required, with new size an upper limit
// ARGS   :
// RETURNS:  TRUE on success, else FALSE
//************************************************************************************
template<class T> GBOOL GTMatrix<T>::resizeM(GINT   new1, GINT   new2)
{
  
  if ( new1*new2 <= n1_*n2_ )  return TRUE;

  DeleteDynamic();

  n1_ = new1;
  n2_ = new2;

  data_.resize(n1_*n2_);

  return TRUE;

} // end of method resizeM


//************************************************************************************
//************************************************************************************
// METHOD : dim
// DESC   : Array dimension (usable)
//          in direction idir 
// ARGS   :
// RETURNS: GINT   size
//************************************************************************************
template<class T> GINT GTMatrix<T>::dim(GINT   idir) 
{

  if      ( idir == 1 )
    return n1_;
  else if ( idir == 2 )
    return n2_;
  else
    return 0;
} // end of method dim


//************************************************************************************
//************************************************************************************
// METHOD : Zero
// DESC   : Zeros out data elemnts
// ARGS   :
// RETURNS: none
//************************************************************************************
template<class T> void GTMatrix<T>::Zero()
{ 
  data = 0;
}  // end of method Zero



//************************************************************************************
//************************************************************************************
// METHOD : Transpose (1)
// DESC   : Computes transpose of *this, but
//          does not destroy data.
// ARGS   :
// RETURNS: transpose of this
//************************************************************************************
template<class T> GBOOL  GTMatrix<T>::Transpose(GTMatrix<T> &trans)
{

  GINT   i, j, k;
  T      *tdata;

  if ( trans.dim(2) !=  n1_ || trans.dim(1) != n2_ ) {
    cout << "GTMatrix<T>::Transpose: incompatible matrix"<< endl;
    exit(1);
  }

  tdata = trans.data();
  for ( j=0; j<n1_; j++ ) {
    k = j*n2_;
    for ( i=0; i<n2_; i++ ) {
       tdata[i+k] = data[j+i*n1_];
    }
  }
  return TRUE;
 
} // end of method Transpose (1)


//************************************************************************************
//************************************************************************************
// METHOD : Transpose (2)
// DESC   : Computes transpose of *this, but
//          does not destroy data. Computes in box of size nx x ny
// ARGS   :
// RETURNS: transpose of this
//************************************************************************************
template<class T> GBOOL  GTMatrix<T>::Transpose(GTMatrix<T> &trans, GINT nx, GINT ny)
{

  GINT   i, j, k;
  T      *tdata;

  if ( n1_ < nx || n2_ < ny ) {
    cout << "GTMatrix<T>::Transpose: incompatible matrix"<< endl;
    exit(1);
  }

  tdata = trans.data();
  for ( j=0; j<nx; j++ ) {
    k = j*ny;
    for ( i=0; i<ny; i++ ) {
       tdata[i+k] = data[j+i*nx];
    }
  }
  return TRUE;
 
} // end of method Transpose (2)


//************************************************************************************
//************************************************************************************
// METHOD : Transpose (3)
// DESC   : Computes transpose of this, but
//          does not destroy data. A copy is made and
//          returned.
// ARGS   :
// RETURNS: transpose of this
//************************************************************************************
template<class T> GTMatrix<T>  GTMatrix<T>::Transpose()
{

  GTMatrix<T> *t;

  t = new GTMatrix<T>(n1_,n2_);
  if ( !Transpose(*t) ) {
    cout << "GTMatrix<T>::Transpose(3): failed" << endl;
    exit(1);
  }

  return *t;

} // end of method Transpose (3)


//************************************************************************************
//************************************************************************************
// METHOD : Inverse (1)
// DESC   : Computes inverse of this, copying the
//          result to mret
// ARGS   :
// RETURNS: inverse of this
//************************************************************************************
template<class T> GBOOL  GTMatrix<T>::Inverse(GTMatrix<T> &mret)
{

  GINT        i, j, *indx;
  GBOOL       bRet=TRUE;
  char      *serr = "GTMatrix<T>::Inverse(1): ";
  

  if ( mret.dim(1) !=  n1_ || mret.dim(2) != n2_ ) {
    cout << serr << "incompatible matrix"<< endl;
    exit(1);
  }
  if ( n1_ != n2_ ) {
    cout << serr << "matrix not square"<< endl;
    exit(1);
  }

#if 0
  T         **A, **V, *W, *col,  *b, wmax, wmin;
  A     = new T * [n1_];
  V     = new T * [n1_];
  W     = new T   [n1_];
  b     = new T   [n1_];
  col   = new T   [n1_];
  for ( i=0; i<n1_; i++ ) {
     A[i] = new T [n2_];
     V[i] = new T [n2_];
  }
  for ( j=0; j<n2_; j++ ) {
    for ( i=0; i<n1_; i++ ) {
      A[i][j] = (T)data[i+j*n1_];
    }
  }

  if ( !svdcmp(A, n1_, n2_, W, V) ) {
     cout << serr << "svdcmp failed" << endl;
     exit(1);
  }

  // Edit the singular (eigen-)values:
  for ( i=0, wmax=0.0; i<n1_; i++) wmax = MAX(wmax,W[i]); 
  wmin = singzero_*wmax;
  for ( i=0; i<n1_; i++) W[i] = W[i] < wmin ? 0.0 : W[i];

  // Do back substitution to find columns of inverse matrix:
  for ( j=0; j<n2_ && bRet; j++ ) {
    for ( i=0; i<n2_; i++ ) {  b[i] = 0.0; }
    b[j] = 1.0; 
    svbksub(A, W, V, n1_, n2_, b, col);
    for ( i=0; i<n1_; i++ ) mret(i,j) = (T) col[i];
  }

  for ( i=0; i<n1_; i++ ) {
    delete [] A [i];
    delete [] V [i];
  }
  delete [] A;
  delete [] V;
  delete [] W;
  delete [] col;
  delete [] b;
#else
 
  T         **A0, **A, *col,  d, *b;
  mret = 0.0;

  A0    = new T * [n1_];
  A     = new T * [n1_];
  col   = new T   [n1_];
  b     = new T   [n1_];
  indx  = new GINT [n2_];
  for ( i=0; i<n1_; i++ ) {
     A0[i] = new T [n2_];
     A [i] = new T [n2_];
  }
  for ( j=0; j<n2_; j++ ) {
    for ( i=0; i<n1_; i++ ) {
      A0[i][j] = (T)data[i+j*n1_];
      A [i][j] = A0[i][j];
    }
  }

  if ( !wludcmp(A, n1_, indx, &d) ) {
     cout << serr << "wludcmp failed" << endl;
     bRet = FALSE;
  }

  for ( j=0; j<n2_ && bRet; j++ ) {
    for ( i=0; i<n1_; i++ ) {  col[i] = 0.0; b[i] = 0.0; }
    col[j] = 1.0; b[j] = col[j];
    if ( !(bRet=lubksb(A, n1_, indx, col)) ) {
       cout << serr << "lubjsb failed" << endl;
       bRet = FALSE;
       break;
    }
    if ( !(bRet=improve(A0, A, n1_, indx, b, col)) ) {
       cout << serr << "improve failed" << endl;
       bRet = FALSE;
       break;
    }
    for ( i=0; i<n1_; i++ ) mret(i,j) = (T) col[i];
  }

  for ( i=0; i<n1_; i++ ) {
    delete [] A0[i];
    delete [] A [i];
  }
  delete [] A0;
  delete [] A;
  delete [] col;
  delete [] b;
  delete [] indx;
 #endif

  return bRet;

} // end of method Inverse (1)


//************************************************************************************
//************************************************************************************
// METHOD : Inverse (2)
// DESC   : Computes inverse of this, copying the result to mret;
//          Inverse is made on a box of size nx x ny.
// ARGS   :
// RETURNS: inverse of this
//************************************************************************************
template<class T> GBOOL  GTMatrix<T>::Inverse(GTMatrix<T> &mret, GINT nx, GINT ny)
{

  GINT        i, j, *indx;
  GBOOL       bRet=TRUE;
  T         **A0, **A, *col,  d, *b;
  char       *serr = "GTMatrix<T>::Inverse(2): ";
  

  if ( n1_ < nx || n2_ < ny ) {
    cout << serr << "incompatible matrix"<< endl;
    exit(1);
  }
  if ( nx != ny ) {
    cout << serr << "matrix not square"<< endl;
    exit(1);
  }

  mret       = 0.0;

  A0    = new T * [nx];
  A     = new T * [nx];
  col   = new T   [nx];
  b     = new T   [nx];
  indx  = new GINT [ny];
  for ( i=0; i<nx; i++ ) {
     A0[i] = new T [ny];
     A [i] = new T [ny];
  }
  for ( j=0; j<ny; j++ ) {
    for ( i=0; i<nx; i++ ) {
      A0[i][j] = (T)data[i+j*nx];
      A [i][j] = A0[i][j];
    }
  }

  if ( !wludcmp(A, nx, indx, &d) ) {
     cout << serr << "wludcmp failed" << endl;
     bRet = FALSE;
  }

  for ( j=0; j<ny && bRet; j++ ) {
    for ( i=0; i<nx; i++ ) {  col[i] = 0.0; b[i] = 0.0; }
    col[j] = 1.0; b[j] = col[j];
    if ( !(bRet=lubksb(A, nx, indx, col)) ) {
       cout << serr << "lubjsb failed" << endl;
       bRet = FALSE;
       break;
    }
    if ( !(bRet=improve(A0, A, nx, indx, b, col)) ) {
       cout << serr << "improve failed" << endl;
       bRet = FALSE;
       break;
    }
    for ( i=0; i<nx; i++ ) mret(i,j) = (T) col[i];
  }

  for ( i=0; i<nx; i++ ) {
    delete [] A0[i];
    delete [] A [i];
  }
  delete [] A0;
  delete [] A;
  delete [] col;
  delete [] b;
  delete [] indx;

  return bRet;

} // end of method Inverse (2)


//************************************************************************************
//************************************************************************************
// METHOD : Inverse (2)
// DESC   : Computes inverse of this, but
//          does not destroy data. A copy is made and
//          returned.
// ARGS   :
// RETURNS: inverse of this
//************************************************************************************
template<class T> GTMatrix<T>  GTMatrix<T>::Inverse()
{

  GTMatrix<T> *mi;

  mi = new GTMatrix<T>(n1_,n2_);

  if ( !Inverse(*mi) ) {
    cout << "GTMatrix<T>::Inverse(2): failed" << endl;
    exit(1);
  }

  return *mi;

} // end of method Inverse (2)


//************************************************************************************
//************************************************************************************
// METHOD : isSymmetric 
// DESC   : determines if matrix is symmetric
// ARGS   :
// RETURNS: TRUE or FALSE 
//************************************************************************************
template<class T> GBOOL  GTMatrix<T>::isSymmetric()
{
  GINT   i, j, k, m;
  GBOOL  bRet;

  if ( n1_ != n2_ ) return FALSE;

  // NOTE: should be symmetric w.r.t some tolerance!!!
  for ( j=1; j<n2_; j++ ) {
    m = j*n2_;
    for ( i=j+1,bRet=TRUE; i<n1_-1; i++ ) {
      k = n2_*i;
      bRet = bRet && data[i+m] == data[i+k];
    }
  }
  return bRet;
 
} // end of method isSymmetric


//************************************************************************************
//************************************************************************************
// METHOD : improve
// DESC   : Taken largely from Numerical Recipes in C++, p. 59: Iterate to improve the
//          solution to Ax = b for x, given A, its LU decomposition, and b, and the 
//          initial solution, x
// ARGS   :
// RETURNS: GBOOL TRUE on success; else FALSE
//************************************************************************************
template<class T> GBOOL  GTMatrix<T>::improve(T **&a, T **alud, GINT n, GINT *&indx, T b[], T x[])
{
  T    sdp, *r=new T[n];
  GINT i, j, m=0, nloop=1;
  if ( a == NULL || alud == NULL || indx == NULL ) {
    delete [] r; 
    return FALSE;
  }
  while ( m < nloop ) {
    for ( i=0; i<n; i++ ) {
      sdp = -b[i];
      for ( j=0; j<n; j++ ) 
        sdp += (GQUAD)a[i][j] * (GQUAD)x[j];
      r[i] = sdp;
    }  
    lubksb(alud, n, indx, r);
    for ( i=0; i<n; i++ ) x[i] -= r[i];
   m++;
  }
    
  delete [] r; 
  return TRUE;
} // end of method improve


//************************************************************************************
//************************************************************************************
// METHOD : wludcmp
// DESC   : Taken largely from Numerical Recipes
// ARGS   :
// RETURNS: GBOOL flag
//************************************************************************************
template<class T> GBOOL  GTMatrix<T>::wludcmp(T **&a, GINT   n, GINT   *&indx, T *d)
{
  if ( a == NULL || indx == NULL ) return FALSE;

  GINT  i, imax, j, k;
  GBOOL bRet=TRUE;
  T     big, dum, sum, temp, gtiny=1.0e-20;
  T     *vv;

  vv = new T [n]; 
  *d=1.0;
  for ( i=0; i<n && bRet; i++ ) {
    big=0.0;
    for ( j=0; j<n; j++ )
      if ( (temp=fabs((GQUAD)a[i][j])) > big ) big=temp;
    if ( big == 0.0 ) {
      bRet = FALSE; 
      break;
    }
    vv[i]=1.0 / big;
  }

  for ( j=0; j<n && bRet; j++ ) {
    for ( i=0; i<j; i++ ) {
      sum=a[i][j];
      for ( k=0; k<i; k++ ) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for ( i=j; i<n; i++ ) {
      sum=a[i][j];
      for ( k=0; k<j; k++ ) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs((GQUAD)sum)) >= big ) {
        big=dum;
        imax=i;
      }
    }
    if ( j != imax ) {
      for ( k=0; k<n; k++ ) {
        dum=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
     *d = -(*d);
     vv[imax]=vv[j];
    }
    indx[j]=imax;
    if ( fabs((GQUAD)a[j][j]) < gtiny ) a[j][j]=gtiny;
//  if ( a[j][j] == 0.0 ) a[j][j]=gtiny;
    if ( j != n-1 ) {
      dum=1.0/(a[j][j]);
      for ( i=j+1; i<n; i++ ) a[i][j] *= dum;
    }
  }
  delete [] vv;


  return bRet;

} // end of method wludcmp


//************************************************************************************
//************************************************************************************
// METHOD : lubksb
// DESC   : Taken largely from Numerical Recipes
// ARGS   :
// RETURNS: GBOOL flag
//************************************************************************************
template<class T> GBOOL  GTMatrix<T>::lubksb(T **&a, GINT   nd, GINT   *&indx, T b[])
{

    GINT  n=nd-1; 
    GINT  i, ii = -1, ip, j;
    T     sum;

#if 0
    for ( i=0; i <= n; i++ ) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii > -1) {
            for (j = ii; j < i; j++) sum -= a[i][j]*b[j];
        }
        else { 
            if (sum) ii = i;
        }
        b[i] = sum;
    }
    for ( i=n; i>=0; i-- ) {
        sum=b[i];
        if ( i < n )
        for (j = i+1; j <= n; j++) sum -= a[i][j]*b[j];
        b[i] = sum /(a[i][i]);
    }
#endif
    n = nd; 
    ii = 0;
    for ( i=0; i<n; i++ ) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if ( ii != 0 ) {
            for ( j=ii-1; j<i; j++ ) sum -= a[i][j]*b[j];
        }
        else  if ( sum != 0.0 ) ii = i+1;
        b[i] = sum;
    }
    for ( i=n-1; i>=0; i-- ) {
        sum=b[i];
        for ( j=i+1; j<n; j++ ) sum -= a[i][j]*b[j];
        b[i] = sum /(a[i][i]);
    }


  return TRUE;

} // end of method lubksb


//************************************************************************************
//************************************************************************************
// METHOD : ludcmp
// DESC   : LU decomposition method provided by
//          Warren Jasper, NC State Univ.
// ARGS   :
// RETURNS: GBOOL flag
//************************************************************************************
template<class T> GBOOL  GTMatrix<T>::ludcmp(T **&a, GINT   nd, GINT   *&indx, T *d)
{
    GINT   i, imax, j, k;
    GINT   n=nd-1;
    GBOOL bRet = TRUE;
    T     big, dum, sum, temp, gtiny=1.0e-20;
    T     *vv = new T [nd];

    *d = 1.0;
    imax = -1;

    for ( i=0; i<=n; i++ ) {
        big = 0.0;
        for ( j=0; j<=n; j++ )
            if ( (temp = fabs((GQUAD) a[i][j] ) ) > big) big = temp;
            if ( big == 0.0 ) {
              bRet = FALSE;  
              cout << "GTMatrix::ludcmp: big = 0" << endl;
              break;
            }
        vv[i] = 1.0 / big;
    }
    for ( j=0; j<=n; j++ ) {
        for ( i=0; i<j; i++ ) {
            sum = a[i][j];
                for ( k=0; k<i; k++ ) sum -= a[i][k]*a[k][j];
                a[i][j] = sum;
        }
        big = 0.0;
        for ( i=j; i<=n; i++ ) {
            sum = a[i][j];
            for ( k=0; k<j; k++ )
                sum -= a[i][k]*a[k][j];
            a[i][j] = sum;
            if ( (dum = vv[i]*fabs((GQUAD)sum)) >= big ) {
                big = dum;
                imax = i;
            }
        }
        if (j != imax) {
            if ( imax < 0 ) { cout << "GTMatrix::ludcmp: imax < 0 " << endl; return FALSE;}
            for ( k=0; k<=n; k++ ) {
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            *d = -(*d);
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if ( a[j][j] == 0.0 ) a[j][j] = gtiny;
        if ( j != n ) {
            dum = 1.0 / (a[j][j]);
            for ( i=j+1; i<=n; i++ ) a[i][j] *= dum;
        }
    }
    delete [] vv;

    return bRet;

} // end of method ludcmp

//************************************************************************************
//************************************************************************************
// METHOD : isamax
// DESC   : BLAS routine that finds the index of element having max.
//          absolute value.
// ARGS   : n   : Number of elements to check.
//          sx  : Vector to be checked.
//          incx: Every incx-th element is checked.
// RETURNS: int
//************************************************************************************
GINT  isamax( GINT  n, GDOUBLE *sx, GINT  incx )
{
  GDOUBLE smax = 0.0e0;
  GINT   i, istmp = 0;

  if( n <= 1 ) return( istmp );
  if( incx != 1 ) {
    /* Code for increment not equal to 1. */
    if( incx < 0 ) sx = sx + ((-n+1)*incx + 1);
    istmp = 0;
    smax  = fabs((GQUAD) *sx );
    sx += incx;
    for( i=1; i<n; i++, sx+=incx )
      if( fabs((GQUAD) *sx ) > smax ) {
        istmp = i;
        smax  = fabs((GQUAD) *sx );
      }
    return( istmp );
  }
  /* Code for increment equal to 1. */
  istmp = 0;
  smax  = fabs((GQUAD)*sx);
  sx++;
  for( i=1; i<n; i++, sx++ )
    if( fabs((GQUAD) *sx ) > smax ) {
      istmp = i;
      smax  = fabs((GQUAD) *sx );
    }
  return( istmp );
}  // end of method isamax


//************************************************************************************
//************************************************************************************
// METHOD : svdcmp (1)
// DESC   : Taken from Numerical Recipes, p.70-72, find 
//          solution to Ax = b for x, given A, using Singular Value Decomp,
//          s.t. A = U W V^T. U replaces A on output. W is a diagonal
//          matrix of singular values represented as a vector in [0, ..., n-1].
//          V (not V^T) is output as an array v[0,...,n-1][0,...,n-1]
// ARGS   : 
// RETURNS: TRUE on success; else FALSE
//************************************************************************************
template<class T> GBOOL GTMatrix<T>::svdcmp(T **a, GINT m, GINT n, T w[], T **v)
{
	GINT    i,its,j,jj,k,l, niter=50, nm;
        GBOOL   bRet=TRUE, flag;
	T       anorm,c,f,g,h,s,scale,x,y,z,*rv1;
        char   *serr = " GBOOL GTMatrix<T>::svdcmp (1): ";

        rv1 = new T [n];
	g=scale=anorm=0.0;
	for (i=0;i<n;i++) {
		l=i+2;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i < m) {
			for (k=i;k<m;k++) scale += fabs((GQUAD)a[k][i]);
			if (scale != 0.0) {
				for (k=i;k<m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -GSIGN(sqrt((GQUAD)s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l-1;j<n;j++) {
					for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i+1 <= m && i != n) {
			for (k=l-1;k<n;k++) scale += fabs((GQUAD)a[i][k]);
			if (scale != 0.0) {
				for (k=l-1;k<n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l-1];
				g = -GSIGN(sqrt((GQUAD)s),f);
				h=f*g-s;
				a[i][l-1]=f-g;
				for (k=l-1;k<n;k++) rv1[k]=a[i][k]/h;
				for (j=l-1;j<m;j++) {
					for (s=0.0,k=l-1;k<n;k++) s += a[j][k]*a[i][k];
					for (k=l-1;k<n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l-1;k<n;k++) a[i][k] *= scale;
			}
		}
		anorm=MAX(anorm,(fabs((GQUAD)w[i])+fabs((GQUAD)rv1[i])));
	}
//
	for (i=n-1;i>=0;i--) {
		if (i < n-1) {
			if (g != 0.0) {
				for (j=l;j<n;j++) v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=MIN(m,n)-1;i>=0;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<n;j++) a[i][j]=0.0;
		if (g != 0.0) {
			g=1.0/g;
			for (j=l;j<n;j++) {
//
				for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<m;j++) a[j][i] *= g;
		} else for (j=i;j<m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n-1;k>=0;k--) {
		for (its=0;its<niter;its++) {
//-----
			flag=TRUE;
			for (l=k;l>=0;l--) {
				nm=l-1;
				if ((GQUAD)(fabs((GQUAD)rv1[l])+anorm) == anorm) {
					flag=FALSE;
					break;
				}
				if ((GQUAD)(fabs((GQUAD)w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l-1;i<k+1;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((GQUAD)(fabs((GQUAD)f)+anorm) == anorm) break;
					g=w[i];
					h=dpythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=0;j<n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == (niter-1)) {
                           cout << serr << "no convergence in " << niter << " iterations" << endl;
                           return FALSE;
                        }
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=dpythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+GSIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=dpythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=dpythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
        delete [] rv1;

  return bRet;
} // end of method svdcmp (1)



//************************************************************************************
//************************************************************************************
// METHOD : svdcmp (2)
// DESC   : Taken from Numerical Recipes, p.70-72, find 
//          solution to Ax = b for x, given A, using Singular Value Decomp,
//          s.t. A = U W V^T. U replaces A on output. W is a diagonal
//          matrix of singular values represented as a vector in [0, ..., n-1].
//          V (not V^T) is output as an array v[0,...,n-1][0,...,n-1]
// ARGS   : 
//          w   : vector in which to store e-values; size of at least dim(diag(this))
//          v   : matrix of same size as a containing values of V
//          rv1 : GVector of dimension diag(this) 
// RETURNS: TRUE on success; else FALSE
//************************************************************************************
template<class T> GBOOL GTMatrix<T>::svdcmp(GTVector<T> &w, GTMatrix<T> &v,  GTVector<T> &rv1)
{
	GINT    i,its,j,jj,k,l, m, n, niter=50, nm;
        GBOOL   bRet=TRUE, flag;
	T       anorm,c,f,g,h,s,scale,x,y,z;
        char   *serr = " GBOOL GTMatrix<T>::svdcmp (2): ";

        m = this->dim(1);
        n = this->dim(2);
	g=scale=anorm=0.0;
	for (i=0;i<n;i++) {
		l=i+2;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i < m) {
			for (k=i;k<m;k++) scale += fabs((GQUAD)(*this)(k,i)); 
			if (scale != 0.0) {
				for (k=i;k<m;k++) {
					(*this)(k,i) /= scale;
					s += (*this)(k,i)*(*this)(k,i);
				}
				f=(*this)(i,i);
				g = -GSIGN(sqrt((GQUAD)s),f);
				h=f*g-s;
				(*this)(i,i)=f-g;
				for (j=l-1;j<n;j++) {
					for (s=0.0,k=i;k<m;k++) s += (*this)(k,i)*(*this)(k,j);
					f=s/h;
					for (k=i;k<m;k++) (*this)(k,j) += f*(*this)(k,i);
				}
				for (k=i;k<m;k++) (*this)(k,i) *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i+1 <= m && i != n) {
			for (k=l-1;k<n;k++) scale += fabs((GQUAD)(*this)(i,k));
			if (scale != 0.0) {
				for (k=l-1;k<n;k++) {
					(*this)(i,k) /= scale;
					s += (*this)(i,k)*(*this)(i,k);
				}
				f=(*this)(i,l-1);
				g = -GSIGN(sqrt((GQUAD)s),f);
				h=f*g-s;
				(*this)(i,l-1)=f-g;
				for (k=l-1;k<n;k++) rv1[k]=(*this)(i,k)/h;
				for (j=l-1;j<m;j++) {
					for (s=0.0,k=l-1;k<n;k++) s += (*this)(j,k)*(*this)(i,k);
					for (k=l-1;k<n;k++) (*this)(j,k) += s*rv1[k];
				}
				for (k=l-1;k<n;k++) (*this)(i,k) *= scale;
			}
		}
		anorm=MAX(anorm,(fabs((GQUAD)w[i])+fabs((GQUAD)rv1[i])));
	}
//
	for (i=n-1;i>=0;i--) {
		if (i < n-1) {
			if (g != 0.0) {
				for (j=l;j<n;j++) v(j,i)=((*this)(i,j)/(*this)(i,l))/g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += (*this)(i,k)*v(k,j);
					for (k=l;k<n;k++) v(k,j) += s*v(k,i);
				}
			}
			for (j=l;j<n;j++) v(i,j)=v(j,i)=0.0;
		}
		v(i,i)=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=MIN(m,n)-1;i>=0;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<n;j++) (*this)(i,j)=0.0;
		if (g != 0.0) {
			g=1.0/g;
			for (j=l;j<n;j++) {
//
				for (s=0.0,k=l;k<m;k++) s += (*this)(k,i)*(*this)(k,j);
				f=(s/(*this)(i,i))*g;
				for (k=i;k<m;k++) (*this)(k,j) += f*(*this)(k,i);
			}
			for (j=i;j<m;j++) (*this)(j,i) *= g;
		} else for (j=i;j<m;j++) (*this)(j,i)=0.0;
		++(*this)(i,i);
	}
	for (k=n-1;k>=0;k--) {
		for (its=0;its<niter;its++) {
//-----
			flag=TRUE;
			for (l=k;l>=0;l--) {
				nm=l-1;
				if ((GQUAD)(fabs((GQUAD)rv1[l])+anorm) == anorm) {
					flag=FALSE;
					break;
				}
				if ((GQUAD)(fabs((GQUAD)w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l-1;i<k+1;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((fabs((GQUAD)f)+anorm) == anorm) break;
					g=w[i];
					h=dpythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y=(*this)(j,nm);
						z=(*this)(j,i);
						(*this)(j,nm)=y*c+z*s;
						(*this)(j,i)=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=0;j<n;j++) v(j,k) = -v(j,k);
				}
				break;
			}
			if (its == (niter-1)) {
                           cout << serr << "no convergence in " << niter << " iterations" << endl;
                           return FALSE;
                        }
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=dpythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+GSIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=dpythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x=v(jj,j);
					z=v(jj,i);
					v(jj,j)=x*c+z*s;
					v(jj,i)=z*c-x*s;
				}
				z=dpythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<m;jj++) {
					y=(*this)(jj,j);
					z=(*this)(jj,i);
				        (*this)(jj,j)=y*c+z*s;
					(*this)(jj,i)=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}

  return bRet;
} // end of method svdcmp (2)


//************************************************************************************
//************************************************************************************
// METHOD : svdcmp (3)
// DESC   : Taken from Numerical Recipes, p.70-72, find 
//          solution to Ax = b for x, given A, using Singular Value Decomp,
//          s.t. A = U W V^T. U replaces A on output. W is a diagonal
//          matrix of singular values represented as a vector in [0, ..., n-1].
//          V (not V^T) is output as an array v[0,...,n-1][0,...,n-1]
// ARGS   : 
//          w    : vector in which to store e-values; size of at least dim(diag(this))
//          v    : matrix of same size as a containing values of V
//          rv1  : GVector of dimension diag(this) 
//          nx,ny: valid dims of matrix on which to operate; may be less than nominal
//                 matrix dimensions. 
// RETURNS: TRUE on success; else FALSE
//************************************************************************************
template<class T> GBOOL GTMatrix<T>::svdcmp(GTVector<T> &w, GTMatrix<T> &v,  GTVector<T> &rv1, GINT nx, GINT ny)
{
	GINT    i,its,j,jj,k,l, m, n, niter=50, nm;
        GBOOL   bRet=TRUE, flag;
	T       anorm,c,f,g,h,s,scale,x,y,z;
        char   *serr = " GBOOL GTMatrix<T>::svdcmp (2): ";

        m = nx;
        n = ny;
	g=scale=anorm=0.0;
	for (i=0;i<n;i++) {
		l=i+2;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i < m) {
			for (k=i;k<m;k++) scale += fabs((GQUAD)(*this)(k,i)); 
			if (scale != 0.0) {
				for (k=i;k<m;k++) {
					(*this)(k,i) /= scale;
					s += (*this)(k,i)*(*this)(k,i);
				}
				f=(*this)(i,i);
				g = -GSIGN(sqrt((GQUAD)s),f);
				h=f*g-s;
				(*this)(i,i)=f-g;
				for (j=l-1;j<n;j++) {
					for (s=0.0,k=i;k<m;k++) s += (*this)(k,i)*(*this)(k,j);
					f=s/h;
					for (k=i;k<m;k++) (*this)(k,j) += f*(*this)(k,i);
				}
				for (k=i;k<m;k++) (*this)(k,i) *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i+1 <= m && i != n) {
			for (k=l-1;k<n;k++) scale += fabs((GQUAD)(*this)(i,k));
			if (scale != 0.0) {
				for (k=l-1;k<n;k++) {
					(*this)(i,k) /= scale;
					s += (*this)(i,k)*(*this)(i,k);
				}
				f=(*this)(i,l-1);
				g = -GSIGN(sqrt((GQUAD)s),f);
				h=f*g-s;
				(*this)(i,l-1)=f-g;
				for (k=l-1;k<n;k++) rv1[k]=(*this)(i,k)/h;
				for (j=l-1;j<m;j++) {
					for (s=0.0,k=l-1;k<n;k++) s += (*this)(j,k)*(*this)(i,k);
					for (k=l-1;k<n;k++) (*this)(j,k) += s*rv1[k];
				}
				for (k=l-1;k<n;k++) (*this)(i,k) *= scale;
			}
		}
		anorm=MAX(anorm,(fabs((GQUAD)w[i])+fabs((GQUAD)rv1[i])));
	}
//
	for (i=n-1;i>=0;i--) {
		if (i < n-1) {
			if (g != 0.0) {
				for (j=l;j<n;j++) v(j,i)=((*this)(i,j)/(*this)(i,l))/g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += (*this)(i,k)*v(k,j);
					for (k=l;k<n;k++) v(k,j) += s*v(k,i);
				}
			}
			for (j=l;j<n;j++) v(i,j)=v(j,i)=0.0;
		}
		v(i,i)=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=MIN(m,n)-1;i>=0;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<n;j++) (*this)(i,j)=0.0;
		if (g != 0.0) {
			g=1.0/g;
			for (j=l;j<n;j++) {
//
				for (s=0.0,k=l;k<m;k++) s += (*this)(k,i)*(*this)(k,j);
				f=(s/(*this)(i,i))*g;
				for (k=i;k<m;k++) (*this)(k,j) += f*(*this)(k,i);
			}
			for (j=i;j<m;j++) (*this)(j,i) *= g;
		} else for (j=i;j<m;j++) (*this)(j,i)=0.0;
		++(*this)(i,i);
	}
	for (k=n-1;k>=0;k--) {
		for (its=0;its<niter;its++) {
//-----
			flag=TRUE;
			for (l=k;l>=0;l--) {
				nm=l-1;
				if ((GQUAD)(fabs((GQUAD)rv1[l])+anorm) == anorm) {
					flag=FALSE;
					break;
				}
				if ((GQUAD)(fabs((GQUAD)w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l-1;i<k+1;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((fabs((GQUAD)f)+anorm) == anorm) break;
					g=w[i];
					h=dpythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y=(*this)(j,nm);
						z=(*this)(j,i);
						(*this)(j,nm)=y*c+z*s;
						(*this)(j,i)=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=0;j<n;j++) v(j,k) = -v(j,k);
				}
				break;
			}
			if (its == (niter-1)) {
                           cout << serr << "no convergence in " << niter << " iterations" << endl;
                           return FALSE;
                        }
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=dpythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+GSIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=dpythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x=v(jj,j);
					z=v(jj,i);
					v(jj,j)=x*c+z*s;
					v(jj,i)=z*c-x*s;
				}
				z=dpythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<m;jj++) {
					y=(*this)(jj,j);
					z=(*this)(jj,i);
				        (*this)(jj,j)=y*c+z*s;
					(*this)(jj,i)=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}

  return bRet;
} // end of method svdcmp (3)



//************************************************************************************
//************************************************************************************
// METHOD : jacobi (1)
// DESC   : Taken from Numerical Recipes, p.346-348, find 
//          eigen-decomposition for square matrix. Operates on *this.
// ARGS   : 
//          d   : vector of eigenvalues, returned
//          v   : matrix of same size as *this whose columns contain normalized e-vectors, returned
//                nrot: no. Jacobi rotations used, returned 
//          a   : tmp matrix of same size as *this
//          b   : tmp vector of same size as diagnonal of *this
//          z   : tmp vector of same size as diagnonal of *this
// RETURNS: TRUE on success; else FALSE
//************************************************************************************
 template<class T> void GTMatrix<T>::jacobi(GTVector<T> &d, GTMatrix<T> &v, GINT &nrot, GTMatrix<T> &a, GTVector<T> &b, GTVector<T> &z)
{
	GINT    j,iq,ip,i,m,n;
	GDOUBLE tresh,theta,tau,t,sm,s,h,g,c;
        char   *serr = "GTMatrix<T>::jacobi (1): ";

        m = this->dim(1);
        n = this->dim(2);
        if ( m != n ) {
          cout << serr << "matrix must be square" << endl;
          exit(1);
        }

        a = *this; // make copy

	for (ip=0;ip<n;ip++) {
		for (iq=0;iq<n;iq++) v(ip,iq)=0.0;
		v(ip,ip)=1.0;
	}
	for (ip=0;ip<n;ip++) {
		b[ip]=d[ip]=a(ip,ip);
		z[ip]=0.0;
	}
	nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++)
				sm += fabs((GQUAD)a(ip,iq));
		}
		if (sm == 0.0) {
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++) {
				g=100.0*fabs((GQUAD)a(ip,iq));
				if (i > 4 && (float)(fabs((GQUAD)d[ip])+g) == (float)fabs((GQUAD)d[ip])
					&& (float)(fabs((GQUAD)d[iq])+g) == (float)fabs((GQUAD)d[iq]))
					a(ip,iq)=0.0;
				else if (fabs((GQUAD)a(ip,iq)) > tresh) {
					h=d[iq]-d[ip];
					if ((float)(fabs((GQUAD)h)+g) == (float)fabs((GQUAD)h))
						t=(a(ip,iq))/h;
					else {
						theta=0.5*h/(a(ip,iq));
						t=1.0/(fabs((GQUAD)theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a(ip,iq);
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a(ip,iq)=0.0;
					for (j=0;j<ip;j++) {
						GTMATRIX_ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<iq;j++) {
						GTMATRIX_ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<n;j++) {
						GTMATRIX_ROTATE(a,ip,j,iq,j)
					}
					for (j=0;j<n;j++) {
						GTMATRIX_ROTATE(v,j,ip,j,iq)
					}
					++nrot;
				}
			}
		}
		for (ip=0;ip<n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	cout << serr << "Too many iterations" << endl;

} // end of method jacobi (1)


//************************************************************************************
//************************************************************************************
// METHOD : SetSingularZero
// DESC   : Serves as public method set tolerance at which the singular (eigen-) values
//          from SVD method are zeroed out.
// ARGS   : T  tol 
// RETURNS: none
//************************************************************************************
template<class T> void GTMatrix<T>::SetSingularZero(T tol)
{
  singzero_ = tol;

} // end of method SetSingularZero


//************************************************************************************
//************************************************************************************
// METHOD : dpythag
// DESC   : Taken from Numerical Recipes and serves as a utility for svdcmp method
//          Computes (a^2 + b^2)^1/2 without underflow or overflow.
// ARGS   : 
// RETURNS: T result
//************************************************************************************
template<class T> T GTMatrix<T>::dpythag(T a, T b)
{
	T    absa,absb;
	absa=fabs((GQUAD)a);
	absb=fabs((GQUAD)b);
	if (absa > absb) return absa*sqrt(1.0+pow((GDOUBLE)(absb/absa),2.0));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+pow((GDOUBLE)(absa/absb),2.0)));
} // end of method dpythag


//************************************************************************************
//************************************************************************************
// METHOD     : svbksub
// DESC   : Taken from Numerical Recipes and serves back substitution 
//              method for SVD
// ARGS     : 
// RETURNS    : T result
//************************************************************************************
template<class T> void GTMatrix<T>::svbksub(T **u, T w[], T **v, GINT n, GINT m, T b[], T x[])
{
        GINT     jj,j,i;
        T        s,*tmp;

        tmp = new T [n];
        for (j=0;j<n;j++) {
                s=0.0;
                if (w[j] != 0.0) {
                        for (i=0;i<m;i++) s += u[i][j]*b[i];
                        s /= w[j];
                }
                tmp[j]=s;
        }
        for (j=0;j<n;j++) {
                s=0.0;
                for (jj=0;jj<n;jj++) s += v[j][jj]*tmp[jj];
                x[j]=s;
        }
        delete [] tmp;
} // end of method svbksub


//************************************************************************************
//************************************************************************************
// METHOD : choldc (1)
// DESC   : Taken from Numerical Recipes: serves as method to perform
//          efficient Cholesky decomposition of an input matrix, a, s.t.
//          a = L L^T. Upper diag. of a isn't touched, and L is
//          stored in the lower diag. part of a. The diag elements
//          are stored in p, and set in a. 
// ARGS   : 
// RETURNS:  T result
//************************************************************************************
template<class T> void GTMatrix<T>::choldc(T **a, T p[], GINT n)
{
        GINT     i, j, k;
        T        sum;

        for (i=0;i<n;i++) {
          for ( j=i; j<n;j++ ) {
            for ( sum=a[i][j], k=i-1; k>=0; k--) sum -= a[i][k]*a[j][k];
            if ( i == j ) {
              if ( sum <= 0.0 ) { cout << "GTMatrix<T>::choldc (1): failed" << endl; exit(1); }
              p[i] = sqrt((GQUAD)sum);
            } else a[j][i] = sum / (p[i]+GTINY);
          }
        }

        // By default, add the diagonal back into matrix:
        for (i=0;i<n;i++) a[i][i] = p[i];

} // end of method choldc (1)



//************************************************************************************
//************************************************************************************
// METHOD   : choldc (2)
// DESC     : Taken from Numerical Recipes: serves as method to perform
//            efficient Cholesky decomposition of an input matrix, that is
//            symmetric positive definite. That is, 
//            B (=this)  = L L^T.  The result is returned in the 
//            lower triangular part of A; the diagonal part is returned in p, and
//            this is set in A. The upper diagonal part of A is zeroed.
// ARGS     : 
// RETURNS  : T result
//************************************************************************************
template<class T> void GTMatrix<T>::choldc(GTMatrix<T> &a, GTVector<T> &p)
{
        GINT     i, j, k, n;
        GDOUBLE  sum;

        a = 0.0;
        n = this->dim(1);
        for (i=0;i<n;i++) {
          for ( j=i;j<n;j++ ) {
            for ( sum=(*this)(i,j), k=i-1; k>=0; k--) sum -= (*this)(i,k) * (*this)(j,k); 
            if ( i == j ) {
              if ( sum <= 0.0 ) { 
                cout << "GTMatrix<T>::choldc (2): failed" << endl;  
                cout << "a=" << a << endl; 
                exit(1); }
              p[i] = sqrt((GQUAD)sum);
            } else a(j,i) = sum / (p[i]+GTINY);
          }
        }

        // By default, add the diagonal back into matrix:
        for (i=0;i<n;i++) a(i,i) = p[i];

} // end of method choldc (2)


//************************************************************************************
//************************************************************************************
// METHOD : CreateIdentity
// DESC   : Creates identity from the matrix
// ARGS   : none.
// RETURNS: none.
//************************************************************************************
template<class T> void GTMatrix<T>::CreateIdentity()
{
  GINT     i, j, k, n;
  char    *serr = "GTMatrix<T>::CreateIdentity: ";
  
  if ( n1_ != n2_ ) {
    cout << serr << "matrix not square" << endl;
     exit(1);
  }
  (*this) = 0.0;
  for ( GINT i=0; i<n1_; i++ ) (*this)(i,i) = 1;

} // end of method CreateIdentity

