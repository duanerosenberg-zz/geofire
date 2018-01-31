//==================================================================================
// Module       : gtmatrix.hpp
// Date         : 1/1/18 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a template matrix object, whose data is composed of
//                a regular array of contiguous data, ordered like a Fortran
//                matrix in row-major order (row data changing fastest over column data).
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

template<class T> class GTMatrix;
template<class T> std::ostream &operator<<(std::ostream &, const GTMatrx<T> &);

#if !defined(_GTMATRIX_HPP)
#define _GTMATRIX_HPP
#include <stdlib.h>
#include <memory.h>
#include <iomanip.h>
#include <iostream.h>
#include "gtypes.h"
#include "gtvector.hpp"


#if !defined(GTMATRIX_ROTATE)
#define GTMATRIX_ROTATE(a,i,j,k,l) g=a(i,j);h=a(k,l);a(i,j)=g-s*(h+g*tau);\
        a(k,l)=h+s*(g-h*tau);
#endif

template<class T> class GTMatrix
{
public:
                          GTMatrix();
                          GTMatrix(const GINT , const GINT );
                          GTMatrix(const GINT);
                          GTMatrix(T *array, GINT  n1, GINT  n2);
                          GTMatrix(GTMatrix<T> *m, GINT *ind, GINT  nn);
                          GTMatrix(const GTMatrix<T> &);

                         ~GTMatrix();

         GTMatrix<T>     &operator=(const GTMatrix<T> &);
         void             operator=(T);
         void             operator+=(const GTMatrix<T> &);
         void             operator+=(T);
         GTMatrix<T>     &operator*(GTMatrix<T> &)      ;
         GTVector<T>     &operator*(GTVector<T> & )     ; // right multiplication: Matrix * Array
         void             operator*=(T);
         GTMatrix<T>     &operator*(T)                  ; // right multiplication: Matrix * const
         GTMatrix<T>      operator*=(T)                 ; // multiplication by const Mat *= const

         GTMatrix<T>     &operator+(GTMatrix<T>  ) ; // addition
         GTMatrix<T>     &operator-(GTMatrix<T>  ) ; // addition
         GBOOL            Transpose(GTMatrix<T> &);
         GBOOL            Transpose(GTMatrix<T> &, GINT Nx, GINT ny);
         GTMatrix<T>      Transpose();
         GBOOL            Inverse(GTMatrix<T> &);
         GBOOL            Inverse(GTMatrix<T> &, GINT nx, GINT ny);
         GTMatrix<T>      Inverse();
         GBOOL            isSymmetric();
         void             SetSingularZero(T tol);

         void updatehost() {
           #pragma acc update self( data_[0:n_-1] )
         };
         void updatedev() {
           #pragma acc update device(data_[0:n_-1])
         };


         inline T  &operator()(const GINT  i, const GINT  j ) {
           #if defined(_G_BOUNDS_CHK)
           const char serr[] = "GTMatrix<T>::&(): ";
           if ( i >= n1 || i < 0 || j >= n2 || j < 0  ) {
             cout << serr << "access error: i=" << i << "; j=" << j << endl;
             while(1);
             exit(1);
           }
           #endif
           return data_[i+j*n1]; 
         }; 

         inline T &operator()(const GINT  i){
           #if defined(_G_BOUNDS_CHK)
           const char serr[] = "GTMatrix<T>::&(): ";
           if ( i >= n1*n2 || i < 0 ) {
             cout << serr << "access error: i=" << i << endl;
             while(1);
             exit(1);
           }
           #endif
           return data_[i];
         }; 

         inline T operator()(const GINT  i, const GINT  j ) const {
           #if defined(_G_BOUNDS_CHK)
           const char serr[] = "GTMatrix<T>::() const: ";
           if ( i >= n1 || i < 0 || j >= n2 || j < 0  ) {
             cout << serr << "access error: i=" << i << "; j=" << j << endl;
             exit(1);
           }
           #endif
           return data_[i+j*n1];
         }; 

         inline T operator()(const GINT  i) const {
           #if defined(_G_BOUNDS_CHK)
           const char serr[] = "GTMatrix<T>::() const: ";
           if ( i >= n1*n2 || i < 0 ) {
             cout << serr << "access error: i=" << i << endl;
             exit(1);
           }
           #endif
           return data_[i];
         }; 

         inline GINT dim(GINT  idir) { 
           if      ( idir == 1 ) return n1;
           else if ( idir == 2 ) return n2;
           else                  return 0; 
         };

         inline T *data() { return this->data_.data(); }

         GINT              tsize(GINT  idir);  // _total_ size = n_idir 
         GTVector<T>      &Vdata() ;           // Get vector data container reference
         GBOOL             resize(GINT  Nx, GINT  Ny);
         GBOOL             resizeM(GINT  Nx, GINT  Ny);
         void              Zero();
         void              CreateIdentity();
  
         friend ostream&   operator<< <> (ostream &, GTMatrix<T> & );


void              DeleteDynamic();
GBOOL             ludcmp (T **&a, GINT  n, GINT  *&indx, T *d);
GBOOL             wludcmp(T **&a, GINT  n, GINT  *&indx, T *d);
GBOOL             lubksb (T **&A, GINT  n, GINT  *&indx, T b[]);
GBOOL             improve(T **&a, T **alud, GINT n, GINT *&indx, T b[], T x[]);
GBOOL             svdcmp(T **a, GINT m, GINT n, T w[], T **v);
GBOOL             svdcmp(GTVector<T> &w, GTMatrix<T> &v,  GTVector<T> &rv1);
GBOOL             svdcmp(GTVector<T> &w, GTMatrix<T> &v,  GTVector<T> &rv1, GINT nx, GINT ny);
void              svbksub(T **u, T w[], T **v, GINT n, GINT m, T b[], T x[]);
void              choldc(T **u, T w[], GINT n);
void              choldc(GTMatrix<T> &a, GTVector<T> &p);
void              jacobi(GTVector<T> &d, GTMatrix<T> &v, GINT &nrot, GTMatrix<T> &a, GTVector<T> &b, GTVector<T> &z);


// Private methods:
private:
T                 dpythag(T a, T b);
GINT              isamax(GINT  n, GDOUBLE *sx, GINT  incx);


// Private data:
GINT              n1_;
GINT              n2_;
G_DATATYPE        dtype_;
GC_DATATYPE       cdtype_;

GTVector<T>       data_;
T                 fNULL;
T                 singzero_;
#if defined(DO_BLAS_TIMING)
public:
        GDOUBLE   GetTime(GINT  i){if      ( i==0 ) return time_result;
                                   else if ( i==1 ) return time_result1;
                                   else if ( i==2 ) return time_result2; }
private:
        GDOUBLE   time_result;
        GDOUBLE   time_result1;
        GDOUBLE   time_result2;
#endif


};

#if 0 //defined(_LINUX) || defined(_AIX) || defined(_MACOSX)
template class GTMatrix<GDOUBLE>;
ostream &operator <<(ostream&, GTMatrix<GDOUBLE>&);
template class GTMatrix<GINT>;
ostream &operator <<(ostream&, GTMatrix<GINT>&);
#if  defined(GBUFF_DEF_GQUAD)
template class GTMatrix<GQUAD>;
ostream &operator <<(ostream&, GTMatrix<GQUAD>&);
#    endif
#endif

#endif
