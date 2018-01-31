//==================================================================================
// Module       : glbasis.hpp
// Date         : 1/19/18 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a spectral element nodal Gauss-Legendre basis
//                Note: This basis is the Gauss-Legendre
//                basis, defined as
//                  h(xi)_j = L_N(xi) / ( dL_N(xi_j)/dxi * (xi-xi_j) )
//                where N is the order of the expansion, and L_N is the Legendre
//                polynomial of order N, and the xi_j are the nodal points. L_N
//                is obtained from the generalized Jacobi polynomial computed in
//                ComputeJacobi.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : GLLBasis
//==================================================================================
#if !defined(_GBASIS_HPP)
#define _GBASIS_HPP
#include "gtypes.h"
#include "gllbasis.hpp"


class GLBasis: public GLLBasis
{
public:

                           GLBasis();
                           GLBasis(GINT );
                           GLBasis(GINT  , GINT );
                           GLBasis(const GLBasis &);
virtual                   ~GLBasis();

         void              operator=(const GLBasis &);
                   
         GDOUBLE           EvalBasis (GINT  i, GDOUBLE xi);
         GVector          *EvalBasis (GINT  i, GVector *xi, GVector &vret);
         GMatrix          *EvalBasis (GVector &eta, GMatrix &mret);
         GMatrix          *EvalBasis (GDOUBLE eta[], GINT neta, GMatrix &mret);
         GMatrix          *EvalDBasis(GVector &eta, GMatrix &mret);
         GMatrix          *EvalDBasis(GDOUBLE eta[], GINT neta, GMatrix &mret);


// Private methods:
private:

         GBOOL             ComputeNodes       ();
         GBOOL             ComputeWeights     ();
         GBOOL             ComputeDerivMatrix ();
         GBOOL             ComputeBasisAtNodes();
         GBOOL             ComputeLegendreMatrix();

         GBOOL             Resize(GINT  order);
//       void              DeleteDynamic();

};

#endif
