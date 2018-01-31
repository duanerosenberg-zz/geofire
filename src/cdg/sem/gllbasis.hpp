//==================================================================================
// Module       : gllbasis.hpp
// Date         : 1/19/18 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a template vector object, with contiguous memory
//                for basic types.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : GNBasis
//==================================================================================

#if !defined(_GLLBASIS_HPP)
#define _GLLBASIS_HPP
#include "gtypes.h"
#include "gnbasis.hpp"
#include "gtvector.hpp"
#include "gtmatrix.hpp"

class GLLBasis: public GNBasis
{
public:

                           GLLBasis();
                           GLLBasis(GINT  inorder);
                           GLLBasis(GINT  inorder, GINT  maxorder);
                           GLLBasis(const GLLBasis &);
virtual                   ~GLLBasis();

virtual  void              operator=(const GLLBasis &);
                   
         GDOUBLE           GetXimin ();
         GDOUBLE           GetXimax (); 
         GINT              GetOrder();
         GVector          *GetXiNodes();
         GVector          *GetXiNodes(GVector &);
         GDOUBLE          *GetXiNodes(GDOUBLE   *, GINT  num);
         GTVector<GQUAD>  *GetQXiNodes();
         GVector          *GetWeights();
         GVector          *GetWeights(GVector &);
         GDOUBLE          *GetWeights(GDOUBLE   *, GINT  num);
         GTVector<GQUAD>  *GetQWeights();
         GMatrix          *GetStiffMatrix();
         GMatrix          *GetStiffMatrix(GMatrix &);
//       GMatrix          *GetBasisAtNodes(GMatrix *);
         GMatrix          *GetDerivMatrix(GBOOL bTranspose=FALSE);
         GMatrix          *GetDerivMatrix(GMatrix &);
//       GDOUBLE          *GetDerivMatrix(GDOUBLE*, GINT  num);
         GMatrix          *GetLegMatrix(GMatrix &);
         void              SetOrder(GINT );
virtual  GDOUBLE           EvalBasis (GINT  i, GDOUBLE xi);
virtual  GVector          *EvalBasis (GINT  i, GVector &xi, GVector &vret);
virtual  GMatrix          *EvalBasis (GVector &eta, GMatrix &mret);
virtual  GMatrix          *EvalBasis (GDOUBLE eta[], GINT neta, GMatrix &mret);
virtual  GMatrix          *EvalDBasis(GVector &eta, GMatrix &mret);
virtual  GMatrix          *EvalDBasis(GDOUBLE eta[], GINT neta, GMatrix &mret);
         GBOOL             Solve ();

//       GMatrix          *GetDBasisAtNodes(GMatrix *);
         GMatrix          *GetMassMatrix(GMatrix &);
//       void              SetXiDomain(GQUAD min, GQUAD max);

// Protected methods:
protected:

virtual  GBOOL             ComputeNodes       ();
virtual  GBOOL             ComputeWeights     ();
virtual  GBOOL             ComputeDerivMatrix ();
virtual  GBOOL             ComputeLegendreMatrix ();
//virtual  GBOOL           ComputeBasisAtNodes();
         void              ComputeJacobi(GINT  &,GQUAD  alpha, GQUAD  beta, GQUAD &Pn, 
                                         GQUAD &dPn, GQUAD &Pnm1, GQUAD &dPnm1, GQUAD &Pnm2,
                                         GQUAD &dPnm2, GQUAD &xi);
virtual  GBOOL             ComputeMassMatrix();
virtual  GBOOL             ComputeStiffMatrix();

virtual  GBOOL             Resize(GINT  order);


// Protected data:
GINT             Np;
GINT             NpMax;
GINT             kstop;
GBOOL            bNeedNodes;
GBOOL            bNeedWeights;
GBOOL            bNeedDerivMatrix;
GBOOL            bNeedBasis;
GBOOL            bNeedDBasis;
GBOOL            bNeedLegMat;
GQUAD            alpha;
GQUAD            beta;
GQUAD            ximin;
GQUAD            ximax;
GQUAD            eps;
GVector          xid;
GVector          weightsd;
GBasisVector     xiNodes;
GBasisVector     Weights;
GBasisVector     Pn;
GBasisVector     dPn;
GBasisMatrix     GXi;
GBasisMatrix     Phi;
GBasisMatrix     dPhi;
GMatrix          dPhid;
GMatrix          dPhidT;
GMatrix          MassMatrixd;
GBasisMatrix     StiffMatrix;
GMatrix          StiffMatrixd;
GBasisMatrix     LegMatrix;

};

#endif
