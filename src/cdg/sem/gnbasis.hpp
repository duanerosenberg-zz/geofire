//==================================================================================
// Module       : gnbasis.hpp
// Date         : 1/19/18 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a spectral element basis. This is a pure virtual class, 
//                intended to form the interface to a general nodal basis object.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(_GNBASIS_HPP)
#define _GNBASIS_HPP
#include "gtypes.h"
#include "gtmatrix.hpp"
#include "gtvector.hpp"


class GNBasis
{
public:

                                    GNBasis(){};
                                    GNBasis(GINT){};
                                    GNBasis(GINT, GINT){};
                                   ~GNBasis(){};

//virtual         void              operator=(const GNBasis &)=0;
                   
  virtual         GDOUBLE           GetXimin ()=0;
  virtual         GDOUBLE           GetXimax ()=0; 
  virtual         GINT              GetOrder()=0;
  virtual         GVector          *GetXiNodes()=0;
  virtual         GVector          *GetXiNodes(GVector *)=0;
  virtual         GDOUBLE          *GetXiNodes(GDOUBLE   *, GINT  num)=0;
  virtual         QVector          *GetQXiNodes()=0;
  virtual         GTVector<GQUAD>  *GetQWeights()=0;
  virtual         GVector          *GetWeights()=0;
  virtual         GVector          *GetWeights(GVector *)=0;
  virtual         GDOUBLE          *GetWeights(GDOUBLE   *, GINT  num)=0;

  virtual         GMatrix          *GetStiffMatrix(GBOOL bflag=FALSE)=0;
  virtual         GMatrix          *GetStiffMatrix(GMatrix *)=0;
  virtual         GMatrix          *GetMassMatrix(GMatrix *)=0;
//virtual         GMatrix          *GetBasisAtNodes(GMatrix *)=0;
  virtual         GMatrix          *GetDerivMatrix(GBOOL bflag=FALSE)=0;
  virtual         GMatrix          *GetDerivMatrix(GMatrix *)=0;
//virtual         GDOUBLE          *GetDerivMatrix(GDOUBLE    *, GINT  num)=0;
  virtual         void              SetOrder(GINT  )=0;
  virtual         GDOUBLE           EvalBasis (GINT  i, GDOUBLE xi)=0;
  virtual         GVector          *EvalBasis (GINT  i, GVector *xi, GVector *vret)=0;
  virtual         GMatrix          *EvalBasis (GVector *eta, GMatrix *mret)=0;
  virtual         GMatrix          *EvalBasis (GDOUBLE eta[], GINT neta, GMatrix *mret)=0;
  virtual         GMatrix          *EvalDBasis(GVector *eta, GMatrix *mret)=0;
  virtual         GMatrix          *EvalDBasis(GDOUBLE eta[], GINT neta, GMatrix *mret)=0;
  virtual         GMatrix          *GetLegMatrix(GMatrix *mret)=0;

  virtual         GBOOL             Solve ()=0;


};

#endif
