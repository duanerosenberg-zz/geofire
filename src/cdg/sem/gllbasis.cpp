//==================================================================================
// Module       : gllbasis.cpp
// Date         : 1/19/18 (DLR)
//                Research
// Description  : Encapsulates the methods and data associated with
//                a spectral element nodal basis
//                Note: This basis is the Gauss-Lobatto-Legendre
//                basis, defined as 
//                  h(xi)_j = -1/(N(N+1)*L_N(xi_j))  * (1-xi^2)*dL_N(xi)/dxi / (xi-xi_j)
//                where N is the order of the expansion, and L_N is the Legendre
//                polynomial of order N, and the xi_j are the nodal points. L_N
//                is obtained from the generalized Jacobi polynomial computed in
//                ComputeJacobi.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : GNBasis
//==================================================================================
#include "gllbasis.hpp"
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <stdio.h>


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method
// DESC   : Instantitate with explicit order and max order
// ARGS   : GINT inOrder
//          GINT maxORdern
// RETURNS: none
//**********************************************************************************
GLLBasis::GLLBasis(GINT  inOrder, GINT  MaxOrder)
:
Np                    (MIN(inOrder,MaxOrder)),
NpMax                 (MAX(inOrder,MaxOrder)),
kstop                 (128),
bNeedNodes            (TRUE),
bNeedWeights          (TRUE),
bNeedDerivMatrix      (TRUE),
bNeedBasis            (TRUE),
bNeedDBasis           (TRUE),
bNeedLegMat           (TRUE),
alpha                 (0.0),
beta                  (0.0),
ximin                 (-1.0),
ximax                 (1.0),
eps                   (1.0e-16)
{
  if ( Np < 1 ) {
    cout << "GLLBasis::GLLBasis: invalid expansion order Np=" << Np << endl;
    exit(1);
  }
  xid         .Resize(1);
  weightsd    .Resize(1);
  xiNodes     .Resize(1);
  Weights     .Resize(1);
  Pn          .Resize(1);
  dPn         .Resize(1);
  Phi         .Resize(1,1);
  dPhi        .Resize(1,1);
  dPhid       .Resize(1,1);
  dPhidT      .Resize(1,1);
  MassMatrixd .Resize(1,1);
  StiffMatrix .Resize(1,1);
  StiffMatrixd.Resize(1,1);
  LegMatrix   .Resize(1,1);

  Resize(NpMax);
} // end of constructor method


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method
// DESC   : Default constructor
// ARGS   : GINT inOrder
//          GINT maxORdern
// RETURNS: none
//**********************************************************************************
GLLBasis::GLLBasis()
:
Np                    (4),
NpMax                 (32),
kstop                 (128),
bNeedNodes            (TRUE),
bNeedWeights          (TRUE),
bNeedDerivMatrix      (TRUE),
bNeedBasis            (TRUE),
bNeedDBasis           (TRUE),
bNeedLegMat           (TRUE),
alpha                 (0.0),
beta                  (0.0),
ximin                 (-1.0),
ximax                 (1.0),
eps                   (1.0e-16)
{
  if ( Np < 1 ) {
    cout << "GLLBasis::GLLBasis: invalid expansion order Np=" << Np << endl;
    exit(1);
  }
  xid         .Resize(1);
  weightsd    .Resize(1);
  xiNodes     .Resize(1);
  Weights     .Resize(1);
  Pn          .Resize(1);
  dPn         .Resize(1);
  Phi         .Resize(1,1);
  dPhi        .Resize(1,1);
  dPhid       .Resize(1,1);
  dPhidT      .Resize(1,1);
  MassMatrixd .Resize(1,1);
  StiffMatrix .Resize(1,1);
  StiffMatrixd.Resize(1,1);
  LegMatrix   .Resize(1,1);
  Resize(NpMax);
} // end of constructor method


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method
// DESC   : Instantiate with explicit order, taken to be max order as well
// ARGS   : GINT inOrder
// RETURNS: none
//**********************************************************************************
GLLBasis::GLLBasis(GINT Order)
:
Np                    (Order),
NpMax                 (Order),
kstop                 (128),
bNeedNodes            (TRUE),
bNeedWeights          (TRUE),
bNeedDerivMatrix      (TRUE),
bNeedBasis            (TRUE),
bNeedDBasis           (TRUE),
bNeedLegMat           (TRUE),
alpha                 (0.0),
beta                  (0.0),
ximin                 (-1.0),
ximax                 (1.0),
eps                   (1.0e-16)
{
  if ( Np < 1 ) {
    cout << "GLLBasis::GLLBasis: invalid expansion order Np=" << Np << endl;
    exit(1);
  }
  xid         .Resize(1);
  weightsd    .Resize(1);
  xiNodes     .Resize(1);
  Weights     .Resize(1);
  Pn          .Resize(1);
  dPn         .Resize(1);
  Phi         .Resize(1,1);
  dPhi        .Resize(1,1);
  dPhid       .Resize(1,1);
  dPhidT      .Resize(1,1);
  MassMatrixd .Resize(1,1);
  StiffMatrix .Resize(1,1);
  StiffMatrixd.Resize(1,1);
  LegMatrix   .Resize(1,1);
  Resize(NpMax);

} // end of constructor method


//**********************************************************************************
//**********************************************************************************
// METHOD : Copy Constructor method
// DESC   : 
// ARGS   : GLLBasis &
// RETURNS: none
//**********************************************************************************
GLLBasis::GLLBasis(const GLLBasis &b)
{
   // copy data:
    alpha    = b.alpha;
    beta     = b.beta;
    ximin    = b.ximin;
    ximax    = b.ximax;
    eps      = b.eps;
    Np       = b.Np;
    NpMax    = b.NpMax;
    kstop    = b.kstop;
    bNeedNodes = b.bNeedNodes;
    bNeedWeights = b.bNeedWeights;
    bNeedDerivMatrix = b.bNeedDerivMatrix;
    bNeedBasis   = b.bNeedBasis; 
    bNeedDBasis  = b.bNeedDBasis;
    bNeedLegMat  = b.bNeedLegMat;


    //  matrices, and node data:
    xid         = b.xid;
    weightsd    = b.weightsd;
    xiNodes     = b.xiNodes;
    Weights     = b.Weights;
    Pn          = b.Pn;
    dPn         = b.dPn;
    Phi         = b.Phi;
    dPhi        = b.dPhi;
    dPhid       = b.dPhid;
    dPhidT      = b.dPhidT;
    MassMatrixd = b.MassMatrixd;
    StiffMatrix = b.StiffMatrix;
    StiffMatrixd= b.StiffMatrixd;
    LegMatrix   = b.LegMatrix;

}

//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
GLLBasis::~GLLBasis()
{
}

//**********************************************************************************
//**********************************************************************************
// METHOD : Assignment method
// DESC   : 
// ARGS   : GLLBasis &
// RETURNS: none
//**********************************************************************************
void GLLBasis::operator=(const GLLBasis &b)
{
  if ( &b != this ) 
  {
   // copy data:
    ximin        = b.ximin;
    ximax        = b.ximax;
    eps          = b.eps;
    Np           = b.Np;
    NpMax        = b.NpMax;
    kstop        = b.kstop;
    bNeedNodes   = b.bNeedNodes;
    bNeedWeights = b.bNeedWeights;
    bNeedDerivMatrix = b.bNeedDerivMatrix;
    bNeedBasis   = b.bNeedBasis; 
    bNeedDBasis  = b.bNeedDBasis;
    bNeedLegMat  = b.bNeedLegMat;


    //  matrices, and node data:
    xid         = b.xid;
    weightsd    = b.weightsd;
    xiNodes     = b.xiNodes;
    Weights     = b.Weights;
    Pn          = b.Pn;
    dPn         = b.dPn;
    Phi         = b.Phi;
    dPhi        = b.dPhi;
    dPhid       = b.dPhid;
    dPhidT      = b.dPhidT;
    MassMatrixd = b.MassMatrixd;
    StiffMatrix = b.StiffMatrix;
    StiffMatrixd= b.StiffMatrixd;
    LegMatrix   = b.LegMatrix;
  }

}

//************************************************************************************
//************************************************************************************
// METHOD : GetXimin
// DESC   : Get minimim of reference interval 
// ARGS   : none
// RETURNS: GQUAD element left boundary
//************************************************************************************
GDOUBLE GLLBasis::GetXimin()
{
  return (GDOUBLE) ximin;
} // end of method GetXimin


//************************************************************************************
//************************************************************************************
// METHOD : GetXimax
// DESC   : Get maximum of reference interval 
// ARG    : none
// RETURNS: GQUAD element right boundary
//************************************************************************************
GDOUBLE GLLBasis::GetXimax()
{
  return (GDOUBLE) ximax;
} // end of method GetXimax



//************************************************************************************
//************************************************************************************
// METHOD : GetOrder
// DESC   : Get Lag. interp. polynomial order  
// ARG    : none
// RETURNS: GINT  element expansion order
//************************************************************************************
GINT  GLLBasis::GetOrder()
{
  return Np;
} // end of method GetOrder


//************************************************************************************
//************************************************************************************
// METHOD : GetXiNodes
// DESC   : Get list (GDOUBLE vector) of reference interval nodes
// ARGS   : none
// RETURNS: GVector*
//************************************************************************************
GVector *GLLBasis::GetXiNodes()
{
  return &xid;
} // end of method GetXiNodes


//************************************************************************************
//************************************************************************************
// METHOD : GetXiNodes (2)
// DESC   : Get list (GDOUBLE array *) or reference interval nodes
// ARGS   : ret : array of GDOUBLE
//          num : num elements in ret array
// RETURNS: GQUAD *
//************************************************************************************
GDOUBLE *GLLBasis::GetXiNodes(GDOUBLE *ret, GINT  num)
{

  if ( ret == NULLPTR || num < xiNodes.dim() ) return NULLPTR;
  if ( bNeedNodes || bNeedWeights )
    if ( !ComputeNodes() ) return NULLPTR;
  GINT  i;
  GQUAD *qptr=xiNodes.data();


  for ( i=0; i<xiNodes.dim(); i++ )
    ret[i] = (GDOUBLE)(*(qptr+i));

  return ret;
} // end of method GetXiNodes (2)


//************************************************************************************
//************************************************************************************
//************************************************************************************
// METHOD : GetXiNodes (3)
// DESC   : Get list of GDOUBLE reference nodes, returning in specified array GTVector
// ARGS   : GTVector &ret
// RETURNS: pointer to input vector on success; else NULLPTR
//************************************************************************************
GVector *GLLBasis::GetXiNodes(GVector &ret)
{
  GetXiNodes(ret.data(), ret.dim());
  return &ret;
} // end of method GetXiNodes (3)


//************************************************************************************
//************************************************************************************
// METHOD : GetQ1XiNodes 
// DESC   : GetGQUAD nodes
// ARGS   : none
// RETURNS: GTVector<GQUAD> *
//************************************************************************************
GTVector<GQUAD> *GLLBasis::GetQXiNodes()
{
  return &xiNodes;
} // end of method GeaQXiNodes 



//************************************************************************************
//************************************************************************************
// METHOD : GetWeights
// DESC   : Get GTVector<GDOUBLE> member data vector weighters 
// ARGS   : none
// RETURNS: pointer to GVector member data
//************************************************************************************
GVector *GLLBasis::GetWeights()
{
  return &weightsd;
} // end of method GetWeights 


//************************************************************************************
//************************************************************************************
// METHOD : GetWeights (2)
// DESC   : Get GDOUBLE weights in specified return array
// ARGS   : ret : GDOUBLE array
//          num : num elements in ret array
// RETURNS: GQUAD *
//************************************************************************************
GDOUBLE *GLLBasis::GetWeights(GDOUBLE *ret, GINT  num)
{

  if ( ret == NULLPTR || num < Weights.dim() ) return NULLPTR;
  if ( bNeedNodes || bNeedWeights )
    if ( !ComputeNodes() ) return NULLPTR;

  GINT  i;
  GQUAD *qptr=Weights.data();

  for ( i=0; i<Weights.dim(); i++ )
    *(ret+i) = (GDOUBLE) (*(qptr+i));

  return ret;
} // end of method GetWeights (2)


//************************************************************************************
//************************************************************************************
// METHOD : GetWeights (3)
// DESC   : Get GTVector<GDOUBLE> weights in specified return object
// ARGS   : ret : GVector<GDOUBLE> array
// RETURNS: pointer to ret GTVector
//************************************************************************************
GVector *GLLBasis::GetWeights(GVector &ret)
{

  GetWeights(ret.data(), ret.dim());

  return &ret;
} // end of method GetWeights (3)


//************************************************************************************
//************************************************************************************
// METHOD : GetQWeights 
// DESC   : Get GTVector<GDOUBLE> weights in specified return object
// ARGS   : none
// RETURNS: pointer to member GTVector<GQUAD>
//************************************************************************************
GTVector<GQUAD> *GLLBasis::GetQWeights()
{
  return &Weights;
} // end of method GetDWeights 


#if 1
//************************************************************************************
//************************************************************************************
// METHOD : GetMassMatrix
// DESC   : Get 1d MassMatrix in specified return array
// ARGS   : ret: GTMatrix<GDOUBLE> return matrix
// RETURNS: pointer to ret GMatrix
//************************************************************************************
GMatrix *GLLBasis::GetMassMatrix(GMatrix &ret)
{
  if ( ret.dim(1) < MassMatrixd.dim(1) || ret.dim(2) < MassMatrixd.dim(2) ) return NULLPTR;

  GINT  i, j;

  for ( i=0; i<ret->dim(1); i++ ) {
    for ( j=0; j<ret->dim(2); j++ ) {
      ret(i,j) = MassMatrixd(i,j);
    }
  }

  return &ret;
} // end of method GetMassMatrix
#endif


//************************************************************************************
//************************************************************************************
// METHOD : GetStiffMatrix
// DESC   : Get pointer to member data 1d stiffness GTMatrix<GDOUBLE> 
// ARGS   : none
// RETURNS: GMatrix *
//************************************************************************************
GMatrix *GLLBasis::GetStiffMatrix()
{
  return &StiffMatrixd;
} // end of method GetStiffMatrix


//************************************************************************************
//************************************************************************************
// METHOD : GetStiffMatrix
// DESC   : Copies 1d stiffness matrix to to supplied object
// ARGS   : ret: GTMatrix<GDOUBLE> 
// RETURNS: pointer to return matrix on success; else NULLPTR
//************************************************************************************
GMatrix *GLLBasis::GetStiffMatrix(GMatrix &ret)
{
  if ( ret.dim(1) < StiffMatrix.dim(1) || ret.dim(2) < StiffMatrix.dim(2) ) return NULLPTR;

  GINT  i, j; 

  for ( i=0; i<StiffMatrix.dim(1); i++ ) 
    for ( j=0; j<StiffMatrix.dim(2); j++ ) 
      ret(i,j) = (GDOUBLE)StiffMatrix(i,j);

  return &ret;
} // end of method GetStiffMatrix


//************************************************************************************
//************************************************************************************
// METHOD : GetLegMatrix
// DESC   : Gets/fills Legendre poly. matrix: Lp = Lp(nLegOrder,iNodeIndex)
// ARGS   : ret: GTMatrix<GDOUBLE> 
// RETURNS: pointer to return GMatrix on success; else NULPTR
//************************************************************************************
GMatrix *GLLBasis::GetLegMatrix(GMatrix &ret)
{
  if ( ret.dim(1) < LegMatrix.dim(1) || ret.dim(2) < LegMatrix.dim(2) ) return NULLPTR;
  if ( bNeedLegMat && !ComputeLegendreMatrix() ) return NULLPTR;

  GINT  i, j;

  for ( i=0; i<LegMatrix.dim(1); i++ )
    for ( j=0; j<LegMatrix.dim(2); j++ )
      ret(i,j) = (GDOUBLE)LegMatrix(i,j);

  return &ret;
} // end of method GetLegMatrix


//************************************************************************************
//************************************************************************************
// METHOD : GetDerivMatrix
// DESC   : Get derivative matrix member data
// ARGS   : bTranspose: TRUE==>return transpose; else don't
// RETURNS: pointer to member data GMatrix
//************************************************************************************
GMatrix *GLLBasis::GetDerivMatrix(GBOOL bTranspose)
{
  if ( bTranspose ) return &dPhidT;
  else              return &dPhid;
} // end of method GetDerivMatrix


//************************************************************************************
//************************************************************************************
// METHOD : GetDerivMatrix
// DESC   : Get deriv matrix in specified return object
// ARGS   : ret: GMatrix to return data
// RETURNS: pointer to ret GMatrix on success; else NULLPTR
//************************************************************************************
GMatrix *GLLBasis::GetDerivMatrix(GMatrix &ret)
{
  if ( ret.dim(1) < dPhi.dim(1) || ret.dim(2) < dPhi.dim(2) ) return NULLPTR;

  if ( bNeedNodes )
    if ( !ComputeNodes() ) return NULLPTR;
  if ( bNeedDerivMatrix )
    if ( !ComputeDerivMatrix() ) return NULLPTR;

  GINT  i, j; 

  for ( i=0; i<dPhi.dim(1); i++ ) 
    for ( j=0; j<dPhi.dim(2); j++ ) 
      ret(i,j) = (GDOUBLE)dPhi(i,j);
  
  return &ret;
} // end of method GetDerivMatrix


#if 0
//************************************************************************************
//************************************************************************************
// METHOD : GetDerivMatrix (3)
// DESC   : Get deriv matrix in specified return array
// ARGS   : ret: array o GDOUBLE
//          num: size of ret array
// RETURNS: pointer to ret array on success; else NULLPTR
//************************************************************************************
GDOUBLE *GLLBasis::GetDerivMatrix(GDOUBLE *ret, GINT  num)
{
  GINT  i,j, n=0;

  if ( ret == NULLPTR || num < dPhi.dim(1)*dPhi.dim(2) ) return NULLPTR;
  if ( bNeedNodes )
    if ( !ComputeNodes() ) return NULLPTR;
  if ( bNeedDerivMatrix )
    if ( !ComputeDerivMatrix() ) return NULLPTR;

  for ( i=0; i< dPhi.dim(1); i++ ) {
    for ( j=0; j< dPhi.dim(2); j++ ) {
      ret[n] = (GDOUBLE) dPhi(i,j);
      n++;
    }
  }
  return ret;
} // end of method GetDerivMatrix (3)


//************************************************************************************
//************************************************************************************
// METHOD : GetBasisAtNodes
// DESC   : Computes Gauss-Lobatto Legendre basis
//              as a function of  nodes in the parent domain. 
//              The returned quantity is a pointer to a matrix
//              Phi_i(xi_j)==Phi(i,j),
//              where xi_j is the jth node element, and
//              Phi_i is given by (see Canuto, et. al.):
//
//              Phi_j (x) = (N*(N+1)Pn(xi_j))^-1  * (1-x^2)*dPn/dx(x)/(x - xi_j)
//              
// ARGS   : ret: array to GDOUBLE
// RETURNS: pointer to ret data on success; else ; else NULLPTR 
//************************************************************************************
GBasisMatrix *GLLBasis::GetBasisAtNodes(GBasisMatrix &ret)
{
  if ( bNeedNodes || bNeedWeights )
    if ( !ComputeNodes() ) return NULLPTR;
 
  if ( bNeedBasis       && !ComputeBasisAtNodes() ) return NULLPTR;


  if ( ret.dim(1) < Phi.dim(1) || ret.dim(2) < Phi.dim(2) ) return NULLPTR;

  GINT  i;
  GDOUBLE *fptr=ret.data();
  GQUAD *qptr=Phi.data();
   
  for ( i=0; i<Phi.dim(1) * Phi.dim(2); i++ )
    *(fptr+i) = (GDOUBLE)(*(qprt+i));

  return &ret;
} // end of method GetBasisAtNodes


//************************************************************************************
//************************************************************************************
// METHOD : SetXiDomain
// DESC   : Set reference domain 
// ARGS   : min, max GQUAD 1d domain extent
// RETURNS: none
//************************************************************************************
void GLLBasis::SetXiDomain(GQUAD min, GQUAD max)
{
  ximin = MIN(min,max);
  ximax = MAX(min,max);

} // end of method SetXiDomain
#endif


//************************************************************************************
//************************************************************************************
// METHOD : SetOrder
// DESC   : 
// ARGS   : order: interp. polynomial order
// RETURNS: none
//************************************************************************************
void GLLBasis::SetOrder(GINT  order)
{
  if ( order != Np )
  {
     Resize(order);
     bNeedNodes  = TRUE;
     bNeedWeights = TRUE;
     bNeedBasis  = TRUE;
     bNeedDBasis = TRUE;
     bNeedDerivMatrix = TRUE;
  }
  Np = order;
} // end of method SetOrder


//************************************************************************************
//************************************************************************************
// METHOD : Resize
// DESC   : Resizes dynamically allocated quantities
//          if required
// ARGS   : newOrder: new interp. polynomial order
// RETURNS: TRUE on success, else FALSE
//************************************************************************************
GBOOL GLLBasis::Resize(GINT  newOrder)
{

  // No resizing necessary if already less than
  // previously allocated quantities:

  NpMax = Np = newOrder;

  //  Resize xiNodes:
  xiNodes .Resize(Np+1);
  xid     .Resize(Np+1);
  weightsd.Resize(Np+1);

  //  Resize Weights:
  Weights.Resize(Np+1);

  //  Resize Pn:
  Pn.Resize(Np+1);

  //  Resize dPn:
  dPn.Resize(Np+1);

  //  Resize basis Phi:
  Phi.Resize(Np+1,Np+1);

  //  Resize basis dPhi:
  dPhi  .Resize(Np+1,Np+1);
  dPhid .Resize(Np+1,Np+1);
  dPhidT.Resize(Np+1,Np+1);

  //  Resize MassMatrix:
  MassMatrixd.Resize(Np+1,Np+1);

  //  Resize StiffMatrix:
  StiffMatrix  .Resize(Np+1,Np+1);
  StiffMatrixd .Resize(Np+1,Np+1);

  //  Resize LegMatrix:
  LegMatrix.Resize(Np+1,Np+1);

  return TRUE;

} // end of method Resize


//************************************************************************************
//************************************************************************************
// METHOD : ComputeNodes 
// DESC   : Computes nodes and weights based on a Gauss-Lobatto Legendre integration scheme
//              Note: nodes are computed from smalles to largest xi
//              Note: this method was taken largely from Canuto et al, 1987, Appendix C.
//              Note: this method will be generalized to accept this as a base class,
//                    and derive other basis types from it....
//              Note: Legendre polynomials of order Np, Np-1 and Np-2 and their
//                    derivatives are computed here for use here and in 
//                    computing integrals. These polynomials are used in
//                    forming basis functions, Phi, and their derivatives,
//                    dPhi, evaluated at each node.
// ARGS   : none
// RETURNS: TRUE on success, else FALSE
//************************************************************************************
GBOOL GLLBasis::ComputeNodes()
{
  GINT  i, j, k, nh, np1;
//GQUAD rv;
  GQUAD det, pnp, pnp1p, pnp1m, rp, rm, dth, cd, sd, cs, ss, 
        pn , pnp1, pdnp1, pdn, pnm1, pdnm1, a, b, poly, pder,
        pdnp1p, pnm1p,  pdnp1m, pnm, pdnm, pnm1m, pdnp,
        recsum, x, delx, cssave, error;
  char  *serr = "GLLBasis::ComputeNodes: ";

  if ( !bNeedNodes && !bNeedWeights ) return TRUE;

  if ( Np < 1 ) return FALSE;
  if ( Np == 1 ) {
    xiNodes(0) = ximin;
    xiNodes(1) = ximax;
    bNeedNodes = FALSE;
    for ( i=0; i<Np+1; i++ ) xid[i] = (GDOUBLE)xiNodes[i];
    return(ComputeWeights());
  }
 
  alpha = 0.0;
  beta  = 0.0;  //alpha=beta=0 ==> Legendre fcn basis 

//rv    = 1.0 + alpha;
  np1   = Np + 1;

  ComputeJacobi(np1,alpha, beta, pnp1p, pdnp1p, pnp, pdnp, pnm1p, pdnm1, ximax);
  ComputeJacobi(np1,alpha, beta, pnp1m, pdnp1m, pnm, pdnm, pnm1m, pdnm1, ximin);
  det  = pnp*pnm1m - pnm*pnm1p;
  rp   = -pnp1p;
  rm   = -pnp1m;
  a    = ( rp*pnm1m - rm*pnm1p ) / det;
  b    = ( rm*pnp   - rp*pnm   ) / det;

  // order nodes from largest to smallest:
  xiNodes (0) = ximax;
  xiNodes(Np) = ximin;
  nh = ( Np + 1 ) / 2;

  // set up recursion relation for the initial guesses for roots:
  dth = 3.1415926535898/( 2.0*Np + 1.0 );
  cd  = cos(2.0*dth);
  sd  = sin(2.0*dth);
  cs  = cos    (dth);
  ss  = sin    (dth);

  // compute the first half of the roots by polynomial deflation:
  for ( j=1; j<nh+1; j++ ) {
    x = cs;
    for ( k=0; k<kstop; k++ ) {
      ComputeJacobi(np1, alpha, beta, pnp1, pdnp1, pn, pdn, pnm1, pdnm1, x);
      poly = pnp1 + a*pn + b*pnm1;
      pder = pdnp1 + a*pdn + b*pdnm1;
      recsum = 0.0;
      for ( i=0; i<j; i++ ) {
        recsum += (  1.0/(x - xiNodes(i)) );
      }
      delx = -poly/(pder - recsum*poly);
      x += delx;
      error = fabs(delx)/ fabs(x);    // MIN(fabs(ximin),fabs(ximax));
      if ( error < eps ) break;
    }
    xiNodes(j) = x;
    cssave     = cs*cd - ss*sd;
    ss         = cs*sd + ss*cd;
    cs         = cssave;
  }

  
  // Symmetrize the nodes:  NOTE: this is correct only for the cases
  // where (xmin,xmax) = (-y, y); ie, is symmetric interval.
  for ( i=0; i<nh; i++ ) {
    xiNodes(np1-i-1) = -xiNodes(i);
  }

  bNeedNodes   = FALSE;
  bNeedWeights = TRUE;

  if ( Np == ( 2*(Np/2) ) ) {
    xiNodes(nh) = 0.0;
  }


  // re-order from smallest to largest:
  GQUAD txi;
  for ( i=0; i<nh; i++ ) {
    txi           = xiNodes(Np-i);
    xiNodes(Np-i) = xiNodes(i);
    xiNodes   (i) = txi;
  }

  for ( i=0; i<Np+1; i++ ) xid[i] = (GDOUBLE)xiNodes[i];

  return(ComputeWeights());

} // end of method ComputeNodes


//************************************************************************************
//************************************************************************************
// METHOD : ComputeWeights
// DESC   : 
//            NOTE: method not really intended to be called publicly; 
//            it should be called only when ComputeNodes is called. 
//            For other bases, this method will change, however.
// ARGS   : none
// RETURNS: TRUE on success, else FALSE
//************************************************************************************
GBOOL GLLBasis::ComputeWeights()
{
  if ( Np < 1 ) return FALSE;
 
  GINT  i;
  GQUAD fact = 2.0/(Np*(Np + 1.0));
  GQUAD ppn, pder, pm1, pdm1, pm2, pdm2;

  for ( i=0; i<Np+1; i++ ) {
    ComputeJacobi(Np, alpha, beta, ppn, pder,pm1, pdm1, pm2, pdm2, xiNodes(i));
    Pn (i) = ppn;    
    dPn(i) = pder;    
    Weights(i) = (Pn(i)==0.0)?0.0:fact/(Pn(i)*Pn(i)); // Note: Pn computed in ComputNodes
  }
  for ( i=0; i<Np+1; i++ ) {
    weightsd    [i] = (GDOUBLE)Weights[i];
  }

  ComputeMassMatrix();

  bNeedWeights = FALSE;
  return TRUE;
} // end of method ComputeWeights


#if 0
//************************************************************************************
//************************************************************************************
// METHOD : ComputeBasisAtNodes
// DESC   :
// ARGS   : none
// RETURNS: TRUE on success, else FALSE
//************************************************************************************
GBOOL GLLBasis::ComputeBasisAtNodes()
{
  if ( bNeedNodes || bNeedWeights )
    if ( !ComputeNodes() ) return NULLPTR;


  GINT  i, j;
  GQUAD ppn_i, pder_i, pm1, pdm1, pm2, pdm2, ppn_j, pder_j;
  GQUAD fact=1.0/(Np*(Np+1.0)), gfact;

  Phi = 0.0;

  for ( i=0; i< Np+1; i++ ) {
    ComputeJacobi(Np, alpha, beta, ppn_i, pder_i,pm1, pdm1, pm2, pdm2, xiNodes(i));
    for ( j=0; j< Np+1; j++ ) {
      Phi(i,j) = 1.0;
      if ( i != j ) {
        ComputeJacobi(Np, alpha, beta, ppn_j, pder_j,pm1, pdm1, pm2, pdm2, xiNodes(j));
        gfact = (1-xiNodes(j)*xiNodes(j))/(xiNodes(j)-xiNodes(i));
        Phi(i,j) = fact/ppn_i * gfact * pder_j;
      }
    }
  }
  bNeedBasis = FALSE;

} // end of method ComputeBasisAtNodes
#endif

//************************************************************************************
//************************************************************************************
// METHOD : ComputeDerivMatrix
// DESC   :
//            NOTE: ComputeWeights (ComputeNodes) must have been called
//                  prior to entering this method, s.t. the Pn_i have been
//                  calculated. 
// ARGS   : none
// RETURNS: TRUE on success, else FALSE
//************************************************************************************
GBOOL GLLBasis::ComputeDerivMatrix()
{
  if ( bNeedWeights ) return FALSE;

  GINT  l, j;
  GQUAD delxi;

  dPhi = 0.0;
  // Note: index j sweeps over basis number; l sweeps over node number

  for ( j=0; j<Np+1; j++ ) {
    for ( l=0; l<j; l++ ) {
      delxi      = xiNodes(l) - xiNodes(j);
      dPhi (l,j) =  Pn(l)/(Pn(j)*delxi+GTINY);
    }
    for ( l=j+1; l<Np+1; l++ ) {
      delxi       = xiNodes(l) - xiNodes(j);
      dPhi  (l,j) =  Pn(l)/(Pn(j)*delxi+GTINY);
    }
  }
  dPhi (0 ,0 ) = -0.25*Np*(Np + 1.0);
  dPhi (Np,Np) = -dPhi(0,0);

  for ( j=0; j<Np+1; j++ )
    for ( l=0; l<Np+1; l++ )
      dPhid(j,l) = dPhi(j,l);

  dPhid.Transpose(dPhidT);
  bNeedDerivMatrix = FALSE;
  return TRUE;

} // end of method ComputeDerivMatrix


//************************************************************************************
//************************************************************************************
// METHOD : ComputeJacobi
// DESC   : Compute Jacobi polynomial nodes and derivatives for polynomial
//          type specified by alpha, beta. Taken from: 
// ARGS   : 
// RETURNS: none  
//************************************************************************************
void GLLBasis::ComputeJacobi(GINT  &N, GQUAD  alpha , GQUAD  beta  , 
                                       GQUAD &poly  , GQUAD &pder  , GQUAD &polym1, 
                                       GQUAD &pderm1, GQUAD &polym2, GQUAD &pderm2, GQUAD &x)
{

  GINT  k; 
  GQUAD apb  , polylist, pderlist, rv                ;
  GQUAD polyn, pdern   , psave   , pdsave  , fk      ;
  GQUAD a1   , a2      , a3      , a4      , b3      ;


  apb = alpha + beta;
  rv = 1.0 + alpha;

  poly = 1.0;
  pder = 0.0;

  if ( N == 0 ) return;

  polylist = poly;
  pderlist = pder;
  poly     = rv * x;
  pder     = rv;
  
  if ( N == 1 ) return;

  for ( k=2; k<=N; k++ ) {
    fk = double(k);
    a1 = 2.0*fk*(fk+apb)*(2.0*fk+apb-2.0);
    a2 = (2.0*fk+apb-1.0)*(alpha*alpha -beta*beta);
    b3 = 2.0*fk+apb-2.0;
    a3 = b3*(b3+1.0)*(b3+2.0);
    a4 = 2.0*(fk+alpha-1.0)*(fk+beta-1.0)*(2.0*fk+apb);
    polyn    = ((a2+a3*x)*poly - a4*polylist)/a1;
    pdern    = ((a2+a3*x)*pder - a4*pderlist+a3*poly)/a1;
    psave    = polylist;
    pdsave   = pderlist;
    polylist = poly;
    poly     = polyn;
    pderlist = pder;
    pder     = pdern;
  }

  polym1 = polylist;
  pderm1 = pderlist;
  polym2 = psave;
  pderm2 = pdsave; 

} // end of method ComputeJacobi


//************************************************************************************
//************************************************************************************
// METHOD : ComputeMassMatrix
// DESC   : Computes mass matrix for this basis.
//          Method uses Gaussian integration, so the
//          nodes and the weights for the current Np order
//          must be calculated.
//
//          The mass matrix is defined as:
//              
//          M_(i,j )= Integral(ximin,ximax) { Phi_i(xi) Phi_j(xi) dxi }
//                  ~ Sum { w_k * Phi_i (xi_k) Phi_j(xi_k) }   (Gaussian quadrature)
//          where w_k are the weights for the basis associated with each node, xi_k.
// ARGS   : none 
// RETURNS: TRUE on success, else FALSE
//************************************************************************************
GBOOL GLLBasis::ComputeMassMatrix()
{
  GINT  i;

  MassMatrixd = 0.0;
 
  for ( i=0; i<Np+1; i++ ) {
    MassMatrixd(i,i) = (GDOUBLE)Weights(i);
  }

  return TRUE;
} // end of method ComputeMassMatrix
 

//************************************************************************************
//************************************************************************************
// METHOD : ComputeStiffMatrix
// DESC   : Computes stiffness matrix for this basis.
//          Method uses Gaussian integration, so the
//          nodes and the weights for the current Np order
//          must be calculated.
//
//          The stiffness matrix is defined as:
//              
//          A_(i,j )= Integral(ximin,ximax) { dPhi_i(xi)/dxi dPhi_j(xi)/dxi dxi }
//                  ~ Sum { w_k * dPhi_i(xi_k)/dxi dPhi_j(xi_k)/dxi }   (Gaussian quadrature)
//          where w_k are the weights for the basis associated with each node, xi_k.
// ARGS   : none
// RETURNS: TRUE on success, else FALSE
//************************************************************************************
GBOOL GLLBasis::ComputeStiffMatrix()
{
  if ( bNeedNodes || bNeedWeights ) {
    if ( !ComputeNodes() ) return FALSE;
  }
  if ( bNeedDerivMatrix )
    if ( !ComputeDerivMatrix() ) return FALSE;

  GINT  i, j, k;

  StiffMatrix = 0.0;

  // Could use matrix algebra here....
  for ( i=0; i<Np+1; i++ ) {
    for ( j=0; j<Np+1; j++ ) {
       for ( k=0; k<Np+1; k++ ) {
         StiffMatrix (i,j) += Weights(k)*dPhi(k,i)*dPhi(k,j) ;  
       }
       StiffMatrixd(i,j)  = StiffMatrix(i,j);
    } 
  } 
  return TRUE;
} // end of method ComputeStiffMatrix


//************************************************************************************
//************************************************************************************
// METHOD : Solve
// DESC   : Computes all quantities, nodes, weights, mass and
//          stiffness matrices, and derivative matrices Phi, and dPhi.
// ARGS   : none
// RETURNS: TRUE on success, else FALSE
//************************************************************************************
GBOOL GLLBasis::Solve()
{
  if ( bNeedNodes || bNeedWeights ) {
    if ( !ComputeNodes() ) return FALSE;
  }

  if ( bNeedDerivMatrix )
    if ( !ComputeDerivMatrix() ) return FALSE;

  if ( !ComputeStiffMatrix() ) return FALSE; // stiffness matrix computed; dPhi also computed.

  return TRUE;

} // end of method Solve



//************************************************************************************
//************************************************************************************
// METHOD : EvalBasis (1)
// DESC   : 
//                h(xi)_i = 1/(N(N+1)*L_N(xi_i))  * (1-xi^2)*dL_N(xi)/dxi / (xi-xi_i)
// ARGS   : i  : which polynomial to evaluat (0,... Np)
//          eta: reference interval value at which to evaluate 
// RETURNS: scalar result of evaluation   
//************************************************************************************
GDOUBLE GLLBasis::EvalBasis (GINT  i, GDOUBLE eta)
{
  GQUAD ppn_i, pm1, pdm1, pm2, pdm2, ppn_xi, pder_xi, del;
  GQUAD fact=-1.0/(Np*(Np+1.0)), gfact, fRet, xi=(GQUAD)eta;

  if ( (bNeedNodes || bNeedWeights) &&  !ComputeNodes() ) {
    cout << "GLLBasis::EvalBasis (1): basis data incomplete" << endl;
    exit(1);
  }

  fRet = 1.0;
  del  = xi-xiNodes(i);
  if ( (GDOUBLE)xi < (GDOUBLE)ximin || (GDOUBLE)xi > (GDOUBLE)ximax ) fRet = 0.0;
  else if ( fabs(del) > 10.0*GTINY ) {
    ppn_i = Pn(i);
    ComputeJacobi(Np, alpha, beta, ppn_xi, pder_xi,pm1, pdm1, pm2, pdm2, xi);
    gfact = (1.0-xi*xi)/del; 
    fRet  = fact * gfact * pder_xi/ppn_i;
  }
  return (GDOUBLE)fRet;

} // end of method EvalBasis (1)


//************************************************************************************
//************************************************************************************
// METHOD     : EvalBasis (2)
// DESC   : 
//                h(xi)_j = 1/(N(N+1)*L_N(xi_j))  * (1-xi^2)*dL_N(xi)/dxi / (xi-xi_j)
// ARGS   : i   : which polynomial to evaluat (0,... Np)
//          eta : deref GVector to hold input ref interval points
//          vret: GVector to hold results of evaluations at eta
// RETURNS: pointer to vret on success; else NULLPTR  
//************************************************************************************
GVector *GLLBasis::EvalBasis (GINT  i, GVector &eta, GVector &vret)
{
  GINT  j;
  GQUAD ppn_i, pm1, pdm1, pm2, pdm2, ppn_xi, pder_xi;
  GQUAD fact=-1.0/(Np*(Np+1.0)), gfact, fRet, xi;

  if ( (bNeedNodes || bNeedWeights) &&  !ComputeNodes() ) {
    cout << "GLLBasis::EvalBasis (2): basis data incomplete" << endl;
    exit(1);
  }

  for ( j=0; j<eta->dim(); j++ ) {
    xi = eta[j];
    fRet = 1.0;
    if ( (GDOUBLE)xi < (GDOUBLE)ximin || (GDOUBLE)xi > (GDOUBLE)ximax ) fRet = 0.0;
    else if ( fabs(xi-xiNodes(i)) > 100.0*GTINY ) {
      ppn_i = Pn(i);
      ComputeJacobi(Np, alpha, beta, ppn_xi, pder_xi,pm1, pdm1, pm2, pdm2, xi);
      gfact = (1.0-xi*xi)/(xi-xiNodes(i)+GTINY); 
      fRet  = fact * gfact * pder_xi/(ppn_i+GTINY);
    }
    vret(j) = (GDOUBLE)fRet;
  }
  return &vret;

} // end of method EvalBasis (2)


//************************************************************************************
//************************************************************************************
// METHOD : EvalBasis (3)
// DESC   : Evaluates basis at input parent domain points , eta_i
//          For Gauss-Lobatto, the basis is:
//          h_j(eta) = -1/(Np*(Np+1)) * (1-eta**2) dL_Np (eta)dxi / (L_Np(xi_j) (eta-xi_j))
// ARGS   : eta  : ref interval points at which to eval all polynomials
//          mret : GMatrix return for evaluated polynomials
// RETURNS: pointer to GMatrix, M_ij = h_j(eta_i) on success; else NULLPTR
//************************************************************************************
GMatrix *GLLBasis::EvalBasis (GVector &eta, GMatrix &mret)
{
  GINT  i, j;
  GQUAD ppn_j, pm1, pdm1, pm2, pdm2, ppn_xi, pder_xi;
  GQUAD fact=-1.0/(Np*(Np+1.0)), gfact, fRet, xi;

  if ( eta == NULLPTR || mret == NULLPTR ) return NULLPTR;

  if ( (bNeedNodes || bNeedWeights) &&  !ComputeNodes() ) {
    cout << "GLLBasis::EvalBasis (3): basis data incomplete" << endl;
    exit(1);
  }

  for ( i=0; i<eta->dim(); i++ ) {
    for ( j=0; j<Np+1; j++ )  {
      fRet = 1.0;
      xi    = (GQUAD) eta(i);
      if ( (GDOUBLE)xi < (GDOUBLE)ximin || (GDOUBLE)xi > (GDOUBLE)ximax ) fRet = 0.0;
      else if ( fabs(xi-xiNodes(j)) > 100.0*GTINY ) {
        ppn_j = Pn(j);
        ComputeJacobi(Np, alpha, beta, ppn_xi, pder_xi,pm1, pdm1, pm2, pdm2, xi);
        gfact = (1.0-xi*xi)/(xi-xiNodes(j)); 
        fRet  = fact * gfact * pder_xi/ppn_j;
      }
      mret(i,j) = fRet;
    }
  }
  return &mret;

} // end of method EvalBasis (3)


//************************************************************************************
//************************************************************************************
// METHOD : EvalBasis (4)
// DESC   : Evaluates basis at input parent domain points , eta_i
//          For Gauss-Lobatto, the basis is:
//          h_j(eta) = -1/(Np*(Np+1)) * (1-eta**2) dL_Np (eta)dxi / (L_Np(xi_j) (eta-xi_j))
// ARGS   : eta  : GDOUBLE array of ref interval points, xi_j
//          neta : num elements in eta array
//          mret : GMatrix return for evaluated polynomials
// RETURNS: pointer to return matrix, mret on success; else NULPTR 
//************************************************************************************
GMatrix *GLLBasis::EvalBasis (GDOUBLE eta[], GINT neta, GMatrix &mret)
{
  GINT  i, j;
  GQUAD ppn_j, pm1, pdm1, pm2, pdm2, ppn_xi, pder_xi;
  GQUAD fact=-1.0/(Np*(Np+1.0)), gfact, fRet, xi;

  if ( (bNeedNodes || bNeedWeights) &&  !ComputeNodes() ) {
    cout << "GLLBasis::EvalBasis (4): basis data incomplete" << endl;
    exit(1);
  }

  for ( i=0; i<neta; i++ ) {
    for ( j=0; j<Np+1; j++ )  {
      fRet = 1.0;
      xi    = (GQUAD) eta[i];
      if ( (GDOUBLE)xi < (GDOUBLE)ximin || (GDOUBLE)xi > (GDOUBLE)ximax ) fRet = 0.0;
      else if ( fabs(xi-xiNodes(j)) > 100.0*GTINY ) {
        ppn_j = Pn(j);
        ComputeJacobi(Np, alpha, beta, ppn_xi, pder_xi,pm1, pdm1, pm2, pdm2, xi);
        gfact = (1.0-xi*xi)/(xi-xiNodes(j)); 
        fRet  = fact * gfact * pder_xi/ppn_j;
      }
      mret(i,j) = fRet;
    }
  }
  return mret;

} // end of method EvalBasis (4)


//************************************************************************************
//************************************************************************************
// METHOD : EvalDBasis
// DESC   : Evaluates basis derivative at input parent domain points , eta_i
//              Deriv. is derived from :
//              h_j(eta) =  -1/(Np*(Np-1)) * (1-eta**2) dL_Np (eta)dxi / (L_Np(xi_j) (eta-xi_j))
// ARGS   : eta  : GVector of ref inteval points at which to evaluate derivative basis
//          mret : GMatrix return object to hld evaluation results
// RETURNS: pointer to mret on success; else NULLPTR
//************************************************************************************
GMatrix *GLLBasis::EvalDBasis (GVector &eta, GMatrix &mret)
{
  GINT  i, j, mm, nn;
  GQUAD ppn_j, pm1, pdm1, pm2, pdm2, ppn_xi, pder_xi, pdd;
  GQUAD fact=-1.0/(Np*(Np+1.0)), gfact, g1, g1i, g2, fRet, xi ;
  char *serr = "GLLBasis::EvalDBasis: ";
  
  if ( (bNeedNodes || bNeedWeights) &&  !ComputeNodes() ) {
    cout << "GLLBasis::EvalDBasis: basis data incomplete" << endl;
    exit(1);
  }

  nn = MIN(eta.dim(),mret.dim(1));
  mm = MIN(Np+1,mret.dim(2)); 
  for ( i=0; i<nn; i++) {
    for ( j=0; j<mm;  j++) {
      fRet = 0.0;
      xi     = (GQUAD) eta(i);
      g1     = xi - xiNodes(j);
      if ( (GDOUBLE)xi < (GDOUBLE)ximin || (GDOUBLE)xi > (GDOUBLE)ximax ) fRet = 0.0;
      else if ( fabs(g1) > 100.0*GTINY ) {
//      ComputeJacobi(Np, alpha, beta, ppn_j , pder_j ,pm1, pdm1, pm2, pdm2, xiNodes(j));
        ppn_j = Pn(j);
        ComputeJacobi(Np, alpha, beta, ppn_xi, pder_xi,pm1, pdm1, pm2, pdm2, xi);
        gfact = fact / ppn_j;
        g1i   = 1.0/g1;
        g2    = (1.0 - xi*xi)*g1i;
        pdd   = 2.*xi*pder_xi - Np*(Np + 1.)*ppn_xi; 
        fRet  = gfact * g1i * ( pdd - (2.0*xi + g2 )*pder_xi);
      } 
      mret(i,j) = fRet;
    }
  }
  return &mret;
} // end of method EvalDBasis

//************************************************************************************
//************************************************************************************
// METHOD : EvalDBasis (2)
// DESC   : Evaluates basis derivative at input parent domain points , eta_i
//          Deriv. is derived from :
//          h_j(eta) =  -1/(Np*(Np-1)) * (1-eta**2) dL_Np (eta)dxi / (L_Np(xi_j) (eta-xi_j))
// ARGS   : eta : array of reference interval points at which to evaluate
//          n   : size of eta array
//          mret: Matrix of basis evaluations, M_ij = dh_j(eta_i)/dxi
// RETURNS: Pointer to mret on success; else NULLPTR
//************************************************************************************
GMatrix *GLLBasis::EvalDBasis (GDOUBLE eta[], GINT n, GMatrix &mret)
{
  GINT  i, j, mm, nn;
  GQUAD ppn_j, pm1, pdm1, pm2, pdm2, ppn_xi, pder_xi, pdd;
  GQUAD fact=-1.0/(Np*(Np+1.0)), gfact, g1, g1i, g2, fRet, xi ;
  char *serr = "GLLBasis::EvalDBasis(2): ";
  
  if ( (bNeedNodes || bNeedWeights) &&  !ComputeNodes() ) {
    cout << "GLLBasis::EvalDBasis: basis data incomplete" << endl;
    exit(1);
  }

  nn = MIN(n,mret.dim(1));
  mm = MIN(Np+1,mret.dim(2)); 
  for ( i=0; i<nn; i++) {
    for ( j=0; j<mm;  j++) {
      fRet = 0.0;
      xi     = (GQUAD) eta[i];
      g1     = xi - xiNodes(j);
      if ( (GDOUBLE)xi < (GDOUBLE)ximin || (GDOUBLE)xi > (GDOUBLE)ximax ) fRet = 0.0;
      else if ( fabs(g1) > 100.0*GTINY ) {
//      ComputeJacobi(Np, alpha, beta, ppn_j , pder_j ,pm1, pdm1, pm2, pdm2, xiNodes(j));
        ppn_j = Pn(j);
        ComputeJacobi(Np, alpha, beta, ppn_xi, pder_xi,pm1, pdm1, pm2, pdm2, xi);
        gfact = fact / ppn_j;
        g1i   = 1.0/g1;
        g2    = (1.0 - xi*xi)*g1i;
        pdd   = 2.*xi*pder_xi - Np*(Np + 1.)*ppn_xi; 
        fRet  = gfact * g1i * ( pdd - (2.0*xi + g2 )*pder_xi);
      } 
      mret(i,j) = fRet;
    }
  }
  return &ret;
} // end of method EvalDBasis


//************************************************************************************
//************************************************************************************
// METHOD : ComputeLegendreMatrix
// DESC   : Computes matrix M_ij = P_i (xi_j),
//          where P_i is the Legendre polynomial of order i, and
//          xi_j is the j-th nodal point
// ARGS   : none.
// RETURNS: TRUE on success; else FALSE 
//************************************************************************************
GBOOL GLLBasis::ComputeLegendreMatrix()
{

  if ( (bNeedNodes || bNeedWeights) &&  !ComputeNodes() ) return FALSE;

  GINT  i, j, p;
  GQUAD ppn_i, pder_i, pm1, pdm1, pm2, pdm2;

  for ( i=0; i<Np+1; i++ ) {
    for ( j=0; j<Np+1; j++ ) {
      p = i;
      ComputeJacobi(p, 0.0, 0.0, ppn_i , pder_i ,pm1, pdm1, pm2, pdm2, xiNodes(j));
//    LegMatrix(i,j) = ppn_i * Weights(j);
      LegMatrix(i,j) = ppn_i ;
    }
  }

  return TRUE;

} // end of method ComputeLegendreMatrix


