//==================================================================================
// Module       : glbasis.cpp
// Date         : 1/18/18 (DLR)
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
#include "glbasis.hpp"
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <stdio.h>


//************************************************************************************
//************************************************************************************
// METHOD : Contructor (1)
// DESC   : Instantiate with polynomial order, and maxOrder
//          
// ARGS   :
// RETURNS: TRUE on success, else FALSE
//************************************************************************************
GLBasis::GLBasis(GINT  inOrder, GINT  MaxOrder)
: GLLBasis(inOrder,MaxOrder)
{
  kstop = 128;
  eps   = 1.0e-16;
  alpha = 0.0;
  beta  = 0.0;
  Resize(NpMax);

} // end of constructor method


//************************************************************************************
//************************************************************************************
// METHOD : Contructor (2)
// DESC   : Default constructor
//          
// ARGS   :
// RETURNS:
//************************************************************************************
GLBasis::GLBasis()
: GLLBasis()
{
  kstop = 128;
  eps   = 1.0e-16;
  alpha = 0.0;
  beta  = 0.0;
  Resize(NpMax);

} // end of constructor method

//************************************************************************************
//************************************************************************************
// METHOD : Contructor (3)
// DESC   : Instantiate with poly order (same for maxOrder)
//          
// ARGS   :
// RETURNS:
//************************************************************************************
GLBasis::GLBasis(GINT  Order)
: GLLBasis(Order)
{
  kstop = 128;
  eps   = 1.0e-16;
  alpha = 0.0;
  beta  = 0.0;
  Resize(NpMax);
} // end of constructor method



//************************************************************************************
//************************************************************************************
// METHOD : Copy consructor
// DESC   : 
//          
// ARGS   :
// RETURNS: 
//************************************************************************************
GLBasis::GLBasis(const GLBasis &b)
{
   // copy data:
    alpha    = b.alpha;
    beta     = b.beta;
    ximin    = b.ximin;
    ximax    = b.ximax;
    Np       = b.Np;
    NpMax    = b.NpMax;
    bNeedNodes = b.bNeedNodes;
    bNeedWeights = b.bNeedWeights;
    bNeedDerivMatrix = b.bNeedDerivMatrix;
    bNeedBasis   = b.bNeedBasis; 
    bNeedDBasis  = b.bNeedDBasis;
    bNeedLegMat  = b.bNeedLegMat;



    //  matrices, and node data:
    xiNodes     = b.xiNodes;
    xid         = b.xid;
    Weights     = b.Weights;
    weightsd    = b.weightsd;
    Pn          = b.Pn;
    dPn         = b.dPn;
    Phi         = b.Phi;
    dPhi        = b.dPhi;
    MassMatrixd = b.MassMatrixd;
    StiffMatrix = b.StiffMatrix;
    StiffMatrixd= b.StiffMatrixd;
    LegMatrix   = b.LegMatrix;

}

//************************************************************************************
//************************************************************************************
// METHOD : Destructor
// DESC   : 
// ARGS   :
// RETURNS: 
//************************************************************************************
GLBasis::~GLBasis()
{
}

//************************************************************************************
//************************************************************************************
// METHOD : Assignment operator
// DESC   : 
// ARGS   :
// RETURNS: 
//************************************************************************************
void GLBasis::operator=(const GLBasis &b)
{
  if ( &b != this ) 
  {
   // copy data:
    ximin    = b.ximin;
    ximax    = b.ximax;
    Np       = b.Np;
    NpMax    = b.NpMax;
    bNeedNodes = b.bNeedNodes;
    bNeedWeights = b.bNeedWeights;
    bNeedDerivMatrix = b.bNeedDerivMatrix;
    bNeedBasis   = b.bNeedBasis; 
    bNeedDBasis  = b.bNeedDBasis;
    bNeedLegMat  = b.bNeedLegMat;

    //  matrices, and node data:
    xiNodes     = b.xiNodes;
    xid         = b.xid;
    Weights     = b.Weights;
    weightsd    = b.weightsd;
    Pn          = b.Pn;
    dPn         = b.dPn;
    Phi         = b.Phi;
    dPhi        = b.dPhi;
    MassMatrixd = b.MassMatrixd;
    StiffMatrix = b.StiffMatrix;
    StiffMatrixd= b.StiffMatrixd;
    LegMatrix   = b.LegMatrix;
  }

}


//************************************************************************************
//************************************************************************************
// METHOD : Resize
// DESC   : Resizes dynamically allocated quantities
//          if required
// ARGS   : GINT newOrder
// RETURNS: TRUE on success, else FALSE
//************************************************************************************
GBOOL GLBasis::Resize(GINT  newOrder)
{

  // No resizing necessary if already less than
  // previously allocated quantities:

  NpMax = Np = newOrder;

  GINT  NN = Np+1;

  //  Resize xiNodes:
  xiNodes.Resize(NN);

  //  Resize Weights:
  Weights.Resize(NN);

  //  Resize Pn, derivatives::
  Pn.Resize(NN);
  dPn.Resize(NN);


  //  Resize basis Phi:
  Phi.Resize(NN,NN);

  //  Resize basis dPhi:
  dPhi.Resize(NN,NN);

  //  Resize MassMatrix:
  MassMatrixd.Resize(NN,NN);

  //  Resize StiffMatrix:
  StiffMatrix.Resize(NN,NN);
  
  return TRUE;

} // end of method Resize


//************************************************************************************
//************************************************************************************
// METHOD : ComputeNodes 
// DESC   : Computes nodes and weights based on a Gauss-Legendre integration scheme
//          Note: nodes are computed from smalles to largest xi
//          Note: this method was taken largely from Canuto et al, 1987, Appendix C.
//          Note: Legendre polynomials of order Np, Np-1 and Np-2 and their
//                derivatives are computed here for use here and in 
//                computing integrals. These polynomials are used in
//                forming basis functions, Phi, and their derivatives,
//                dPhi, evaluated at each node.
// ARGS   : none
// RETURNS: TRUE on success, else FALSE
//************************************************************************************
GBOOL GLBasis::ComputeNodes()
{
  if ( !bNeedNodes && !bNeedWeights ) return TRUE;
  if ( Np < 1 ) return FALSE;
#if 0
  if ( Np == 1 ) {
    xiNodes(0) = ximin;
    xiNodes(1) = ximax;
    bNeedNodes = FALSE;
    return(ComputeWeights());
  }
#endif
 
  GINT  i, j, k, np1, nh;
//GQUAD rv;
  GQUAD det, pnp, pnp1p, pnp1m, rp, rm, dth, cd, sd, cs, ss, 
        pn , pnp1, pdnp1, pdn, pnm1, pdnm1, poly, pder,
        pdnp1p, pnm1p,  pdnp1m, pnm, pdnm, pnm1m, pdnp,
        recsum, x, delx, cssave, error, a, b;
 
  alpha = 0.0;
  beta  = 0.0;  //alpha=beta=0 ==> Legendre fcn basis 

//rv    = 1.0 + alpha;
  np1   = Np + 1;
  nh    = np1 / 2;

  ComputeJacobi(np1,alpha, beta, pnp1p, pdnp1p, pnp, pdnp, pnm1p, pdnm1, ximax);
  ComputeJacobi(np1,alpha, beta, pnp1m, pdnp1m, pnm, pdnm, pnm1m, pdnm1, ximin);
  det  = pnp*pnm1m - pnm*pnm1p;
  rp   = -pnp1p;
  rm   = -pnp1m;
  a    = ( rp*pnm1m - rm*pnm1p ) / det;
  b    = ( rm*pnp   - rp*pnm   ) / det;

  // set up recursion relation for the initial guesses for roots:
  dth = 3.1415926535898/( 2.0*Np + 1.0 );
  cd  = cos(2.0*dth);
  sd  = sin(2.0*dth);
  cs  = cos    (dth);
  ss  = sin    (dth);

  // compute the roots by polynomial deflation:
  // together with Newton-Raphson....
  for ( j=0; j<nh; j++ ) {
    x = cs;
    error = 1.0;
    for ( k=0; k<kstop && error>eps; k++ ) {
      ComputeJacobi(np1, alpha, beta, pnp1, pdnp1, pn, pdn, pnm1, pdnm1, x);
//    poly = pnp1 + a*pn + b*pnm1;
//    pder = pdnp1 + a*pdn + b*pdnm1;
      poly = pnp1 ;
      pder = pdnp1;
      recsum = 0.0;
      for ( i=0; i<j; i++ ) {
        recsum += (  1.0/(x - xiNodes(i)) );
      }
      delx = -poly/(pder - recsum*poly);
      x += delx;
      error = fabs(delx)/ fabs(x);    // MIN(fabs(ximin),fabs(ximax));
    }
    xiNodes(j) = x;
    cssave     = cs*cd - ss*sd;
    ss         = cs*sd + ss*cd;
    cs         = cssave;
  }


#if 1
  // Symmetrize the nodes:  NOTE: this is correct only for the cases
  // where (xmin,xmax) = (-y, y), or something like this.
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
#endif

  for ( i=0; i<Np+1; i++ ) xid[i] = (GDOUBLE) xiNodes[i];

  return(ComputeWeights());

} // end of method ComputeNodes


//************************************************************************************
//************************************************************************************
// METHOD : ComputeWeights
// DESC   : 
//          NOTE: method not really intended to be called publicly; 
//          it should be called only when ComputeNodes is called. 
//          For other bases, this method will change, however.
// ARGS   : none
// RETURNS: TRUE on success, else FALSE
//************************************************************************************
GBOOL GLBasis::ComputeWeights()
{
  if ( Np < 1 ) return FALSE;
#if 0
  if ( Np == 1 ) {
    Weights(0) = 0.5;
    Weights(1) = 0.5;
    bNeedWeights = FALSE;
    return TRUE;
  }
#endif
 
  GINT  i, np1=Np+1;
  GQUAD fact;
  GQUAD ppn, pder, pm1, pdm1, pm2, pdm2;

  for ( i=0; i<np1; i++ ) {
    ComputeJacobi(np1, alpha, beta, ppn, pder,pm1, pdm1, pm2, pdm2, xiNodes(i));
    Pn (i) = ppn;    
    dPn(i) = pder;    
    fact   = 2.0/( 1.0 - xiNodes(i)*xiNodes(i) );
    Weights(i) = (pder==0.0)?0.0:fact/(pder*pder); 
  }
  for ( i=0; i<Np+1; i++ ) weightsd[i] = (GDOUBLE) Weights[i];
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
GBOOL GLBasis::ComputeBasisAtNodes()
{
  if ( bNeedNodes || bNeedWeights )
    if ( !ComputeNodes() ) return NULLPTR;


  GINT  i, j, np1=Np+1;
  GQUAD ppn_i, pder_i, pm1, pdm1, pm2, pdm2, ppn_j, pder_j;
  GQUAD fact=1.0/(Np*(Np.0)), gfact;

  Phi = 0.0;

  for ( i=0; i<np1; i++ ) {
    Phi(i,i) = 1.0;
    for ( j=0; j< i; j++ ) {
      gfact = 1.0/(xiNodes(j)-xiNodes(i));
      Phi(i,j) = 1.0/dPn(i) * gfact * Pn(j);
    }
    for ( j=i+1; j<np1; j++ ) {
      gfact = 1.0/(xiNodes(j)-xiNodes(i));
      Phi(i,j) = 1.0/dPn(i) * gfact * Pn(j);
    }
  }
  bNeedBasis = FALSE;

} // end of method ComputeBasisAtNodes
#endif

//************************************************************************************
//************************************************************************************
// METHOD : ComputeDerivMatrix
// DESC   :
//          NOTE: ComputeWeights (ComputeNodes) must have been called
//                prior to entering this method, s.t. the Pn_i have been
//                calculated. 
// ARGS   : none
// RETURNS: TRUE on success, else FALSE
//************************************************************************************
GBOOL GLBasis::ComputeDerivMatrix()
{
  if ( bNeedWeights ) return FALSE;

  GINT  l, j, np1=Np+1;
  GQUAD delxi, xi, xisq, d2Pn;

  dPhi = 0.0;
  // Note: index j sweeps over basis number; l sweeps over node number

  for ( j=0; j<np1; j++ ) {
    for ( l=0; l<j; l++ ) {
      delxi      = xiNodes(l) - xiNodes(j);
      dPhi (l,j) =  dPn(l)/(dPn(j)*delxi);
    }
    for ( l=j+1; l<np1; l++ ) {
      delxi       = xiNodes(l) - xiNodes(j);
      dPhi  (l,j) =  dPn(l)/(dPn(j)*delxi);
    }
  }
  for ( l=0; l<np1; l++ ) {
    xi          = xiNodes(l);
    xisq        = xi * xi;
    d2Pn        = xi*dPn(l)/(1.0-xisq);
    dPhi  (l,l) =  d2Pn  / dPn(l);
  }


  bNeedDerivMatrix = FALSE;
  return TRUE;

} // end of method ComputeDerivMatrix


//************************************************************************************
//************************************************************************************
// METHOD : EvalBasis (1)
// DESC   : 
// ARGS   : i  : which polynomial to evaluate
//          eta: scalar ref interval value at which to evaluate ith polynomial
// RETURNS:  
//************************************************************************************
GDOUBLE GLBasis::EvalBasis (GINT  i, GDOUBLE eta)
{
  GINT  np1=Np+1;
  GQUAD pm1, pdm1, pm2, pdm2, ppn_xi, pder_xi;
  GQUAD gfact, fRet, xi=(GQUAD)eta;


  if ( bNeedNodes || bNeedWeights &&  !ComputeNodes() ) {
    cout << "GLBasis::EvalBasis (1): basis data incomplete" << endl;
    exit(1);
  }

  fRet = 1.0;
//if ( (GDOUBLE)xi < (GDOUBLE)xiNodes.Min() || (GDOUBLE)xi > (GDOUBLE)xiNodes.Max()) fRet = 0.0;
  if ( (GDOUBLE)xi < (GDOUBLE)ximin || (GDOUBLE)xi > (GDOUBLE)ximax ) fRet = 0.0;
  else if ( fabs(xi-xiNodes(i)) > 10.0*GTINY ) {
  //ComputeJacobi(np1, alpha, beta, ppn_i , pder_i ,pm1, pdm1, pm2, pdm2, xiNodes(i));
    ComputeJacobi(np1, alpha, beta, ppn_xi, pder_xi,pm1, pdm1, pm2, pdm2, xi);
    gfact = 1.0/(xi-xiNodes(i)); 
    fRet  = gfact * ppn_xi/(dPn(i)+GTINY);
  }
  return (GDOUBLE)fRet;

} // end of method EvalBasis (1)


//************************************************************************************
//************************************************************************************
// METHOD : EvalBasis (2)
// DESC   :
// ARGS   :
// ARGS   : i   : which polynomial to evaluate
//          eta : array of ref interval values at which to evaluate ith polynomial
//          vret: GVector return object for evluations
// RETURNS: pointer to vret on success; else NULLPTRPTR
//************************************************************************************
GVector *GLBasis::EvalBasis (GINT  i, GVector &eta, GVector &vret)
{
  GINT  j, np1=Np+1;
  GQUAD pm1, pdm1, pm2, pdm2, ppn_xi, pder_xi;
  GQUAD gfact, fRet, xi;

  if ( bNeedNodes || bNeedWeights &&  !ComputeNodes() ) {
    cout << "GLBasis::EvalBasis (2): basis data incomplete" << endl;
    exit(1);
  }

  for ( j=0; j<eta->dim(); j++ )
  {
    xi = (*eta)(j);
    fRet = 1.0;
    if ( (GDOUBLE)xi < (GDOUBLE)xiNodes.Min() || (GDOUBLE)xi > (GDOUBLE)xiNodes.Max()) fRet = 0.0;
    else if ( fabs(xi-xiNodes(i)) > 10.0*GTINY ) {
      ComputeJacobi(np1, alpha, beta, ppn_xi, pder_xi,pm1, pdm1, pm2, pdm2, xi);
      gfact = 1.0/(xi-xiNodes(i)); 
      fRet  = gfact * ppn_xi/(dPn(i)+GTINY);
    }
    vret(j) = (GDOUBLE)fRet;
  }
  return &vret;

} // end of method EvalBasis (2)

//************************************************************************************
//************************************************************************************
// METHOD : EvalBasis (3)
// DESC   : Evaluates basis at input parent domain points , eta_i
// ARGS   :
//          eta : GVevtor of ref interval values at which to evaluate all polynomials
//          neta: num elements in eta array
//          mret: return matrix M_ij = dh_j(eta_i)/dxi
// RETURNS: pointer to mret on success; else NULLPTR
//************************************************************************************
GMatrix *GLBasis::EvalBasis (GVector &eta, GMatrix &mret)
{

  GINT  i, j, np1=Np+1, imax, jmax;
  GQUAD pm1, pdm1, pm2, pdm2, ppn_xi, pder_xi;
  GQUAD gfact, fRet, xi;

  if ( bNeedNodes || bNeedWeights &&  !ComputeNodes() ) {
    cout << "GLBasis::EvalBasis (3): basis data incomplete" << endl;
    exit(1);
  }

  imax = MIN(eta.dim(),mret.dim(1));
  jmax = MIN(np1,mret.dim(2));
  for ( i=0; i<imax; i++) { 
    xi    = (GQUAD) eta(i);
    for ( j=0; j<jmax; j++) {
      fRet = 1.0;
//    if ( (GDOUBLE)xi < (GDOUBLE)xiNodes.Min() || (GDOUBLE)xi > (GDOUBLE)xiNodes.Max()) fRet = 0.0;
      if ( (GDOUBLE)xi < (GDOUBLE)ximin || (GDOUBLE)xi > (GDOUBLE)ximax ) fRet = 0.0;
      else if ( fabs(xi-xiNodes(j)) > 10.0*GTINY ) {
        ComputeJacobi(np1, alpha, beta, ppn_xi, pder_xi,pm1, pdm1, pm2, pdm2, xi);
        gfact = 1.0/(xi-xiNodes(j)); 
        fRet  = gfact * ppn_xi/(dPn(j)+GTINY);
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
// ARGS   :
//          eta : array of ref interval values at which to evaluate all polynomials
//          neta: num elements in eta array
//          mret: return matrix M_ij = dh_j(eta_i)/dxi
// RETURNS: pointer to mret on success; else NULLPTR
//************************************************************************************
GMatrix *GLBasis::EvalBasis (GDOUBLE eta[], GINT neta, GMatrix &mret)
{
  if ( mret == NULLPTR ) return NULLPTR;

  GINT  i, j, np1=Np+1, imax, jmax;
  GQUAD pm1, pdm1, pm2, pdm2, ppn_xi, pder_xi;
  GQUAD gfact, fRet, xi;

  if ( bNeedNodes || bNeedWeights &&  !ComputeNodes() ) {
    cout << "GLBasis::EvalBasis (4): basis data incomplete" << endl;
    exit(1);
  }

  imax = MIN(neta,mret.dim(1));
  jmax = MIN(np1,mret.dim(2));
  for ( i=0; i<imax; i++) { 
    xi    = (GQUAD) eta[i];
    for ( j=0; j<jmax; j++) {
      fRet = 1.0;
//    if ( (GDOUBLE)xi < (GDOUBLE)xiNodes.Min() || (GDOUBLE)xi > (GDOUBLE)xiNodes.Max()) fRet = 0.0;
      if ( (GDOUBLE)xi < (GDOUBLE)ximin || (GDOUBLE)xi > (GDOUBLE)ximax ) fRet = 0.0;
      else if ( fabs(xi-xiNodes(j)) > 10.0*GTINY ) {
        ComputeJacobi(np1, alpha, beta, ppn_xi, pder_xi,pm1, pdm1, pm2, pdm2, xi);
        gfact = 1.0/(xi-xiNodes(j)); 
        fRet  = gfact * ppn_xi/(dPn(j)+GTINY);
      }
      mret(i,j) = fRet;
    }
  }
  return mret;

} // end of method EvalBasis (3)


//************************************************************************************
//************************************************************************************
// METHOD : EvalDBasis
// DESC   : Evaluates basis derivative at input parent domain points , eta_i
// ARGS   :
//          eta : GVector ref interval values at which to evaluate all polynomials
//          mret: return matrix M_ij = dh_j(eta_i)/dxi
// RETURNS: pointer to mret on success; else NULLPTR
//************************************************************************************
GMatrix *GLBasis::EvalDBasis (GVector &eta, GMatrix &mret)
{ 
  if ( eta == NULLPTR || mret == NULLPTR ) return NULLPTR;
  
  GINT  i, j, np1=Np+1, imax, jmax;
  GQUAD pm1, pdm1, pm2, pdm2, ppn_xi, pder_xi;
  GQUAD fRet, xi, dxi;
  
  if ( bNeedNodes || bNeedWeights &&  !ComputeNodes() ) {
    cout << "GLBasis::EvalDBasis: basis data incomplete" << endl;
    exit(1);
  }

  imax = MIN(eta.dim(),mret.dim(1));
  jmax = MIN(np1,mret.dim(2));
  for ( i=0; i<imax; i++) { 
    xi    = (GQUAD) eta(i);
    for ( j=0; j<jmax; j++) {
      dxi  = xi - xiNodes(j);
      fRet = xi / ( 1.0 - xi*xi);
//    if ( (GDOUBLE)xi < (GDOUBLE)xiNodes.Min() || (GDOUBLE)xi > (GDOUBLE)xiNodes.Max()) fRet = 0.0;
      if ( (GDOUBLE)xi < (GDOUBLE)ximin || (GDOUBLE)xi > (GDOUBLE)ximax ) fRet = 0.0;
      else if ( fabs(dxi) > 10.0*GTINY ) {
        ComputeJacobi(np1, alpha, beta, ppn_xi, pder_xi,pm1, pdm1, pm2, pdm2, xi);
        fRet  = ( pder_xi*dxi - ppn_xi ) / (dPn(j)*dxi*dxi);
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
// ARGS   :
//          eta : array of ref interval values at which to evaluate all polynomials
//          neta: num elements in eta array
//          mret: return matrix M_ij = dh_j(eta_i)/dxi
// RETURNS: pointer to mret on success; else NULLPTR
//************************************************************************************
GMatrix *GLBasis::EvalDBasis (GDOUBLE eta[], GINT neta, GMatrix &mret)
{ 
  if ( mret == NULLPTR ) return NULLPTR;
  
  GINT  i, j, np1=Np+1, imax, jmax;
  GQUAD pm1, pdm1, pm2, pdm2, ppn_xi, pder_xi;
  GQUAD fRet, xi, dxi;
  
  if ( bNeedNodes || bNeedWeights &&  !ComputeNodes() ) {
    cout << "GLBasis::EvalDBasis (2): basis data incomplete" << endl;
    exit(1);
  }

  imax = MIN(neta,mret.dim(1));
  jmax = MIN(np1,mret.dim(2));
  for ( i=0; i<imax; i++) { 
    xi    = (GQUAD) eta[i];
    for ( j=0; j<jmax; j++) {
      dxi  = xi - xiNodes(j);
      fRet = xi / ( 1.0 - xi*xi);
//    if ( (GDOUBLE)xi < (GDOUBLE)xiNodes.Min() || (GDOUBLE)xi > (GDOUBLE)xiNodes.Max()) fRet = 0.0;
      if ( (GDOUBLE)xi < (GDOUBLE)ximin || (GDOUBLE)xi > (GDOUBLE)ximax ) fRet = 0.0;
      else if ( fabs(dxi) > 10.0*GTINY ) {
        ComputeJacobi(np1, alpha, beta, ppn_xi, pder_xi,pm1, pdm1, pm2, pdm2, xi);
        fRet  = ( pder_xi*dxi - ppn_xi ) / (dPn(j)*dxi*dxi);
      } 
      mret(i,j) = fRet;
    }
  }
  return &mret;
} // end of method EvalDBasis (2)


//************************************************************************************
//************************************************************************************
// METHOD : ComputeLegendreMatrix
// DESC   : Computes matrix M_ij = P_i (xi_j),
//          where P_i is the Legendre polynomial of order i, and
//          xi_j and W_j are the j-th nodal point, and weight,
//          respectively.
// ARGS   : none
// RETURNS: TRUE on success; else FALSE
//************************************************************************************
GBOOL GLBasis::ComputeLegendreMatrix()
{
  if ( !bNeedLegMat ) return TRUE;

  if ( bNeedNodes || bNeedWeights &&  !ComputeNodes() ) return FALSE;
  
  GINT  i, j;
  GQUAD ppn_i, pder_i, pm1, pdm1, pm2, pdm2;

  for ( i=0; i<Np+1; i++ ) {
    for ( j=0; j<Np+1; j++ ) {
      ComputeJacobi(i, 0.0, 0.0, ppn_i , pder_i ,pm1, pdm1, pm2, pdm2, xiNodes(j));
//    LegMatrix(i,j) =  ppn_i * Weights(j);
      LegMatrix(i,j) =  ppn_i;
    }
  }
  bNeedLegMat = FALSE;
  return TRUE;

} // end of method ComputeLegendreMatrix
