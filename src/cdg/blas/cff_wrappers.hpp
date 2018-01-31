//************************************************************************************
// Module       : cff_wrappers.hpp
// Date         : 1/1/18 (DLR)
// Description  : Wrappers for the Linear algebra Fortran functions
// Copyright    : Copyright 2018-2018. Colorado State University. All rights reserved
// Derived From : none.
//************************************************************************************
#if defined(MPI)
# error 'MTK:: MPI defined here--1'
#endif
#if defined(MPI)
# error 'MTK:: MPI defined here--2'
#endif
#include "cff_wrappers_extern.h"
#if defined(MPI)
# error 'MTK:: MPI defined here--3'
#endif
#include "blas_dinterface.hpp"
#if defined(MPI)
# error 'MTK:: MPI defined here--4'
#endif

#if !defined(CFF_WRAPPERS_HPP)
#define CFF_WRAPPERS_HPP

void w_dmxm     (GDOUBLE C[], GDOUBLE A[], GINT nai, GINT naj, GDOUBLE B[], GINT nbi, GINT nbj, GINT isz);
void w_dmxmcf   (GDOUBLE C[], GDOUBLE A[], GINT nai, GINT naj, GDOUBLE B[], GINT nbi, GINT nbj, GINT isz);
void w_dmxv     (GDOUBLE y[], GINT ny, GDOUBLE A[], GDOUBLE x[], GINT nx, GINT isz);
void w_dmxDm    (GDOUBLE C[], GDOUBLE A[],GINT nai, GINT naj, GDOUBLE b[], GINT nb, GINT isz);
void w_dDmxm    (GDOUBLE C[], GDOUBLE a[], GINT na, GDOUBLE B[], GINT nbi, GINT nbj, GINT isz);
void w_dzaxpby  (GDOUBLE z[], GDOUBLE x[], GDOUBLE a, GDOUBLE y[], GDOUBLE b, GINT nxy, GINT isz);
void w_dxaxpby  (GDOUBLE x[], GDOUBLE a, GDOUBLE y[], GDOUBLE b, GINT nxy, GINT isz);
void w_dzvxvpt  (GDOUBLE z[], GDOUBLE x[], GDOUBLE y[], GINT nxy, GINT isz);
void w_dvvxvpt  (GDOUBLE x[], GDOUBLE y[], GINT nxy, GINT isz);
void w_dvvxvptpv(GDOUBLE x[], GDOUBLE y[], GDOUBLE z[], GDOUBLE cz, GINT nxy, GINT isz);
void w_dvpvxvpt (GDOUBLE z[], GDOUBLE x[], GDOUBLE y[], GDOUBLE cz, GINT nxy, GINT isz);
void w_ddot     (GDOUBLE *dot, GDOUBLE x[], GDOUBLE y[], GINT nxy, GINT isz);
void w_dcopy    (GDOUBLE x[], GDOUBLE y[], GINT nxy, GINT isz);
void w_doop     (GDOUBLE u[], GDOUBLE g[], GINT ig[], GINT nm, GINT nop);  

inline void w_doop(GDOUBLE u[], GDOUBLE g[], GINT ig[], GINT nm, GINT nop)
{
#if defined(_LINUX) 
  doop_(u, g, ig, &nm, &nop);
#elif defined(_AIX)
  doop (u, g, ig, &nm, &nop);
#elif defined(_MACOSX)
  doop_ (u, g, ig, &nm, &nop);
#else
  DOOP_(u, g, ig, &nm, &nop);
#endif
} // end of function wrappers w_doop

inline void w_dmxm(GDOUBLE C[], GDOUBLE A[], GINT nai, GINT naj, GDOUBLE B[], GINT nbi, GINT nbj, GINT isz)
{
#if defined(_GBLAS)
#  if defined(_LINUX) 
  dmxm_(C, A, &nai, &naj, B, &nbi, &nbj, &isz);
#  elif defined(_AIX)
  dmxm (C, A, &nai, &naj, B, &nbi, &nbj, &isz);
#  elif defined(_MACOSX)
  dmxm_(C, A, &nai, &naj, B, &nbi, &nbj, &isz);
#  else
  DMXM_(C, A, &nai, &naj, B, &nbi, &nbj, &isz);
#  endif
#else
  blas_dmxm (C, A,  nai,  naj, B,  nbi,  nbj);
#endif
} // end of function wrappers w_dmxm

inline void w_dmxmcf(GDOUBLE C[], GDOUBLE A[], GINT nai, GINT naj, GDOUBLE B[], GINT nbi, GINT nbj, GINT isz)
{
#if defined(_GBLAS)
#  if defined(_LINUX) 
  dmxmcf_(C, A, &nai, &naj, B, &nbi, &nbj, &isz);
#  elif defined(_AIX)
  dmxmcf (C, A, &nai, &naj, B, &nbi, &nbj, &isz);
#  elif defined(_MACOSX)
  dmxmcf_(C, A, &nai, &naj, B, &nbi, &nbj, &isz);
#  else
  DMXMCF_(C, A, &nai, &naj, B, &nbi, &nbj, &isz);
#  endif
#else
  blas_dmxm (C, A,  nai,  naj, B,  nbi,  nbj);
#endif
} // end of function wrappers w_dmxmcf

inline void w_dmxv   (GDOUBLE y[], GINT ny, GDOUBLE A[], GDOUBLE x[], GINT nx, GINT isz)
{
#if defined(_LINUX)
  dmxv_  (y, &ny, A, x, &nx, &isz);
#elif defined(_AIX)
  dmxv   (y, &ny, A, x, &nx, &isz);
#elif defined(_MACOSX)
  dmxv_  (y, &ny, A, x, &nx, &isz);

#else
  DMXV_  (y, &ny, A, x, &nx, &isz);
#endif
} // end of function wrappers w_dmxv

inline void w_dmxDm  (GDOUBLE C[], GDOUBLE A[],GINT nai, GINT naj, GDOUBLE b[], GINT nb, GINT isz)
{
#if defined(_LINUX) 
   dmxdm_ (C, A, &nai, &naj, b, &nb, &isz);
#elif defined(_AIX)
   dmxdm  (C, A, &nai, &naj, b, &nb, &isz);
#elif defined(_MACOSX)
   dmxdm_ (C, A, &nai, &naj, b, &nb, &isz);

#else
   MXDM_ (C, A, &nai, &naj, b, &nb, &isz);
#endif

}

inline void w_daApbB  (GDOUBLE C[], GDOUBLE A[], GDOUBLE B[], GINT n, GINT m, GDOUBLE a, GDOUBLE b, GINT isz)
{
#if defined(_LINUX) 
   daapbb_ (C, A, B, &m, &n, &a, &b, &isz);
#elif defined(_AIX)
   daapbb  (C, A, B, &m, &n, &a, &b, &isz);
#elif defined(_MACOSX)
   daapbb_ (C, A, B, &m, &n, &a, &b, &isz);
#else
   DAAPBB_ (C, A, B, &m, &n, &a, &b, &isz);
#endif

}

inline void w_dDmxm  (GDOUBLE C[], GDOUBLE a[], GINT na, GDOUBLE B[], GINT nbi, GINT nbj, GINT isz)
{
#if defined(_LINUX)
  ddmxm_ (C, a, &na, B, &nbi, &nbj, &isz);
#elif defined(_AIX)
  ddmxm  (C, a, &na, B, &nbi, &nbj, &isz);
#elif defined(_MACOSX)
  ddmxm_ (C, a, &na, B, &nbi, &nbj, &isz);
#else
  DDMXM_ (C, a, &na, B, &nbi, &nbj, &isz);
#endif
} 

inline void w_dzaxpby(GDOUBLE z[], GDOUBLE x[], GDOUBLE a, GDOUBLE y[], GDOUBLE b, GINT nxy, GINT isz)
{
#if defined(_LINUX) 
  dzaxpby_(z, x, &a, y, &b, &nxy, &isz);
#elif defined(_AIX)
  dzaxpby (z, x, &a, y, &b, &nxy, &isz);
#elif defined(_MACOSX)
  dzaxpby_(z, x, &a, y, &b, &nxy, &isz);
#else
  DZAXPBY_(z, x, &a, y, &b, &nxy, &isz);
#endif
} 

inline void w_dxaxpby(GDOUBLE x[], GDOUBLE a, GDOUBLE y[], GDOUBLE b, GINT nxy, GINT isz)
{
#if defined(_LINUX) 
  dxaxpby_(x, &a, y, &b, &nxy, &isz);
#elif defined(_AIX)
  dxaxpby (x, &a, y, &b, &nxy, &isz);
#elif defined(_MACOSX)
  dxaxpby_(x, &a, y, &b, &nxy, &isz);
#else
  DXAXPBY)(x, &a, y, &b, &nxy, &isz);
#endif
} 

inline void w_dzvxvpt(GDOUBLE z[], GDOUBLE x[], GDOUBLE y[], GINT nxy, GINT isz)
{
#if defined(_LINUX) 
  dzvxvpt_(z, x, y, &nxy, &isz);
#elif defined(_AIX)
  dzvxvpt (z, x, y, &nxy, &isz);
#elif defined(_MACOSX)
  dzvxvpt_(z, x, y, &nxy, &isz);
#else
  DZVXVPT_(z, x, y, &nxy, &isz);
#endif
} 

inline void w_dvvxvpt(GDOUBLE x[], GDOUBLE y[], GINT nxy, GINT isz)
{
#if defined(_LINUX)
  dvvxvpt_(x, y, &nxy, &isz);
#elif defined(_AIX)
  dvvxvpt (x, y, &nxy, &isz);
#elif defined(_MACOSX)
  dvvxvpt_(x, y, &nxy, &isz);
#else
  DVVXVPT_(x, y, &nxy, &isz);
#endif
} 


inline void w_dvvxvptpv(GDOUBLE x[], GDOUBLE y[], GDOUBLE z[], GDOUBLE cz, GINT nxy, GINT isz)
{
#if defined(_LINUX)
  dvvxvptpv_(x, y, z, &cz, &nxy, &isz);
#elif defined(_AIX)
  dvvxvptpv (x, y, z, &cz, &nxy, &isz);
#elif defined(_MACOSX)
  dvvxvptpv_(x, y, z, &cz, &nxy, &isz);
#else
  DVVXVPTPV_(x, y, z, &cz, &nxy, &isz);
#endif
} 

inline void w_dvpvxvpt(GDOUBLE z[], GDOUBLE x[], GDOUBLE y[], GDOUBLE cz, GINT nxy, GINT isz)
{
#if defined(_LINUX)
  dvpvxvpt_(z, x, y, &cz, &nxy, &isz);
#elif defined(_AIX)
  dvpvxvpt (z, x, y, &cz, &nxy, &isz);
#elif defined(_MACOSX)
  dvpvxvpt_(z, x, y, &cz, &nxy, &isz);
#else
  DVPVXVPT_(z, x, y, &cz, &nxy, &isz);
#endif
}


inline void w_ddot   (GDOUBLE *dot, GDOUBLE x[], GDOUBLE y[], GINT nxy, GINT isz)
{
#if defined(_LINUX) 
  ddotg_ (dot, x, y, &nxy, &isz);
#elif defined(_AIX)
  ddotg  (dot, x, y, &nxy, &isz);
#elif defined(_MACOSX)
  ddotg_ (dot, x, y, &nxy, &isz);
#else
  DDOTG_ (dot, x, y, &nxy, &isz);
#endif
} 

inline void w_dcopy(GDOUBLE x[], GDOUBLE y[], GINT nxy, GINT isz)
{
#if defined(_LINUX) 
  dcopy_(x, y, &nxy, &isz);
#elif defined(_AIX)
  dcopy (x, y, &nxy, &isz);
#elif defined(_MACOSX)
  dcopy_(x, y, &nxy, &isz);
#else
  DCOPY_(x, y, &nxy, &isz);
#endif
} 

inline void w_isassign(GDOUBLE u[], GDOUBLE g[], GINT ig[], GINT nm, GINT nop)
{
#if defined(_LINUX) 
  isassign_(u, g, ig, &nm, &nop);
#elif defined(_AIX)
  isassign (u, g, ig, &nm, &nop);
#elif defined(_MACOSX)
  isassign_(u, g, ig, &nm, &nop);
#else
  ISASSIGN_(u, g, ig, &nm, &nop);
#endif
} // end of function wrappers w_isassign

inline void w_isum(GDOUBLE &res, GDOUBLE u[], GINT ig[], GINT ne)
{
#if defined(_LINUX) 
  isum_(&res, u, ig, &ne);
#elif defined(_AIX)
  isum (&res, u, ig, &ne);
#elif defined(_MACOSX)
  isum_(&res, u, ig, &ne);
#else
  ISUM_(u, g, ig, &nm, &nop);
#endif
} // end of function wrappers w_isum

inline void w_iprod(GDOUBLE &res, GDOUBLE u[], GINT ig[], GINT ne)
{
#if defined(_LINUX) 
  iprod_(&res, u, ig, &ne);
#elif defined(_AIX)
  iprod (&res, u, ig, &ne);
#elif defined(_MACOSX)
  iprod_(&res, u, ig, &ne);
#else
  IPROD_(u, g, ig, &nm, &nop);
#endif
} // end of function wrappers w_iprod

inline void w_imax(GDOUBLE &res, GDOUBLE u[], GINT ig[], GINT ne)
{
#if defined(_LINUX) 
  imax_(&res, u, ig, &ne);
#elif defined(_AIX)
  imax (&res, u, ig, &ne);
#elif defined(_MACOSX)
  imax_(&res, u, ig, &ne);
#else
  IMAX_(u, g, ig, &nm, &nop);
#endif
} // end of function wrappers w_imax

inline void w_imin(GDOUBLE &res, GDOUBLE u[], GINT ig[], GINT ne)
{
#if defined(_LINUX) 
  imin_(&res, u, ig, &ne);
#elif defined(_AIX)
  imin (&res, u, ig, &ne);
#elif defined(_MACOSX)
  imin_(&res, u, ig, &ne);
#else
  IMIN_(u, g, ig, &nm, &nop);
#endif
} // end of function wrappers w_imin

#endif
