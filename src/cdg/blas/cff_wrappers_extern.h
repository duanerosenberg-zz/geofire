//************************************************************************************
// Module       : cff_wrappers_extern.h
// Date         : 1/1/18 (DLR)
// Description  : External declaration header for cff Fortran functions.
//                There should be a one-to-one correspondence between these
//                headers and their calls in cff_wrappers.cpp module.
// Copyright    : Copyright 2018-2018. Colorado State University. All rights reserved
// Derived From : none.
//************************************************************************************
#include "gtypes.h"

extern "C" {
#if defined(_LINUX) 
void doop_(GDOUBLE u[], GDOUBLE g[], GINT ig[], GINT *nm, GINT *nop);
#elif defined(_AIX)
void doop (GDOUBLE u[], GDOUBLE g[], GINT ig[], GINT *nm, GINT *nop);
#elif defined(_MACOSX)
void doop_(GDOUBLE u[], GDOUBLE g[], GINT ig[], GINT *nm, GINT *nop);
#else
void DOOP_(GDOUBLE u[], GDOUBLE g[], GINT ig[], GINT *nm, GINT *nop);
#endif

#if defined(_LINUX) 
void dmxm_(GDOUBLE C[], GDOUBLE A[], GINT *nai, GINT *naj, GDOUBLE B[], GINT *nbi, GINT *nbj, GINT *isz);
#elif defined(_AIX)
void dmxm (GDOUBLE C[], GDOUBLE A[], GINT *nai, GINT *naj, GDOUBLE B[], GINT *nbi, GINT *nbj, GINT *isz);
#elif defined(_MACOSX)
void dmxm_(GDOUBLE C[], GDOUBLE A[], GINT *nai, GINT *naj, GDOUBLE B[], GINT *nbi, GINT *nbj, GINT *isz);
#else
void DMXM_(GDOUBLE C[], GDOUBLE A[], GINT *nai, GINT *naj, GDOUBLE B[], GINT *nbi, GINT *nbj, GINT *isz);
#endif

#if defined(_LINUX) 
void dmxmcf_(GDOUBLE C[], GDOUBLE A[], GINT *nai, GINT *naj, GDOUBLE B[], GINT *nbi, GINT *nbj, GINT *isz);
#elif defined(_AIX)
void dmxmcf (GDOUBLE C[], GDOUBLE A[], GINT *nai, GINT *naj, GDOUBLE B[], GINT *nbi, GINT *nbj, GINT *isz);
#elif defined(_MACOSX)
void dmxmcf_(GDOUBLE C[], GDOUBLE A[], GINT *nai, GINT *naj, GDOUBLE B[], GINT *nbi, GINT *nbj, GINT *isz);
#else
void DMXMCF_(GDOUBLE C[], GDOUBLE A[], GINT *nai, GINT *naj, GDOUBLE B[], GINT *nbi, GINT *nbj, GINT *isz);
#endif

#if defined(_LINUX) 
void dmxv_(GDOUBLE y[], GINT *ny, GDOUBLE A[], GDOUBLE x[], GINT *nx, GINT *isz);
#elif defined(_AIX)
void dmxv (GDOUBLE y[], GINT *ny, GDOUBLE A[], GDOUBLE x[], GINT *nx, GINT *isz);
#elif defined(_MACOSX)
void dmxv_(GDOUBLE y[], GINT *ny, GDOUBLE A[], GDOUBLE x[], GINT *nx, GINT *isz);
#else
void DMXV_(GDOUBLE y[], GINT *ny, GDOUBLE A[], GDOUBLE x[], GINT *nx, GINT *isz);
#endif

#if defined(_LINUX) 
void dmxdm_(GDOUBLE C[], GDOUBLE A[], GINT *nai, GINT *naj, GDOUBLE b[], GINT *nb, GINT *isz);
#elif defined(_AIX)
void dmxdm (GDOUBLE C[], GDOUBLE A[], GINT *nai, GINT *naj, GDOUBLE b[], GINT *nb, GINT *isz);
#elif defined(_MACOSX)
void dmxdm_(GDOUBLE C[], GDOUBLE A[], GINT *nai, GINT *naj, GDOUBLE b[], GINT *nb, GINT *isz);
#else
void DMXDM_(GDOUBLE C[], GDOUBLE A[], GINT *nai, GINT *naj, GDOUBLE b[], GINT *nb, GINT *isz);
#endif

#if defined(_LINUX) 
void ddmxm_(GDOUBLE C[], GDOUBLE a[], GINT *na, GDOUBLE B[], GINT *nbi, GINT *nbj, GINT *isz);
#elif defined(_AIX)
void ddmxm (GDOUBLE C[], GDOUBLE a[], GINT *na, GDOUBLE B[], GINT *nbi, GINT *nbj, GINT *isz);
#elif defined(_MACOSX)
void ddmxm_(GDOUBLE C[], GDOUBLE a[], GINT *na, GDOUBLE B[], GINT *nbi, GINT *nbj, GINT *isz);
#else
void DDMXM_(GDOUBLE C[], GDOUBLE a[], GINT *na, GDOUBLE B[], GINT *nbi, GINT *nbj, GINT *isz);
#endif

#if defined(_LINUX) 
void daapbb_(GDOUBLE C[], GDOUBLE A[], GDOUBLE B[], GINT *n, GINT *m, GDOUBLE *a, GDOUBLE *b, GINT *isz);
#elif defined(_AIX)
void daapbb (GDOUBLE C[], GDOUBLE A[], GDOUBLE B[], GINT *n, GINT *m, GDOUBLE *a, GDOUBLE *b, GINT *isz);
#elif defined(_MACOSX)
void daapbb_(GDOUBLE C[], GDOUBLE A[], GDOUBLE B[], GINT *n, GINT *m, GDOUBLE *a, GDOUBLE *b, GINT *isz);
#else
void DAAPBB_(GDOUBLE C[], GDOUBLE A[], GDOUBLE B[], GINT *n, GINT *m, GDOUBLE *a, GDOUBLE *b, GINT *isz);
#endif

#if defined(_LINUX)
void dzaxpby_(GDOUBLE z[], GDOUBLE x[], GDOUBLE *a, GDOUBLE y[], GDOUBLE *b, GINT *nxy, GINT *isz);
#elif defined(_AIX)
void dzaxpby (GDOUBLE z[], GDOUBLE x[], GDOUBLE *a, GDOUBLE y[], GDOUBLE *b, GINT *nxy, GINT *isz);
#elif defined(_MACOSX)
void dzaxpby_(GDOUBLE z[], GDOUBLE x[], GDOUBLE *a, GDOUBLE y[], GDOUBLE *b, GINT *nxy, GINT *isz);
#else
void DZAXPBY_(GDOUBLE z[], GDOUBLE x[], GDOUBLE *a, GDOUBLE y[], GDOUBLE *b, GINT *nxy, GINT *isz);
#endif

#if defined(_LINUX) 
void dxaxpby_(GDOUBLE x[], GDOUBLE *a, GDOUBLE y[], GDOUBLE *b, GINT *nxy, GINT *isz);
#elif defined(_AIX)
void dxaxpby (GDOUBLE x[], GDOUBLE *a, GDOUBLE y[], GDOUBLE *b, GINT *nxy, GINT *isz);
#elif defined(_MACOSX)
void dxaxpby_(GDOUBLE x[], GDOUBLE *a, GDOUBLE y[], GDOUBLE *b, GINT *nxy, GINT *isz);
#else
void DXAXPBY_(GDOUBLE x[], GDOUBLE *a, GDOUBLE y[], GDOUBLE *b, GINT *nxy, GINT *isz);
#endif

#if defined(_LINUX) 
void dzvxvpt_(GDOUBLE z[], GDOUBLE x[], GDOUBLE y[], GINT *nxy, GINT *isz);
#elif defined(_AIX)
void dzvxvpt (GDOUBLE z[], GDOUBLE x[], GDOUBLE y[], GINT *nxy, GINT *isz);
#elif defined(_MACOSX)
void dzvxvpt_(GDOUBLE z[], GDOUBLE x[], GDOUBLE y[], GINT *nxy, GINT *isz);
#else
void DZVXVPT_(GDOUBLE z[], GDOUBLE x[], GDOUBLE y[], GINT *nxy, GINT *isz);
#endif

#if defined(_LINUX) 
void dvvxvpt_(GDOUBLE x[], GDOUBLE y[], GINT *nxy, GINT *isz);
#elif defined(_AIX)
void dvvxvpt (GDOUBLE x[], GDOUBLE y[], GINT *nxy, GINT *isz);
#elif defined(_MACOSX)
void dvvxvpt_(GDOUBLE x[], GDOUBLE y[], GINT *nxy, GINT *isz);
#else
void DVVXVPT_(GDOUBLE x[], GDOUBLE y[], GINT *nxy, GINT *isz);
#endif

#if defined(_LINUX) 
void dvvxvptpv_(GDOUBLE x[], GDOUBLE y[], GDOUBLE z[], GDOUBLE *cz, GINT *nxy, GINT *isz);
#elif defined(_AIX)
void dvvxvptpv (GDOUBLE x[], GDOUBLE y[], GDOUBLE z[], GDOUBLE *cz, GINT *nxy, GINT *isz);
#elif defined(_MACOSX)
void dvvxvptpv_(GDOUBLE x[], GDOUBLE y[], GDOUBLE z[], GDOUBLE *cz, GINT *nxy, GINT *isz);
#else
void DVVXVPTPV_(GDOUBLE x[], GDOUBLE y[], GDOUBLE z[], GDOUBLE *cz, GINT *nxy, GINT *isz);
#endif

#if defined(_LINUX) 
void dvpvxvpt_(GDOUBLE z[], GDOUBLE x[], GDOUBLE y[], GDOUBLE *cz, GINT *nxy, GINT *isz);
#elif defined(_AIX)
void dvpvxvpt (GDOUBLE z[], GDOUBLE x[], GDOUBLE y[], GDOUBLE *cz, GINT *nxy, GINT *isz);
#elif defined(_MACOSX)
void dvpvxvpt_(GDOUBLE z[], GDOUBLE x[], GDOUBLE y[], GDOUBLE *cz, GINT *nxy, GINT *isz);
#else
void DVPVXVPT_(GDOUBLE z[], GDOUBLE x[], GDOUBLE y[], GDOUBLE *cz, GINT *nxy, GINT *isz);
#endif

#if defined(_LINUX)
void ddotg_(GDOUBLE *dot, GDOUBLE x[], GDOUBLE y[], GINT *nxy, GINT *isz);
#elif defined(_AIX)
void ddotg (GDOUBLE *dot, GDOUBLE x[], GDOUBLE y[], GINT *nxy, GINT *isz);
#elif defined(_MACOSX)
void ddotg_(GDOUBLE *dot, GDOUBLE x[], GDOUBLE y[], GINT *nxy, GINT *isz);
#else
void DDOTG_(GDOUBLE *dot, GDOUBLE x[], GDOUBLE y[], GINT *nxy, GINT *isz);
#endif

#if defined(_LINUX) 
void dcopy_(GDOUBLE x[], GDOUBLE y[], GINT *nxy, GINT *isz);
#elif defined(_AIX)
void dcopy (GDOUBLE x[], GDOUBLE y[], GINT *nxy, GINT *isz);
#elif defined(_MACOSX)
void dcopy_(GDOUBLE x[], GDOUBLE y[], GINT *nxy, GINT *isz);
#else
void DCOPY_(GDOUBLE x[], GDOUBLE y[], GINT *nxy, GINT *isz);
#endif

#if defined(_LINUX) 
void isassign_(GDOUBLE u[], GDOUBLE g[], GINT ig[], GINT *nm, GINT *nop);
#elif defined(_AIX)
void isassign (GDOUBLE u[], GDOUBLE g[], GINT ig[], GINT *nm, GINT *nop);
#elif defined(_MACOSX)
void isassign_(GDOUBLE u[], GDOUBLE g[], GINT ig[], GINT *nm, GINT *nop);
#else
void ISASSIGN_(GDOUBLE u[], GDOUBLE g[], GINT ig[], GINT *nm, GINT *nop);
#endif

#if defined(_LINUX) 
void isum_(GDOUBLE *res, GDOUBLE u[], GINT ig[], GINT *ne);
#elif defined(_AIX)
void isum (GDOUBLE *res, GDOUBLE u[], GINT ig[], GINT *ne);
#elif defined(_MACOSX)
void isum_(GDOUBLE *res, GDOUBLE u[], GINT ig[], GINT *ne);
#else
void ISUM_(GDOUBLE *res, GDOUBLE u[], GINT ig[], GINT *ne);
#endif

#if defined(_LINUX) 
void iprod_(GDOUBLE *res, GDOUBLE u[], GINT ig[], GINT *ne);
#elif defined(_AIX)
void iprod (GDOUBLE *res, GDOUBLE u[], GINT ig[], GINT *ne);
#elif defined(_MACOSX)
void iprod_(GDOUBLE *res, GDOUBLE u[], GINT ig[], GINT *ne);
#else
void IPROD__(GDOUBLE *res, GDOUBLE u[], GINT ig[], GINT *ne);
#endif

#if defined(_LINUX) 
void imax_(GDOUBLE *res, GDOUBLE u[], GINT ig[], GINT *ne);
#elif defined(_AIX)
void imax (GDOUBLE *res, GDOUBLE u[], GINT ig[], GINT *ne);
#elif defined(_MACOSX)
void imax_(GDOUBLE *res, GDOUBLE u[], GINT ig[], GINT *ne);
#else
void IMAX_(GDOUBLE *res, GDOUBLE u[], GINT ig[], GINT *ne);
#endif

#if defined(_LINUX) 
void imin_(GDOUBLE *res, GDOUBLE u[], GINT ig[], GINT *ne);
#elif defined(_AIX)
void imin (GDOUBLE *res, GDOUBLE u[], GINT ig[], GINT *ne);
#elif defined(_MACOSX)
void imin_(GDOUBLE *res, GDOUBLE u[], GINT ig[], GINT *ne);
#else
void IMIN_(GDOUBLE *res, GDOUBLE u[], GINT ig[], GINT *ne);
#endif

}

